"""Native HGT/PDS sampling without GDAL. Only the supplied known formats.
X east, Y north in a small spherical local tangent approximation. Equirectangular
longitude spacing is corrected by cos(latitude), essential near the lunar pole.
"""
import hashlib
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from .config import graded_axis

ROOT=Path(__file__).resolve().parents[1]


def _bilinear(mm,rows,cols,nodata=None):
    i=np.floor(rows).astype(int);j=np.floor(cols).astype(int)
    if i.min()<0 or j.min()<0 or i.max()+1>=mm.shape[0] or j.max()+1>=mm.shape[1]:
        raise ValueError('requested local patch/stencil outside DEM tile; no extrapolation')
    a=rows-i;b=cols-j
    values=np.array([mm[i,j],mm[i+1,j],mm[i,j+1],mm[i+1,j+1]],dtype=float)
    if not np.all(np.isfinite(values)) or (nodata is not None and np.any(values==nodata)):
        raise ValueError('DEM nodata in sampling stencil; no silent void filling')
    z=(1-a)*(1-b)*values[0]+a*(1-b)*values[1]+(1-a)*b*values[2]+a*b*values[3]
    digest=hashlib.sha256(np.stack([rows,cols,*values]).tobytes()).hexdigest()
    return z,dict(row_min=int(i.min()),row_max=int(i.max()+1),col_min=int(j.min()),col_max=int(j.max()+1),sample_sha256=digest)


def native_dem(config,*,prepare_pad=False):
    c=config;body=c['body']
    x=graded_axis(c['width_m'],c['core_width_m'],c['core_spacing_m'],c['max_outer_spacing_m'])
    y=graded_axis(c['length_m'],c['core_length_m'],c['core_spacing_m'],c['max_outer_spacing_m'])
    X,Y=np.meshgrid(x,y)
    if max(c['width_m'],c['length_m'])>5000:raise ValueError('local tangent reader limited to5km domain; use projection for larger sites')
    scale=1.;offset=0.;nodata=None;byte_offset=0
    if body=='Earth':
        image=ROOT/'n36e129.hgt';sources=[image]
        n=int(math.isqrt(image.stat().st_size//2))
        if image.stat().st_size!=2*n*n:raise ValueError('invalid HGT size')
        shape=(n,n);dtype='>i2';north=37.;west=129.;ppd=n-1;radius=6371000.;pixel_center=False;nodata=-32768
        lat0=c.get('site_latitude_deg',36.11);lon0=c.get('site_longitude_deg',129.36)
    elif body=='Moon':
        label=ROOT/'ldem_1024_90s_75s_000_030.xml';tree=ET.parse(label).getroot()
        def text(name,root=tree):
            found=root.find('.//{*}'+name)
            if found is None:raise ValueError(f'missing PDS4 {name}')
            return found.text.strip()
        if text('map_projection_name')!='Equirectangular' or float(text('standard_parallel_1'))!=0:
            raise ValueError('unsupported lunar projection')
        arr=tree.find('.//{*}Array_2D_Image')
        if text('data_type',arr)!='SignedLSB2' or text('axis_index_order',arr)!='Last Index Fastest':raise ValueError('unsupported PDS4 type/order')
        axes={text('axis_name',axis):int(text('elements',axis)) for axis in arr.findall('{*}Axis_Array')}
        shape=(axes['Line'],axes['Sample']);dtype='<i2';byte_offset=int(arr.find('{*}offset').text)
        image=ROOT/text('file_name');sources=[label,image]
        scale=float(text('scaling_factor',arr));offset=float(text('value_offset',arr))
        north=float(text('north_bounding_coordinate'));west=float(text('west_bounding_coordinate'))
        ppd=float(text('pixel_scale_x'));ppy=float(text('pixel_scale_y'))
        if ppd!=ppy:raise ValueError('unsupported nonsquare angular pixels')
        radius=float(text('a_axis_radius'))*1000;pixel_center=True
        lat0=c.get('site_latitude_deg',.5*(north+float(text('south_bounding_coordinate'))))
        lon0=c.get('site_longitude_deg',.5*(west+float(text('east_bounding_coordinate'))))
    elif body=='Mars':
        label=ROOT/'megt00n000hb.lbl';labeltext=label.read_text()
        def value(key):
            m=re.search(r'^\s*'+re.escape(key)+r'\s*=\s*([^\r\n]+)',labeltext,re.M)
            if not m:raise ValueError(f'missing PDS3 {key}')
            return m.group(1).split('<')[0].strip().strip('"')
        if value('SAMPLE_TYPE')!='MSB_INTEGER' or int(value('SAMPLE_BITS'))!=16 or value('MAP_PROJECTION_TYPE')!='SIMPLE CYLINDRICAL':
            raise ValueError('unsupported Mars PDS3 encoding/projection')
        image=ROOT/value('^IMAGE').lower();sources=[label,image]
        shape=(int(value('LINES')),int(value('LINE_SAMPLES')));dtype='>i2'
        north=float(value('MAXIMUM_LATITUDE'));west=float(value('WESTERNMOST_LONGITUDE'));ppd=float(value('MAP_RESOLUTION'))
        radius=float(value('A_AXIS_RADIUS'))*1000;pixel_center=True
        lat0=c.get('site_latitude_deg',.5*(north+float(value('MINIMUM_LATITUDE'))))
        lon0=c.get('site_longitude_deg',.5*(west+float(value('EASTERNMOST_LONGITUDE'))))
    else:raise ValueError('unknown body')
    if not np.all(np.isfinite([lat0,lon0])) or abs(lat0)>=89.9:raise ValueError('site latitude too close to pole for local tangent mapping')
    coslat=math.cos(math.radians(lat0));dlat=Y/radius*180/math.pi;dlon=X/(radius*coslat)*180/math.pi
    shift=.5 if pixel_center else 0.
    rows=(north-(lat0+dlat))*ppd-shift;cols=(lon0+dlon-west)*ppd-shift
    needed=byte_offset+2*math.prod(shape)
    if image.stat().st_size<needed:raise ValueError('truncated raster file')
    mm=np.memmap(image,mode='r',dtype=dtype,shape=shape,offset=byte_offset)
    raw,window=_bilinear(mm,rows,cols,nodata);del mm
    elevation=raw*scale;center=float(elevation[len(y)//2,len(x)//2]);Z=elevation-center;original=Z.copy()
    if prepare_pad:
        hx=c['core_width_m']/2;hy=c['core_length_m']/2
        within=(abs(X)<=hx+1e-9)&(abs(Y)<=hy+1e-9)
        pad=float(np.min(Z[within]));distance=np.maximum(np.maximum(abs(X)-hx,abs(Y)-hy),0.)
        w=np.clip(1-distance/max(c['core_spacing_m']*2,1.),0,1);w=w*w*(3-2*w)
        Z=(1-w)*Z+w*pad
    dy=radius*math.pi/(180*ppd);dx=dy*coslat
    return (X,Y,Z),dict(kind='native_DEM',body=body,site_latitude_deg=lat0,site_longitude_deg=lon0,
        native_spacing_m=[dx,dy],FE_min_spacing_m=[float(min(np.diff(x))),float(min(np.diff(y)))],
        dtype=dtype,shape=shape,scaling_factor=scale,radius_offset_m=offset,center_elevation_relative_reference_m=center,
        coordinate_convention='X east Y north, local spherical tangent; longitude cos(latitude) correction',
        sources=[dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in sources],window=window,
        prepare_pad=prepare_pad,max_grade_change_m=float(np.max(abs(Z-original))),
        scope='native pixel interpolation does not add terrain detail; assumed subsurface; initial grading only, no excavation history')
