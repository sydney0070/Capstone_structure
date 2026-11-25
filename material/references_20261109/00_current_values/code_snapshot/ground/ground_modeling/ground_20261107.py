"""Dated ground/input/kinematics corrections for structural SSI research."""
from pathlib import Path
from collections import defaultdict
import hashlib,json,math,re,os,tempfile
os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"matplotlib"))
import numpy as np
import openseespy.opensees as ops
from insight_motion import preprocess_acc_mseed,load_acceleration_csv


def build_ground_20261107(config,terrain=None):
    """Audited append-only builder; retain established brick/mass implementation."""
    from ground_hybrid.mesh import build_ground
    return validate_ground_20261107(build_ground(config,terrain))


def validate_ground_20261107(g):
    """Apply the same frozen-material audit to Cartesian or conforming builders."""
    import sys
    project=Path(__file__).resolve().parents[2]
    sys.path.insert(0,str(project/'material'))
    from materials_20261107 import validate_record
    for name,record in g['materials'].items():validate_record(name,record)
    g['revision_validation']=dict(date='20261107',materials_checked=list(g['materials']),
        max_element_aspect_ratio=g['mesh']['max_aspect_ratio'],mass_kg=g['mass_kg'],
        site_claim=g['config'].get('site_scope','Synthetic reference unless caller supplies independently registered terrain/geology'),
        constitutive_scope='DP pressure-independent elastic stiffness + frictional plasticity; no G/Gmax/damping curve refit')
    return g


def file_sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def require_finite(values, context):
    """Reject raw nonfinite response values before reductions can mask NaNs."""
    a=np.asarray(values,dtype=float)
    if not np.isfinite(a).all():
        raise RuntimeError(f'Nonfinite response: {context}')
    return a


def insight_reference_site(width_m, nx=31, ny=31, root=None):
    """Anchor an idealized horizontal site to ELYSE; preserve the real DEM samples.

    No cut/fill or construction history is inferred. The flat plane is a geometry
    idealization, not the measured surface or a reconstructed soil profile.
    """
    from insight_site_dem import read_insight_hirise_patch
    root=Path(root) if root is not None else Path(__file__).resolve().parent
    img=root/'DTEEC_037262_1845_036761_1845_U01.IMG'
    xml=root/'stationxml.xb.elyse.2019.151.xml'
    site=station_georeference(img,xml)
    X,Y,_,z=read_insight_hirise_patch(img,width_m=width_m,height_m=width_m,nx=nx,ny=ny,
        center_row=site['row_0based'],center_col=site['col_0based'])
    _,_,_,center=read_insight_hirise_patch(img,width_m=2.,height_m=2.,nx=3,ny=3,
        center_row=site['row_0based'],center_col=site['col_0based'])
    elevation=float(center[1,1])
    import xml.etree.ElementTree as ET
    ns={'s':'http://www.fdsn.org/xml/station/1'}
    station=next(s for s in ET.parse(xml).getroot().findall('.//s:Station',ns) if s.get('code')=='ELYSE')
    station_elevation=float(station.find('s:Elevation',ns).text)
    site.update(kind='insight_station_centered_flat_reference',width_m=width_m,
        local_origin='ELYSE StationXML coordinates; x East, y North, z Up',
        local_z0_dem_elevation_m=elevation,station_metadata_elevation_m=station_elevation,
        dem_minus_station_elevation_m=elevation-station_elevation,
        dem_elevation_range_m=[float(z.min()),float(z.max())],
        max_flattening_difference_m=float(np.max(abs(z-elevation))),
        terrain_mode='horizontal reference at station DEM elevation; actual relief sampled but not applied',
        geology_scope='Unchanged assumed InSight crust proxy + JEZ-1 + basalt; not measured ELYSE stratigraphy',
        source_files={str(p):file_sha256(p) for p in (img,xml)},
        coordinate_caveat='Station metadata anchor near InSight lander, not an independent lander-to-SEIS survey; DEM vertical disagreement retained')
    return site,np.column_stack((X.ravel(),Y.ravel(),z.ravel(),(z-elevation).ravel(),np.zeros(z.size)))


def station_georeference(img_path,stationxml):
    """PDS3 zero-based pixel centres: x=(sample-offset)*scale, y=(offset-line)*scale.
    Exact frame agreement of DEM/StationXML is not independently surveyed here.
    """
    import xml.etree.ElementTree as ET
    with Path(img_path).open('rb') as f:label=f.read(65536).decode('ascii','replace')
    def number(key):
        m=re.search(r'^\s*'+re.escape(key)+r'\s*=\s*([-+0-9.eE]+)',label,re.M)
        if not m:raise ValueError(f'Missing projection metadata: {key}')
        return float(m.group(1))
    if not re.search(r'MAP_PROJECTION_TYPE\s*=\s*"EQUIRECTANGULAR"',label) or number('MAP_PROJECTION_ROTATION')!=0 or number('CENTER_LATITUDE')!=0:raise ValueError('Unsupported DEM orientation')
    ns={'s':'http://www.fdsn.org/xml/station/1'};root=ET.parse(stationxml).getroot()
    st=next(s for s in root.findall('.//s:Station',ns) if s.get('code')=='ELYSE')
    lat=float(st.find('s:Latitude',ns).text);lon=float(st.find('s:Longitude',ns).text)
    resolution=number('MAP_RESOLUTION');row=number('LINE_PROJECTION_OFFSET')-lat*resolution
    col=number('SAMPLE_PROJECTION_OFFSET')+(lon-number('CENTER_LONGITUDE'))*resolution
    rows,cols=int(number('LINES')),int(number('LINE_SAMPLES'))
    if not 0<=row<rows-1 or not 0<=col<cols-1:raise ValueError('Station outside DEM')
    return dict(station='ELYSE',latitude_deg=lat,longitude_deg=lon,row_0based=row,col_0based=col,
        projection='north-up planetocentric equirectangular',pixel_scale_m=number('MAP_SCALE'),
        caveat='Coordinate conversion, not independent DEM/station datum co-registration verification')


def resolution(g,freqmax,dt):
    from ground_hybrid.dynamics import mesh_frequency_guidance
    r=mesh_frequency_guidance(g,freqmax)
    r.update(integration_dt_s=dt,samples_per_input_cycle=1/(dt*freqmax),
        temporal_bandwidth_pass=dt<=1/(20*freqmax)+1e-12,
        scope='Frequency/mesh heuristic using reference G,K; not a time/mesh convergence proof')
    if not r['mesh_target_pass']:raise ValueError('Input bandwidth exceeds 10-elements/wavelength mesh target; reduce --freqmax or refine --spacing')
    if not r['temporal_bandwidth_pass']:raise ValueError('At least 20 integration samples per highest input cycle required')
    return r


def prepare_motion(source,out,*,duration=2.,warmup=5.,start=None,freqmax=1.,rate=10.,allow_short_warmup=False):
    """Process full record, then select warm-up + target; no re-taper of target samples."""
    out=Path(out);full=preprocess_acc_mseed(source,out/'motion_full_enu.csv',freqmax=freqmax,target_rate=rate)
    dt,xyz=load_acceleration_csv(full)
    if not np.isfinite([duration,warmup]).all() or duration<2*dt or warmup<0:raise ValueError('Invalid motion duration/warmup')
    peak=int(np.argmax(np.linalg.norm(xyz,axis=1)))
    count=int(round(duration/dt))
    target=max(0,min(peak-count//2,len(xyz)-count-1)) if start is None else int(round(start/dt))
    if start is not None and (not math.isfinite(start) or start<0):raise ValueError('Invalid source target start')
    end=target+count
    if target<0 or end>=len(xyz):raise ValueError('Requested target interval outside record; no silent truncation')
    requested_first=target-int(round(warmup/dt))
    if requested_first<0 and not allow_short_warmup:
        raise ValueError('Requested warmup precedes record start; explicitly allow_short_warmup to truncate')
    first=max(0,requested_first);selected=xyz[first:end+1].copy()
    warm_n=target-first
    # Start smoothing affects warm-up only. Target interval is byte-for-byte preserved.
    n=min(warm_n,int(round(min(1.,warmup)/dt)))
    if n>=2:selected[:n]*=(.5*(1-np.cos(np.linspace(0,np.pi,n))))[:,None]
    meta=json.loads(full.with_suffix('.json').read_text())
    meta.update(source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),source_peak_s=peak*dt,
        source_start_s=first*dt,target_start_source_s=target*dt,target_end_source_s=end*dt,
        target_start_model_s=warm_n*dt,effective_warmup_s=warm_n*dt,requested_warmup_s=warmup,requested_duration_s=duration,effective_duration_s=count*dt,
        requested_target_start_s=start,warmup_truncated=requested_first<0,allow_short_warmup=allow_short_warmup,
        selection='peak_acceleration_window_NOT_global_max_displacement',target_retapered=False,
        zero_warmup_warning=warm_n==0,units_verification='Supplied ACC product assumed m/s2; original product calibration history still external',
        history_caveat='Truncated preceding history; warmup sensitivity required for nonlinear assessment')
    p=out/'motion_window_enu.csv'
    np.savetxt(p,np.column_stack((np.arange(len(selected))*dt,selected)),delimiter=',',header='time_s,east_m_s2,north_m_s2,up_m_s2',comments='')
    p.with_suffix('.json').write_text(json.dumps(meta,indent=2))
    return dt,selected,meta


def lateral_dashpots(g,tag=1000000):
    """Fixed base retained, static side reactions preserved, viscous side tractions.
    Approximate local radiating sides; not free-field or oblique-wave validation.
    """
    if g['boundary']!='static_box':raise ValueError('Boundary transition allowed once')
    base={t for f in g['faces'] if f['side']=='bottom' for t in f['nodes']}
    coefficients=defaultdict(lambda:np.zeros(3))
    for f in g['faces']:
        if f['side'] not in ('xmin','xmax','ymin','ymax'):continue
        m=g['materials'][f['material']];rho=m['density_kg_m3'];vs=math.sqrt(m['G_Pa']/rho);vp=math.sqrt((m['K_Pa']+4*m['G_Pa']/3)/rho)
        axis=0 if f['side'][0]=='x' else 1
        for t,area in zip(f['nodes'],f['weights_m2']):
            if t in base:continue
            C=np.full(3,rho*vs*area);C[axis]=rho*vp*area;coefficients[t]+=C
    ops.loadConst('-time',0.);ops.reactions()
    pattern=tag+9000;ops.timeSeries('Constant',pattern);ops.pattern('Plain',pattern,pattern)
    released=[]
    for t,fixed in g['fixities'].items():
        if t in base:continue
        r=ops.nodeReaction(t);load=np.zeros(3)
        for i,yes in enumerate(fixed):
            if yes:load[i]=r[i];released.append((t,i+1))
        if any(fixed):ops.load(t,*load)
    for t,dof in released:ops.remove('sp',t,dof)
    ops.model('basic','-ndm',3,'-ndf',3)
    for i,(t,C) in enumerate(sorted(coefficients.items())):
        anchor=tag+i
        if anchor in ops.getNodeTags():raise ValueError('Dashpot tag collision')
        ops.node(anchor,*ops.nodeCoord(t));ops.fix(anchor,1,1,1)
        for d,c in enumerate(C,1):
            e=tag+10000+3*i+d
            ops.uniaxialMaterial('Viscous',e,float(c),1.)
            ops.element('zeroLength',e,t,anchor,'-mat',e,'-dir',d)
    # Commit newly applied balancing loads BEFORE Newmark initializes acceleration.
    # Otherwise initial acceleration can be computed from stale nodal unbalances.
    physical=list(g['nodes']);before=np.array([ops.nodeDisp(t)[:3] for t in physical])
    ops.wipeAnalysis();ops.constraints('Plain');ops.numberer('RCM');ops.system('UmfPack')
    ops.test('NormUnbalance',max(.001,g['mass_kg']*g['config']['gravity_m_s2']*1e-10),80)
    ops.algorithm('Newton');ops.integrator('LoadControl',0.);ops.analysis('Static')
    if ops.analyze(1):raise RuntimeError('Boundary transition equilibrium failed')
    movement=float(np.max(np.linalg.norm(np.array([ops.nodeDisp(t)[:3] for t in physical])-before,axis=1)))
    if movement>1e-8:raise RuntimeError(f'Boundary transition changed static state: {movement}')
    ops.wipeAnalysis()
    g['boundary']='fixed_base_viscous_sides'
    return dict(kind=g['boundary'],released_dofs=len(released),side_nodes=len(coefficients),base_nodes_retained=len(base),transition_movement_m=movement,
        caveat='Reference impedances, fixed base, no free-field boundary; static reactions conserved; rest/attenuation checks required')


def remove_base_rigid_motion(reference_base,current_base,reference_points,current_points):
    """Least-squares finite rigid fit of base translations; removes rotation too."""
    B=np.asarray(reference_base);C=np.asarray(current_base)
    x=B-B.mean(axis=0);y=C-C.mean(axis=0)
    u,s,vt=np.linalg.svd(x.T@y);R=u@vt
    if np.linalg.det(R)<0:u[:,-1]*=-1;R=u@vt
    rigid=(np.asarray(reference_points)-B.mean(axis=0))@R+C.mean(axis=0)
    return np.asarray(current_points)-rigid


def advance(dt,*,max_halvings=6,level=0,attempts=None,on_commit=None):
    """Fail loudly after bounded subdivisions; never skip a failed interval."""
    if attempts is None:attempts=[]
    before=ops.getTime();ok=ops.analyze(1,dt)
    attempts.append(dict(t=before,dt=dt,code=ok))
    if ok:
        if abs(ops.getTime()-before)>1e-9:raise RuntimeError('Failed step changed committed time')
        if level>=max_halvings:raise RuntimeError('Dynamic interval failed after bounded substeps; no physical failure inference')
        advance(dt/2,max_halvings=max_halvings,level=level+1,attempts=attempts,on_commit=on_commit)
        advance(dt/2,max_halvings=max_halvings,level=level+1,attempts=attempts,on_commit=on_commit)
    elif abs(ops.getTime()-before-dt)>1e-8:raise RuntimeError('Unexpected committed step time')
    elif on_commit is not None:on_commit()
    return attempts


def ground_state(g, *, gauss_points=1):
    """Snapshot stress/strain extrema; 1 or all 8 GP, not a capacity verdict."""
    from capmat.materials import dp_constants
    if gauss_points not in (1,8):raise ValueError('gauss_points must be 1 or 8')
    strain_max=stress_max=0.;yield_max=-float('inf');yield_element=yield_ip=None;low_pressure={}
    constants={name:dp_constants(m['phi_deg'],m['c_Pa'],m['psi_deg']) for name,m in g['materials'].items() if m['active_law']=='DruckerPrager'}
    for e,record in g['elements'].items():
        name=record['material'];m=g['materials'][name]
        for ip in range(1,gauss_points+1):
            stress=np.asarray(ops.eleResponse(e,'material',ip,'stress'),float)
            strain=np.asarray(ops.eleResponse(e,'material',ip,'strain'),float)
            if len(stress)!=6 or len(strain)!=6 or not np.isfinite(np.r_[stress,strain]).all():raise RuntimeError('Ground stress/strain response unavailable')
            stress_max=max(stress_max,float(max(abs(stress))));strain_max=max(strain_max,float(max(abs(strain))))
            if m['active_law']=='DruckerPrager':
                p=constants[name];I1=float(sum(stress[:3]));dev=stress[:3]-I1/3
                norm=math.sqrt(float(dev@dev+2*(stress[3:6]@stress[3:6])))
                f=norm+p['rho']*I1-math.sqrt(2/3)*p['sigmaY']
                if f>yield_max:yield_max=f;yield_element=e;yield_ip=ip
            if name.startswith('moon_jsc1a_'):
                # Test evidence is triaxial confinement 15--60 kPa. Mean
                # effective compression is only a screening proxy, not sigma3.
                row=low_pressure.setdefault(name,dict(sampled_GP=0,mean_compression_below_15kPa_GP=0,minimum_mean_compression_Pa=float('inf'),maximum_mean_compression_Pa=-float('inf')))
                mean=-float(sum(stress[:3]))/3;row['sampled_GP']+=1;row['mean_compression_below_15kPa_GP']+=int(mean<15000)
                row['minimum_mean_compression_Pa']=min(row['minimum_mean_compression_Pa'],mean);row['maximum_mean_compression_Pa']=max(row['maximum_mean_compression_Pa'],mean)
    return dict(sampled_elements=len(g['elements']),integration_point=1 if gauss_points==1 else '1..8',sampled_gauss_points=len(g['elements'])*gauss_points,abs_stress_component_max_Pa=stress_max,
        abs_strain_component_max=strain_max,DP_max_yield_function_Pa=yield_max if yield_element is not None else None,
        DP_element_at_max=yield_element,DP_ip_at_max=yield_ip,lunar_low_pressure_screen=low_pressure,
        scope=f'{gauss_points} Gauss point(s) per element at snapshot time. DP shear surface diagnostic only; tensile surface not checked. Mean pressure screen is not triaxial confinement validation. Not time extrema, plastic work, bearing capacity or failure verdict.')


def static_ramp(g,loads,*,pattern=800001,steps=20,max_attempts=150,algorithms=("KrylovNewton","NewtonLineSearch","Newton"),on_commit=None,attempt_observer=None,initial_increment=None,linear_system='UmfPack',first_step_algorithm=None):
    """Bounded DP/elastic gravity or subsequent load ramp; no tolerance relaxation."""
    from ground_hybrid.solver import continue_load,force_balance
    if linear_system not in ('UmfPack','SuperLU'):raise ValueError('validated sparse unsymmetric system required')
    if any(m['active_law']=='PressureDependMultiYield02' for m in g['materials'].values()):raise ValueError('Use separately validated PDMY stage initialization')
    ops.timeSeries('Linear',pattern);ops.pattern('Plain',pattern,pattern)
    for t,v in loads.items():
        ndf=len(ops.nodeDisp(t));v=list(v)
        if len(v)>ndf or len(v)<3 or not np.isfinite(v).all():raise ValueError('Invalid load dimensions/values')
        ops.load(t,*(v+[0.]*(ndf-len(v))))
    ops.wipeAnalysis();ops.constraints('Plain');ops.numberer('RCM');ops.system(linear_system)
    norm=sum(float(np.linalg.norm(v[:3])) for v in loads.values())
    ops.test('NormUnbalance',max(.001,norm*1e-10),80);ops.algorithm('Newton');ops.integrator('LoadControl',0.);ops.analysis('Static')
    count=[0]
    def progress(a):
        count[0]+=1
        print(f"[static {pattern}] factor={a['after']:.4f}, increment={a['increment']:.3g}",flush=True)
        if a['code']==0 and on_commit is not None:on_commit()
        if count[0]>max_attempts:raise RuntimeError('Static load attempt limit exceeded; no strength or tolerance changed')
    attempts=continue_load(pattern,steps=steps,min_increment=1e-5,progress=progress,algorithms=algorithms,attempt_observer=attempt_observer,initial_increment=initial_increment,max_attempts=max_attempts,first_step_algorithm=first_step_algorithm)
    if abs(ops.getLoadFactor(pattern)-1)>1e-8:raise RuntimeError('Incomplete load factor')
    ops.loadConst('-time',0.);ops.integrator('LoadControl',0.)
    if ops.analyze(1):raise RuntimeError('Zero-step final equilibrium failed')
    old=g.setdefault('applied_loads',{})
    for t,v in loads.items():
        a=np.zeros(max(len(old.get(t,[])),len(v)));a[:len(old.get(t,[]))]+=old.get(t,[]);a[:len(v)]+=v;old[t]=a
    balance=force_balance(g)
    applied_M=sum((np.cross(np.array(ops.nodeCoord(t)),np.asarray(v[:3]))+np.pad(np.array(v[3:6]),(0,3-len(v[3:6]))) for t,v in old.items()),np.zeros(3))
    support_F=np.zeros(3);support_M=np.zeros(3)
    for t,fix in g['fixities'].items():
        reaction=np.asarray(ops.nodeReaction(t)[:3])*np.array(fix)
        support_F+=reaction;support_M+=np.cross(np.array(ops.nodeCoord(t)),reaction)
    scale=max(sum(np.linalg.norm(v[:3]) for v in old.values()),1.)
    force_error=float(np.linalg.norm(support_F+np.sum([v[:3] for v in old.values()],axis=0))/scale)
    moment_error=float(np.linalg.norm(support_M+applied_M)/(scale*max(g['config']['width_m'],1.)))
    if force_error>1e-5:raise RuntimeError('Physical ground support force equilibrium failed')
    g['stage']='gravity_equilibrated' if pattern==800001 else 'static_loaded'
    ops.wipeAnalysis()
    return dict(final_factor=float(ops.getLoadFactor(pattern)),attempts=len(attempts),linear_system=linear_system,force_balance=balance,
        physical_support_force_error=force_error,physical_support_moment_error=moment_error,
        moment_caveat='Contact bearing auxiliary rotational reactions excluded; moment residual is diagnostic for contact')


def inspect_site():
    """Read-only site audit; vertical disagreement is reported, never auto-corrected."""
    from insight_site_dem import read_insight_hirise_patch
    import xml.etree.ElementTree as ET
    root=Path(__file__).resolve().parent
    xml=root/'stationxml.xb.elyse.2019.151.xml';img=root/'DTEEC_037262_1845_036761_1845_U01.IMG'
    d=station_georeference(img,xml)
    _,_,_,z=read_insight_hirise_patch(img,width_m=60,height_m=60,nx=5,ny=5,center_row=d['row_0based'],center_col=d['col_0based'])
    ns={'s':'http://www.fdsn.org/xml/station/1'};st=next(s for s in ET.parse(xml).getroot().findall('.//s:Station',ns) if s.get('code')=='ELYSE')
    elevation=float(st.find('s:Elevation',ns).text)
    d.update(station_metadata_elevation_m=elevation,dem_at_station_m=float(z[2,2]),dem_minus_station_elevation_m=float(z[2,2]-elevation),
             dem_patch_range_m=[float(z.min()),float(z.max())],
             caveat='Station-centred pixel conversion works; horizontal/vertical datum or DEM accuracy still needs external verification. Synthetic coupled model does not use this DEM.')
    return d

if __name__=='__main__':
    import argparse,sys
    if '--inspect-site' in sys.argv:
        p=argparse.ArgumentParser(description='Ground-only site audit')
        p.add_argument('--inspect-site',action='store_true');p.add_argument('--output',type=Path)
        a=p.parse_args()
        d=inspect_site()
        if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2))
        print(json.dumps(d,indent=2))
    else:
        from ground_only_20261109 import supervised_cli
        raise SystemExit(supervised_cli())
