"""Unconfined RC fiber section; positive concrete fibers exclude rebar circles.
Geometry and torsion are explicit. Concrete02 softening is not regularized here:
this supplies constituents/section, not a mesh-independent collapse model.
"""
import math
import numpy as np
from .materials import register_uniaxial,material,_tag

def rectangle_fibers(width_m,depth_m,bars,*,ny=36,nz=8,quadrature_order=16):
 """bars=[(y_m,z_m,diameter_m),...]; localy is depth, localz is width.
 Each horizontal Gauss slice excludes exact circular chord intervals.
 Breakpoints include bar extents. Areas positive; quadrature error reported.
 """
 if not all(math.isfinite(v) and v>0 for v in [width_m,depth_m]):raise ValueError('positive section dimensions')
 if any(type(n) is not int or n<2 for n in [ny,nz,quadrature_order]):raise ValueError('integer subdivisions/order>=2')
 for i,(y,z,d) in enumerate(bars):
  if not all(math.isfinite(v) for v in [y,z,d]) or d<=0 or abs(y)+d/2>=depth_m/2 or abs(z)+d/2>=width_m/2:raise ValueError('bar must lie strictly inside section')
  for y2,z2,d2 in bars[:i]:
   if math.hypot(y-y2,z-z2)<(d+d2)/2:raise ValueError('overlapping rebar')
 breaks=sorted(set(np.linspace(-depth_m/2,depth_m/2,ny+1).tolist()+[y+sign*d/2 for y,z,d in bars for sign in [-1,1]]))
 fib=[]
 for lo,hi in zip(breaks[:-1],breaks[1:]):
  order=quadrature_order if any(lo<by+d/2 and hi>by-d/2 for by,bz,d in bars) else 2
  gy,gw=np.polynomial.legendre.leggauss(order)
  for q,w in zip(gy,gw):
   y=(lo+hi)/2+(hi-lo)*q/2;dy=(hi-lo)*w/2;holes=[]
   for by,bz,d in bars:
    if abs(y-by)<d/2:
     dz=math.sqrt((d/2)**2-(y-by)**2);holes.append((bz-dz,bz+dz))
   intervals=[];left=-width_m/2
   for start,end in sorted(holes):
    if start>left:intervals.append((left,start))
    left=max(left,end)
   if left<width_m/2:intervals.append((left,width_m/2))
   for start,end in intervals:
    count=max(1,int(math.ceil((end-start)/(width_m/nz))))
    for k in range(count):
     zc=start+(end-start)*(k+.5)/count;half=(end-start)/(2*count)
     for sign in [-1,1]:
      z=zc+sign*half/math.sqrt(3);area=dy*half
      fib.append((float(y),float(z),float(area)))
 area_steel=sum(math.pi*d*d/4 for y,z,d in bars);target=width_m*depth_m-area_steel;area=sum(a for y,z,a in fib)
 error=abs(area-target)/(width_m*depth_m)
 if error>1e-6:raise ValueError(f'concrete quadrature area error {error}; refine ny/order')
 return fib,dict(gross_area_m2=width_m*depth_m,concrete_area_m2=area,steel_area_m2=area_steel,concrete_area_expected_m2=target,area_relative_error=error,concrete_fiber_count=len(fib),concrete_Iz_m4=sum(a*y*y for y,z,a in fib),concrete_Iy_m4=sum(a*z*z for y,z,a in fib),scope='positive concrete quadrature excludes bars; steel lumped at bar centroid, bar self-inertia neglected')

def register_rc_section(ops,section_tag,concrete_tag,rebar_tag,*,width_m,depth_m,bars,ny=36,nz=8,GJ_N_m2=None,elastic=False):
 for t in [section_tag,concrete_tag,rebar_tag]:_tag(t)
 if concrete_tag==rebar_tag:raise ValueError('distinct constituent tags required')
 if GJ_N_m2 is not None and (not math.isfinite(GJ_N_m2) or GJ_N_m2<=0):raise ValueError('positive explicit torsional rigidity')
 fibers,meta=rectangle_fibers(width_m,depth_m,bars,ny=ny,nz=nz)
 register_uniaxial(ops,'concrete_c30_prototype',concrete_tag,elastic=elastic);register_uniaxial(ops,'rebar_sd400_prototype',rebar_tag,elastic=elastic)
 if GJ_N_m2 is None:ops.section('Fiber',section_tag)
 else:ops.section('Fiber',section_tag,'-GJ',GJ_N_m2)
 for y,z,a in fibers:ops.fiber(y,z,a,concrete_tag)
 for y,z,d in bars:ops.fiber(y,z,math.pi*d*d/4,rebar_tag)
 mc=material('concrete_c30_prototype');ms=material('rebar_sd400_prototype')
 return dict(section_tag=section_tag,**meta,mass_per_length_kg_m=meta['concrete_area_m2']*mc['density_kg_m3']+meta['steel_area_m2']*ms['density_kg_m3'],GJ_N_m2=GJ_N_m2,mass_in_material=False,elastic=elastic,
             status='UNCONFINED_NOMINAL_C30_SD400_PROTOTYPE',restriction='nonlinear softening remains sensitive to element/integration choice; use validated fixed formulation scope or supply regularization/hinge-length evidence; 3D requires GJ and shear model as appropriate')

def register_rc_hinge_beam(ops,element_tag,node_i,node_j,transformation_tag,integration_tag,section_tag,elastic_section_tag,concrete_tag,rebar_tag,*,width_m,depth_m,bars,hinge_i_m,hinge_j_m,ny=36,nz=8,GJ_N_m2=None,include_element_mass=True):
 """Fixed HingeRadau physical formulation; hinge lengths are REQUIRED inputs.
 Both end sections share nominal unconfined constituent laws. Interior uses
 the transformed elastic section. No claim of measured hinge length/fracture
 energy; verify axial ratio/drift and use hinge sensitivity for new structure.
 """
 for t in [element_tag,node_i,node_j,transformation_tag,integration_tag,section_tag,elastic_section_tag,concrete_tag,rebar_tag]:_tag(t)
 if section_tag==elastic_section_tag:raise ValueError('distinct section tags')
 coords_i=np.array(ops.nodeCoord(node_i));coords_j=np.array(ops.nodeCoord(node_j));L=float(np.linalg.norm(coords_j-coords_i))
 if not math.isfinite(L) or L<=0:raise ValueError('positive member length')
 if not all(math.isfinite(v) and v>0 for v in [hinge_i_m,hinge_j_m]):raise ValueError('positive explicit hinge lengths')
 if 4*(hinge_i_m+hinge_j_m)>=L:raise ValueError('this wrapper requires a positive elastic interior:4*(lpI+lpJ)<L')
 ndim=len(coords_i)
 if ndim==3 and GJ_N_m2 is None:raise ValueError('3D beam requires explicit torsional rigidity GJ')
 meta=register_rc_section(ops,section_tag,concrete_tag,rebar_tag,width_m=width_m,depth_m=depth_m,bars=bars,ny=ny,nz=nz,GJ_N_m2=GJ_N_m2)
 Ec=30e9;Es=200e9;Gc=Ec/(2*1.2)
 EA=Ec*meta['concrete_area_m2']+Es*meta['steel_area_m2']
 EIz=Ec*meta['concrete_Iz_m4']+Es*sum(math.pi*d*d/4*y*y for y,z,d in bars)
 EIy=Ec*meta['concrete_Iy_m4']+Es*sum(math.pi*d*d/4*z*z for y,z,d in bars)
 if ndim==2:ops.section('Elastic',elastic_section_tag,Ec,EA/Ec,EIz/Ec)
 elif ndim==3:ops.section('Elastic',elastic_section_tag,Ec,EA/Ec,EIz/Ec,EIy/Ec,Gc,GJ_N_m2/Gc)
 else:raise ValueError('2D/3D beam only')
 ops.beamIntegration('HingeRadau',integration_tag,section_tag,hinge_i_m,section_tag,hinge_j_m,elastic_section_tag)
 args=['forceBeamColumn',element_tag,node_i,node_j,transformation_tag,integration_tag,'-iter',40,1e-12]
 if include_element_mass:args+=['-mass',meta['mass_per_length_kg_m']]
 ops.element(*args)
 return dict(meta,element_tag=element_tag,length_m=L,hinge_i_m=hinge_i_m,hinge_j_m=hinge_j_m,integration='HingeRadau',EA_N=EA,EIz_N_m2=EIz,EIy_N_m2=EIy,mass_in_element=include_element_mass,scope='fixed hinge formulation; LP calibration and structure-level stability remain geometry dependent; no arbitrary Lobatto-point substitution')
