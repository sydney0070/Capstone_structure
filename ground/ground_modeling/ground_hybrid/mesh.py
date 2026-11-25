"""3D bbarBrick ground, consistent gravity, surface quadrature and port metadata."""
from collections import defaultdict
import math
import numpy as np
import openseespy.opensees as ops
from ground_modeling import (BRICK_LOCAL_COORDS, brick_shape_derivatives,
                             brick_geometry_metrics, orient_terrain_grid_for_opensees)
from .config import validate_config, graded_axis, depth_levels, assigned_material
from .materials import register

GAUSS = (-1/math.sqrt(3), 1/math.sqrt(3))
FACE_INDEX = {'bottom':(0,3,2,1),'top':(4,5,6,7),'xmin':(0,4,7,3),
              'xmax':(1,2,6,5),'ymin':(0,1,5,4),'ymax':(3,7,6,2)}


def brick_weights(coords):
    weights = np.zeros(8)
    for xi in GAUSS:
        for eta in GAUSS:
            for zeta in GAUSS:
                point=np.array([xi,eta,zeta]);N=np.prod(1+BRICK_LOCAL_COORDS*point,axis=1)/8
                det=np.linalg.det(brick_shape_derivatives(xi,eta,zeta).T @ coords)
                if det <= 0: raise ValueError('inverted element')
                weights += N*det
    return weights


def face_weights(coords):
    signs=np.array([[-1,-1],[1,-1],[1,1],[-1,1]])
    w=np.zeros(4);normal=np.zeros(3)
    for r in GAUSS:
        for s in GAUSS:
            N=np.prod(1+signs*np.array([r,s]),axis=1)/4
            dr=signs[:,0]*(1+signs[:,1]*s)/4
            ds=signs[:,1]*(1+signs[:,0]*r)/4
            cross=np.cross(dr@coords,ds@coords)
            w+=N*np.linalg.norm(cross);normal+=cross
    return w,normal/np.linalg.norm(normal)


def build_ground(config, terrain=None, *, roller_sides=True):
    """Append in caller-owned domain; NEVER wipe. terrain=(X,Y,Z) local meters.
    Rectilinear XY, surface-following layers and planar bottom. Caller must
    reserve tag_start namespaces; tags collide loudly. Gravity is separate.
    Set roller_sides=False when the caller supplies the lateral constraints.
    """
    c=validate_config(config)
    if terrain is None:
        x=graded_axis(c['width_m'],c['core_width_m'],c['core_spacing_m'],c['max_outer_spacing_m'])
        y=graded_axis(c['length_m'],c['core_length_m'],c['core_spacing_m'],c['max_outer_spacing_m'])
        X,Y=np.meshgrid(x,y);Z=np.zeros_like(X)
    else:
        X,Y,Z=[np.array(v,dtype=float,copy=True) for v in terrain]
    if X.shape != Y.shape or X.shape != Z.shape or X.ndim != 2 or min(X.shape)<2 or not np.all(np.isfinite([X,Y,Z])):
        raise ValueError('finite, matching 2D terrain grids with at least2x2 nodes required; no silent void filling')
    origin=(float(X.min()),float(Y.min()))
    X,Y,Z,orientation=orient_terrain_grid_for_opensees(X,Y,Z)
    # Legacy helper translates the origin. Preserve caller's local coordinates
    # so material lenses, footprint and structure ports stay in the same frame.
    X+=origin[0];Y+=origin[1]
    if not np.allclose(X,X[0:1,:]) or not np.allclose(Y,Y[:,0:1]):
        raise ValueError('rectilinear XY terrain required')
    if not np.all(np.diff(X[0])>0) or not np.all(np.diff(Y[:,0])>0):raise ValueError('increasing XY axes required')
    ny,nx=Z.shape;depths=depth_levels(c);nz=len(depths)-1
    if (nx-1)*(ny-1)*nz>c['max_elements']:raise ValueError('element budget exceeded; adjust mesh explicitly')
    start=c['tag_start'];node_count=nx*ny*(nz+1);element_count=(nx-1)*(ny-1)*nz
    if set(range(start,start+node_count)) & set(ops.getNodeTags()):raise ValueError('ground node tag collision')
    if set(range(start,start+element_count)) & set(ops.getEleTags()):raise ValueError('ground element tag collision')
    names=sorted({l['material'] for l in c['layers']}|{l['material'] for l in c['lenses']})
    ops.model('basic','-ndm',3,'-ndf',3)
    mats={n:register(ops,n,start+i,elastic=c['elastic'],stiffness_factor=c['stiffness_factor']) for i,n in enumerate(names)}
    ops.model('basic','-ndm',3,'-ndf',3)
    nodes={};index={};loads=defaultdict(lambda:np.zeros(3));base_z=float(np.min(Z)-c['depth_m'])
    for k,d in enumerate(depths):
        for j in range(ny):
            for i in range(nx):
                tag=start+len(nodes);xyz=np.array([X[j,i],Y[j,i],base_z if k==nz else Z[j,i]-d])
                ops.node(tag,*xyz);nodes[tag]=xyz;index[i,j,k]=tag
    elements={};faces=[];fractions=defaultdict(float);mass=0.;max_aspect=0.
    for k in range(nz):
        for j in range(ny-1):
            for i in range(nx-1):
                keys=[(i,j,k+1),(i+1,j,k+1),(i+1,j+1,k+1),(i,j+1,k+1),
                      (i,j,k),(i+1,j,k),(i+1,j+1,k),(i,j+1,k)]
                tags=[index[q] for q in keys];coords=np.array([nodes[t] for t in tags]);geom=brick_geometry_metrics(coords)
                if geom['min_jacobian_det']<=0 or geom['aspect_ratio']>c['max_aspect_ratio']:
                    raise ValueError(f'invalid/distorted brick at {i,j,k}: {geom}')
                center=coords.mean(axis=0);surface=Z[j:j+2,i:i+2].mean();depth=float(surface-center[2])
                name=assigned_material(c,*center[:2],depth);m=mats[name];tag=start+len(elements)
                ops.element('bbarBrick',tag,*tags,m['tag'],0.,0.,0.)
                volumes=brick_weights(coords);volume=float(sum(volumes));rho=m['density_kg_m3']
                for t,v in zip(tags,volumes):loads[t][2]-=rho*c['gravity_m_s2']*v
                mass+=volume*rho;fractions[name]+=volume;max_aspect=max(max_aspect,geom['aspect_ratio'])
                elements[tag]=dict(nodes=tags,material=name,volume_m3=volume,center_m=center.tolist(),depth_m=depth,
                                   height_m=float(coords[4:,2].mean()-coords[:4,2].mean()))
                boundary=[]
                if k==0:boundary.append('top')
                if k==nz-1:boundary.append('bottom')
                if i==0:boundary.append('xmin')
                if i==nx-2:boundary.append('xmax')
                if j==0:boundary.append('ymin')
                if j==ny-2:boundary.append('ymax')
                for side in boundary:
                    f=[tags[q] for q in FACE_INDEX[side]];w,n=face_weights(np.array([nodes[t] for t in f]))
                    faces.append(dict(side=side,nodes=f,weights_m2=w.tolist(),area_m2=float(sum(w)),
                                      normal=n.tolist(),element=tag,material=name))
    fixities={t:[0,0,0] for t in nodes}
    for (i,j,k),t in index.items():
        if roller_sides:
            if i in (0,nx-1):fixities[t][0]=1
            if j in (0,ny-1):fixities[t][1]=1
        if k==nz:fixities[t]=[1,1,1]
        if any(fixities[t]):ops.fix(t,*fixities[t])
    total=sum(loads.values());volume=sum(fractions.values())
    if abs(total[2]+mass*c['gravity_m_s2'])>1e-8*max(abs(total[2]),1):raise RuntimeError('gravity assembly balance')
    return dict(config=c,nodes=nodes,index=index,elements=elements,faces=faces,materials=mats,
                gravity_loads=dict(loads),fixities=fixities,mass_kg=mass,volume_m3=volume,
                material_volume_fractions={n:v/volume for n,v in fractions.items()},
                mesh=dict(nx=nx,ny=ny,nz=nz,max_aspect_ratio=max_aspect,depth_levels_m=depths,
                          orientation=orientation),
                boundary='static_box' if roller_sides else 'fixed_base',stage='built',
                displacement_reference={},attachments=[],boundary_force_loads={})
