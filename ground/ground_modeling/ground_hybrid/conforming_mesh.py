"""Square far field with a fixed polygonal disk/annulus footprint.

An angular mesh is joined to a tensor-product square core and a square exterior.
Radial refinement subdivides the SAME polygon edges; angular refinement changes
the circle approximation and must be studied separately. Static box only.
"""
from collections import defaultdict
import math
import numpy as np
import openseespy.opensees as ops
from ground_modeling import brick_geometry_metrics
from .config import validate_config,depth_levels,assigned_material
from .materials import register
from .mesh import brick_weights,face_weights,FACE_INDEX


def footprint_quads(radius, width, domain_width, *, kind='ring', sectors=24, radial_subdivisions=1):
    if kind not in ('ring','raft'):raise ValueError('ring/raft required')
    if type(sectors) is not int or sectors<24 or sectors%24:raise ValueError('sectors must be a positive multiple of 24')
    if type(radial_subdivisions) is not int or radial_subdivisions<1:raise ValueError('positive integer radial refinement required')
    if not all(math.isfinite(v) and v>0 for v in (radius,width,domain_width)) or width>=2*radius:raise ValueError('positive finite geometry with width<diameter')
    inner=radius-width/2;outer=radius+width/2;half=domain_width/2
    if half<=outer:raise ValueError('far boundary must enclose footprint')
    a=min(radius/4,inner/2);theta=np.arange(sectors)*2*math.pi/sectors
    direction=np.c_[np.cos(theta),np.sin(theta)]
    core=direction*(a/np.max(abs(direction),axis=1))[:,None]
    far=direction*(half/np.max(abs(direction),axis=1))[:,None]
    # No cells cross the two exact footprint radii. Additional rings grade the
    # transition, without moving the inner/outer polygon under radial refinement.
    radii=sorted(set([a*1.5,inner,radius,outer,(outer+half)/2]))
    rings=[(core,None)]+[(direction*r,r) for r in radii]+[(far,None)]
    xy=[];lookup={};cells=[];selected=[]
    def node(point):
        key=tuple(np.round(point,11))
        if key not in lookup:lookup[key]=len(xy);xy.append(np.asarray(point))
        return lookup[key]
    # Every polar core-boundary point belongs to this nonuniform square grid.
    axis=np.unique(np.round(np.r_[core[np.isclose(abs(core[:,1]),a),0],-a,a],11))
    grid={(i,j):node((x,y)) for j,y in enumerate(axis) for i,x in enumerate(axis)}
    for j in range(len(axis)-1):
        for i in range(len(axis)-1):
            cells.append([grid[i,j],grid[i+1,j],grid[i+1,j+1],grid[i,j+1]]);selected.append(kind=='raft')
    previous=[node(p) for p in core]
    for index in range(len(rings)-1):
        p0,r0=rings[index];p1,r1=rings[index+1]
        is_selected=(r1 is not None and r1<=outer+1e-10) if kind=='raft' else (r0 is not None and r1 is not None and r0>=inner-1e-10 and r1<=outer+1e-10)
        for sub in range(1,radial_subdivisions+1):
            current=[node(p) for p in p0+(p1-p0)*(sub/radial_subdivisions)]
            for i in range(sectors):
                j=(i+1)%sectors
                cells.append([previous[i],current[i],current[j],previous[j]]);selected.append(is_selected)
            previous=current
    xy=np.array(xy);cells=np.array(cells,dtype=int)
    # All edges have either two adjacent quads or lie on the square far field.
    edges=defaultdict(int)
    for cell in cells:
        q=xy[cell];signed=.5*np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1))
        if signed<=0:raise ValueError('nonpositive conforming quadrilateral')
        for i in range(4):edges[tuple(sorted((cell[i],cell[(i+1)%4])))]+=1
    if any(count not in (1,2) for count in edges.values()):raise ValueError('nonmanifold mesh')
    for (i,j),count in edges.items():
        if count==1 and not any(np.allclose(xy[[i,j],k],side,atol=1e-9) for k in (0,1) for side in (-half,half)):
            raise ValueError('unconnected interior edge in conforming mesh')
    return xy,cells,np.array(selected,dtype=bool)


def subdivide_quads(xy, quads, selected, subdivisions):
    """Bilinear subdivision; new boundary vertices stay on the SAME polygon."""
    if type(subdivisions) is not int or subdivisions<1:raise ValueError('positive integer plan refinement')
    if subdivisions==1:return xy,quads,selected
    points=[];mapping={};cells=[];flags=[]
    def node(p):
        key=tuple(np.round(p,11))
        if key not in mapping:mapping[key]=len(points);points.append(p)
        return mapping[key]
    for cell,flag in zip(quads,selected):
        p=xy[cell];grid={}
        for j in range(subdivisions+1):
            v=j/subdivisions
            for i in range(subdivisions+1):
                u=i/subdivisions
                grid[i,j]=node((1-u)*(1-v)*p[0]+u*(1-v)*p[1]+u*v*p[2]+(1-u)*v*p[3])
        for j in range(subdivisions):
            for i in range(subdivisions):cells.append([grid[i,j],grid[i+1,j],grid[i+1,j+1],grid[i,j+1]]);flags.append(flag)
    return np.array(points),np.array(cells),np.array(flags)


def build_conforming_ground(
    config, *, radius, width, kind, sectors=24, radial_subdivisions=1,
    plan_refinement=1, vertical_subdivisions=1, excavation_depth_m=0.0,
):
    """Extrude the footprint mesh and expose its bearing plane at the excavation level."""
    c=validate_config(config)
    if c['width_m']!=c['length_m']:raise ValueError('square exterior required')
    xy,quads,selected=footprint_quads(radius,width,c['width_m'],kind=kind,sectors=sectors,radial_subdivisions=radial_subdivisions)
    xy,quads,selected=subdivide_quads(xy,quads,selected,plan_refinement)
    if type(vertical_subdivisions) is not int or vertical_subdivisions<1:raise ValueError('positive integer vertical refinement')
    levels=depth_levels(c);depths=[levels[0]]
    for lo,hi in zip(levels,levels[1:]):depths.extend(np.linspace(lo,hi,vertical_subdivisions+1)[1:].tolist())
    nz=len(depths)-1;start=c['tag_start'];nxy=len(xy)
    cut = next(
        (k for k, depth in enumerate(depths[:-1]) if math.isclose(depth, excavation_depth_m, abs_tol=1e-9)),
        None,
    )
    if cut is None or excavation_depth_m < 0:
        raise ValueError('excavation base must coincide with a mesh level above the model base')
    retained_plan_nodes = set(quads[~selected].ravel())
    if len(quads)*nz>c['max_elements']:raise ValueError('conforming element budget exceeded')
    if set(range(start,start+nxy*len(depths)))&set(ops.getNodeTags()):raise ValueError('ground node tag collision')
    if set(range(start,start+len(quads)*nz))&set(ops.getEleTags()):raise ValueError('ground element tag collision')
    ops.model('basic','-ndm',3,'-ndf',3)
    names=sorted({v['material'] for v in c['layers']}|{v['material'] for v in c['lenses']})
    mats={n:register(ops,n,start+i,elastic=c['elastic'],stiffness_factor=c['stiffness_factor']) for i,n in enumerate(names)}
    nodes={};index={};fixities={};half=c['width_m']/2
    for k,d in enumerate(depths):
        for i,point in enumerate(xy):
            if k < cut and i not in retained_plan_nodes:
                continue
            t=start+len(nodes);xyz=np.r_[point,-d];nodes[t]=xyz;index[i,k]=t;ops.node(t,*xyz)
            fix=[int(math.isclose(abs(point[0]),half,abs_tol=1e-8)),int(math.isclose(abs(point[1]),half,abs_tol=1e-8)),0]
            if k==nz:fix=[1,1,1]
            fixities[t]=fix
            if any(fix):ops.fix(t,*fix)
    elements={};faces=[];loads=defaultdict(lambda:np.zeros(3));fractions=defaultdict(float);mass=0.;max_aspect=0.
    for k in range(nz):
        for qi,quad in enumerate(quads):
            if k < cut and selected[qi]:
                continue
            tags=[index[i,k+1] for i in quad]+[index[i,k] for i in quad]
            coords=np.array([nodes[t] for t in tags]);geom=brick_geometry_metrics(coords)
            if geom['min_jacobian_det']<=0 or geom['aspect_ratio']>c['max_aspect_ratio']:raise ValueError(f'invalid conforming brick {qi,k}: {geom}')
            center=coords.mean(axis=0);depth=-center[2];name=assigned_material(c,*center[:2],depth);m=mats[name];e=start+len(elements)
            ops.element('bbarBrick',e,*tags,m['tag'],0.,0.,0.)
            w=brick_weights(coords);volume=float(sum(w));rho=m['density_kg_m3'];mass+=volume*rho;fractions[name]+=volume;max_aspect=max(max_aspect,geom['aspect_ratio'])
            for t,v in zip(tags,w):loads[t][2]-=rho*c['gravity_m_s2']*v
            elements[e]=dict(nodes=tags,material=name,volume_m3=volume,center_m=center.tolist(),depth_m=float(depth),height_m=float(depths[k+1]-depths[k]))
            boundary=[]
            top_level = cut if selected[qi] else 0
            if k == top_level:
                boundary.append(('top',FACE_INDEX['top']))
            if k==nz-1:boundary.append(('bottom',FACE_INDEX['bottom']))
            for i in range(4):
                j=(i+1)%4;p=xy[quad[[i,j]]]
                side=next((label for axis,sgn,label in [(0,-1,'xmin'),(0,1,'xmax'),(1,-1,'ymin'),(1,1,'ymax')] if np.allclose(p[:,axis],sgn*half,atol=1e-8)),None)
                if side:boundary.append((side,(i,j,j+4,i+4)))
            for side,local in boundary:
                f=[tags[i] for i in local];fw,normal=face_weights(np.array([nodes[t] for t in f]))
                record=dict(side=side,nodes=f,weights_m2=fw.tolist(),area_m2=float(sum(fw)),normal=normal.tolist(),element=e,material=name)
                if side=='top':record['foundation_selected']=bool(selected[qi])
                faces.append(record)
    volume=sum(fractions.values())
    cut_area = sum(f['area_m2'] for f in faces if f.get('foundation_selected'))
    expected_volume = c['width_m']**2*c['depth_m'] - cut_area*depths[cut]
    if not math.isclose(volume,expected_volume,rel_tol=1e-9):raise RuntimeError('conforming volume check')
    if not math.isclose(-sum(v[2] for v in loads.values()),mass*c['gravity_m_s2'],rel_tol=1e-10):raise RuntimeError('conforming gravity check')
    return dict(config=c,nodes=nodes,index=index,elements=elements,faces=faces,materials=mats,gravity_loads=dict(loads),fixities=fixities,mass_kg=mass,volume_m3=volume,
        material_volume_fractions={n:v/volume for n,v in fractions.items()},
        mesh=dict(kind='square_core_polygon_rings_square_farfield',sectors=sectors,radial_subdivisions=radial_subdivisions,plan_refinement=plan_refinement,vertical_subdivisions=vertical_subdivisions,nz=nz,max_aspect_ratio=max_aspect,depth_levels_m=depths),
        footprint=dict(radius_m=radius,width_m=width,kind=kind,polygon_sectors=sectors),
        boundary='static_box',stage='built',displacement_reference={},attachments=[],boundary_force_loads={})
