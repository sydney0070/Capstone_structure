"""Draft Mars raft / Moon annular-shell foundations with explicit dome connectors.

SI, local ENU. Geometry/material/anchor values are declared research assumptions.
No pressure-vessel closure, construction-stage or site bearing-capacity certification.
Builders append; caller owns the OpenSees domain and gravity analysis.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import math,sys
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'ground/ground_modeling'),str(ROOT/'material/material_release_20261101')]
from ground_hybrid.mesh import face_weights
from SSI.contact import register_contact
from SSI.parameters import contact_parameters

@dataclass(frozen=True)
class FoundationConfig:
    body:str='Mars'
    kind:str='raft'
    radius_m:float=20.
    width_m:float=4.
    thickness_m:float=.6
    E_Pa:float=30e9
    nu:float=.2
    density_kg_m3:float=2400.
    mu:float=.35
    penalty_factor:float=10.
    max_area_relative_error:float=.1
    tag_start:int=300000
    # Supplied only when a deliberate tension-connector experiment is requested.
    anchor_stiffness_N_m:float|None=None
    anchor_capacity_N:float|None=None


def add_foundation(g,config,ports_xy):
    """Area-resolved ShellMITC4 + friction/uplift at all supporting mesh vertices.

    Boundary follows selected Cartesian faces, not an exact circular perimeter.
    Face-centre mask selects an annulus/disk; footprint area error is reported.
    Shell midsurface/contact plane coincidence omits half-thickness eccentricity.
    Foundation is elastic nominal homogeneous prototype, no inferred lunar RC.
    """
    c=config
    if c.body not in ('Mars','Moon') or c.kind not in ('raft','ring'):raise ValueError('body/kind')
    for key in ['radius_m','width_m','thickness_m','E_Pa','density_kg_m3','penalty_factor']:
        if not math.isfinite(getattr(c,key)) or getattr(c,key)<=0:raise ValueError(f'positive finite {key}')
    if not -1<c.nu<.5 or not math.isfinite(c.mu) or c.mu<0 or c.width_m>=2*c.radius_m:raise ValueError('invalid Poisson/friction/ring width')
    if g['stage']!='built' or g['attachments']:raise ValueError('one foundation before gravity')
    anchor=(c.anchor_stiffness_N_m,c.anchor_capacity_N)
    if any(v is not None for v in anchor) and not all(v is not None and math.isfinite(v) and v>0 for v in anchor):raise ValueError('both positive anchor stiffness and capacity required')
    faces=[];outer=c.radius_m+c.width_m/2;inner=max(0,c.radius_m-c.width_m/2)
    if 'footprint' in g:
        fp=g['footprint']
        if fp['kind']!=c.kind or not math.isclose(fp['radius_m'],c.radius_m) or not math.isclose(fp['width_m'],c.width_m):raise ValueError('conforming mesh footprint/config mismatch')
    for f in g['faces']:
        if f['side']!='top':continue
        xyz=np.array([g['nodes'][n] for n in f['nodes']]);r=np.linalg.norm(xyz[:,:2].mean(axis=0))
        choose=f.get('foundation_selected',r<=outer if c.kind=='raft' else inner<=r<=outer)
        if choose:faces.append(f)
    if not faces:raise ValueError('no resolved footprint cells')
    target=math.pi*outer**2 if c.kind=='raft' else math.pi*(outer**2-inner**2)
    selected_area=sum(f['area_m2'] for f in faces)
    if not math.isfinite(c.max_area_relative_error) or not 0<c.max_area_relative_error<1:raise ValueError('area error limit in (0,1)')
    if abs(selected_area-target)/target>c.max_area_relative_error:
        raise ValueError(f'Footprint area error {abs(selected_area-target)/target:.1%} exceeds {c.max_area_relative_error:.1%}; refine mesh, do not treat coarse cells as exact ring/disk')
    if g['config']['body']!=c.body:raise ValueError('foundation/soil body mismatch')
    nodes=sorted({t for f in faces for t in f['nodes']})
    if np.ptp([g['nodes'][n][2] for n in nodes])>1e-8:raise ValueError('requires explicitly horizontal prepared/reference plane')
    B=c.tag_start
    if set(range(B,B+3*len(nodes)+1))&set(ops.getNodeTags()):raise ValueError('foundation node tag collision')
    if set(range(B,B+len(faces)+3*len(nodes)+1))&set(ops.getEleTags()):raise ValueError('foundation element tag collision')
    ops.model('basic','-ndm',3,'-ndf',6)
    mapping={t:B+i for i,t in enumerate(nodes)}
    for t,n in mapping.items():ops.node(n,*g['nodes'][t])
    # Shell density=0; all mass is explicitly nodal, once, for auditable inertia.
    ops.section('ElasticMembranePlateSection',B,c.E_Pa,c.nu,c.thickness_m,0.)
    weights={t:0. for t in nodes};local={t:np.zeros(3) for t in nodes};elements=[]
    for i,f in enumerate(faces):
        e=B+i;ops.element('ShellMITC4',e,*[mapping[t] for t in f['nodes']],B);elements.append(e)
        mat=g['materials'][f['material']];height=g['elements'][f['element']]['height_m']
        for t,area in zip(f['nodes'],f['weights_m2']):weights[t]+=area;local[t]+=area*np.array([mat['E_Pa'],mat['G_Pa'],height])
    contacts=[];mass={};loads={};anchors=[]
    for i,t in enumerate(nodes):
        n=mapping[t];area=weights[t];E,G,h=local[t]/area
        p=contact_parameters(area,E,G,h,50000.,.001,c.mu,penalty_factor=c.penalty_factor)
        contacts.append(register_contact(ops,B+len(faces)+i,n,t,p,auxiliary_nodes=(B+len(nodes)+2*i,B+len(nodes)+2*i+1),normal_material_tag=B+10000+2*i,rotation_material_tag=B+10001+2*i,friction_tag=B+20000+i,restore_ndf=6))
        m=area*c.thickness_m*c.density_kg_m3;mass[n]=m;ops.mass(n,m,m,m,0.,0.,0.);loads[n]=[0.,0.,-m*g['config']['gravity_m_s2'],0.,0.,0.]
    ports=[]
    for xy in ports_xy:
        matches=[mapping[t] for t in nodes if np.linalg.norm(g['nodes'][t][:2]-xy)<1e-7]
        if len(matches)!=1:raise ValueError('foundation grid must resolve every exact dome port; no snapping')
        ports.append(matches[0])
    if len(set(ports))!=len(ports):raise ValueError('unique ports required')
    if anchor[0] is not None:
        reverse={v:k for k,v in mapping.items()}
        for i,n in enumerate(ports):
            # Positive extension only, bilinear perfect-plastic tie. This is NOT
            # an embedded anchor/pullout model; capacity is an explicit input.
            tag=B+40000+i
            # 3DOF ground cannot directly share a 6DOF zeroLength: use existing
            # ground-side auxiliary node with translations tied to physical soil.
            ci=nodes.index(reverse[n]);soil_aux=contacts[ci]['auxiliary_nodes'][0]
            register_tension_anchor(soil_aux,n,tag=tag,stiffness_N_m=anchor[0],capacity_N=anchor[1])
            anchors.append(dict(element=tag,foundation_node=n,ground_node=reverse[n],capacity_N=anchor[1],stiffness_N_m=anchor[0]))
    area=sum(weights.values());target=math.pi*outer**2 if c.kind=='raft' else math.pi*(outer**2-inner**2)
    result=dict(config=asdict(c),nodes=list(mapping.values()),elements=elements,contacts=contacts,anchors=anchors,dome_ports=ports,
        nodal_mass_kg=mass,gravity_loads=loads,area_m2=area,target_area_m2=target,area_relative_error=(area-target)/target,
        total_mass_kg=sum(mass.values()),weights_m2={mapping[t]:a for t,a in weights.items()},
        scope='draft elastic foundation material; nonlinear soil and friction/uplift; '+('fixed polygon footprint' if 'footprint' in g else 'staircase footprint')+'; zero eccentricity; no shear failure, pressure closure or anchor pullout calibration')
    g['attachments'].append(result);return result


from SSI.connectors_20261109 import connect_dome_port,register_tension_anchor
