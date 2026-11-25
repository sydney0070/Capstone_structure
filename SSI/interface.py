"""Area-weighted 3DOF soil to 6DOF shell raft contacts and structure ports."""
import math
import numpy as np
import openseespy.opensees as ops
from ground_hybrid.mesh import face_weights
from ground_hybrid.materials import RELEASE
from .parameters import contact_parameters
from .contact import register_contact
from capmat.elastic_rc import register_elastic_rc


def footprint_faces(g, bounds):
    xmin,xmax,ymin,ymax=map(float,bounds)
    if not xmin<xmax or not ymin<ymax:raise ValueError('ordered rectangular bounds required')
    selected=[]
    for f in g['faces']:
        if f['side']!='top':continue
        xyz=np.array([g['nodes'][t] for t in f['nodes']])
        inside=(xyz[:,0]>=xmin-1e-8)&(xyz[:,0]<=xmax+1e-8)&(xyz[:,1]>=ymin-1e-8)&(xyz[:,1]<=ymax+1e-8)
        if np.all(inside):selected.append(f)
    expected=(xmax-xmin)*(ymax-ymin)
    if not selected:raise ValueError('no complete mesh faces in footprint; align/refine mesh')
    z=np.array([g['nodes'][t][2] for f in selected for t in f['nodes']])
    if np.ptp(z)>1e-8:raise ValueError('flat foundation pad required; grade terrain explicitly before building')
    area=sum(f['area_m2'] for f in selected)
    if abs(area-expected)>expected*1e-8:raise ValueError('footprint edges must match grid; partial faces cannot be silently omitted')
    return selected


def add_raft(g,bounds=(-1.,1.,-1.,1.),*,thickness_m=.4,mu=None,penalty_factor=10.,
             pressure_reference_Pa=50000.,allowed_penetration_m=1e-4,tag_start=300000,
             reinforcement_ratios=(.01,.01,0.)):
    """Illustrative elastic RC ShellMITC4 raft; dimensions are not a design.
    Reference midsurface is at contact plane (zero eccentricity idealization).
    Must be added BEFORE gravity. All six raft DOFs are valid structure ports.
    """
    if g['stage']!='built':raise ValueError('attach raft before gravity; staged installation not implemented')
    if not math.isfinite(thickness_m) or thickness_m<=0:raise ValueError('positive thickness')
    faces=footprint_faces(g,bounds);ground_nodes=sorted({t for f in faces for t in f['nodes']})
    if g['attachments']:raise ValueError('one raft helper per ground; use explicit port API for multiple foundations')
    n=len(ground_nodes)
    if set(range(tag_start,tag_start+3*n+10))&set(ops.getNodeTags()):raise ValueError('raft node tag collision')
    if set(range(tag_start,tag_start+n+len(faces)+10))&set(ops.getEleTags()):raise ValueError('raft element tag collision')
    ops.model('basic','-ndm',3,'-ndf',6)
    mapping={}
    for i,t in enumerate(ground_nodes):
        new=tag_start+i;ops.node(new,*g['nodes'][t]);mapping[t]=new
    # Initial elastic bonded RC continuum using the previous verified package.
    rc=register_elastic_rc(ops,tag_start,reinforcement_volume_fractions=reinforcement_ratios)
    ops.nDMaterial('PlateFiber',tag_start+1,tag_start)
    ops.section('PlateFiber',tag_start+2,tag_start+1,thickness_m)
    rho=2400*(1-sum(reinforcement_ratios))+7850*sum(reinforcement_ratios)
    shell_tags=[];weights={t:0. for t in ground_nodes};local={t:np.zeros(3) for t in ground_nodes}
    for i,f in enumerate(faces):
        eid=tag_start+i;shell_tags.append(eid)
        ops.element('ShellMITC4',eid,*[mapping[t] for t in f['nodes']],tag_start+2)
        m=g['materials'][f['material']];h=g['elements'][f['element']]['height_m']
        for t,a in zip(f['nodes'],f['weights_m2']):
            weights[t]+=a;local[t]+=a*np.array([m['E_Pa'],m['G_Pa'],h])
    if mu is None:
        # Mars value is measured JEZ-1/smooth-steel; Earth/Moon same numerical
        # value is an explicit interface assumption, not a soil friction angle.
        mu=math.tan(math.radians(18.9)) if g['config']['body']=='Mars' else .35
    contacts=[]
    for i,t in enumerate(ground_nodes):
        area=weights[t];E,G,h=local[t]/area
        p=contact_parameters(area,E,G,h,pressure_reference_Pa,allowed_penetration_m,mu,penalty_factor=penalty_factor)
        eid=tag_start+len(faces)+i
        contact=register_contact(ops,eid,mapping[t],t,p,
            auxiliary_nodes=(tag_start+n+2*i,tag_start+n+2*i+1),
            normal_material_tag=tag_start+10+2*i,rotation_material_tag=tag_start+11+2*i,
            friction_tag=tag_start+i,restore_ndf=6)
        contacts.append(contact)
        g['gravity_loads'][mapping[t]]=np.array([0.,0.,-rho*thickness_m*area*g['config']['gravity_m_s2'],0.,0.,0.])
    raft=dict(nodes=list(mapping.values()),ground_to_raft=mapping,weights_m2={mapping[t]:a for t,a in weights.items()},
              shell_elements=shell_tags,contacts=contacts,bounds_m=list(bounds),thickness_m=thickness_m,
              mass_kg=rho*thickness_m*sum(weights.values()),mu=mu,reinforcement_ratios=list(reinforcement_ratios),
              material=rc,scope='elastic RC, zero offset contact plane, fixed-normal point friction; no slab cracking or staged installation',
              interface_evidence='Mars: JEZ-1 smooth steel friction transferred to raft surface as assumption; Earth/Moon: assumed0.35; no adhesion')
    g['attachments'].append(raft);return raft


def contact_summary(raft):
    rows=[]
    for c in raft['contacts']:
        forces=np.array(ops.eleForce(c['element_tag']));r=forces[6:9]
        ground=np.array(ops.nodeDisp(c['ground_node'])[:3]);foundation=np.array(ops.nodeDisp(c['foundation_node'])[:3])
        auxiliary=[ops.nodeDisp(a)[:3] for a in c['auxiliary_nodes']]
        if not all(np.isfinite(v).all() for v in (forces,ground,foundation,*auxiliary)):
            raise RuntimeError(f"Nonfinite contact response: element={c['element_tag']} time={ops.getTime()}")
        relative=foundation-ground
        rows.append(dict(element=c['element_tag'],ground_node=c['ground_node'],raft_node=c['foundation_node'],
                         area_m2=c['area_m2'],normal_N=float(-r[2]),shear_N=float(np.linalg.norm(r[:2])),
                         pressure_Pa=float(-r[2]/c['area_m2']),relative_displacement_m=relative.tolist(),
                         friction_limit_N=c['mu']*max(-r[2],0.),
                         mpc_error_m=float(max(np.linalg.norm(np.array(ops.nodeDisp(a)[:3])-ops.nodeDisp(p)[:3]) for a,p in zip(c['auxiliary_nodes'],[c['ground_node'],c['foundation_node']])))))
    return dict(total_normal_N=sum(r['normal_N'] for r in rows),
                max_penetration_m=max(max(-r['relative_displacement_m'][2],0.) for r in rows),
                max_mpc_error_m=max(r['mpc_error_m'] for r in rows),rows=rows)


def raft_response(g,raft):
    rows=[]
    for t in raft['nodes']:
        u=np.array(ops.nodeDisp(t));reference=g['displacement_reference'].get(t,np.zeros(6))
        rows.append(dict(node=t,coordinate_m=list(ops.nodeCoord(t)),displacement_m_rad=u.tolist(),
                         additional_settlement_m=float(-(u[2]-reference[2]))))
    values=[r['additional_settlement_m'] for r in rows]
    mean=sum(r['additional_settlement_m']*raft['weights_m2'][r['node']] for r in rows)/sum(raft['weights_m2'].values())
    return dict(mean_additional_settlement_m=float(mean),mean_convention='tributary-area weighted',differential_settlement_m=max(values)-min(values),
                contact=contact_summary(raft),nodes=rows)


def pressure_loads(raft,pressure_Pa):
    if not math.isfinite(pressure_Pa):raise ValueError('finite pressure')
    return {t:[0.,0.,-pressure_Pa*a,0.,0.,0.] for t,a in raft['weights_m2'].items()}


def connect_matching_nodes(raft,pairs,*,translations_only=False):
    """Prefer reusing existing raft nodes. Optional coincident duplicate ports.
    Plain handler only: no rigidLink/rigidDiaphragm or chained MP constraints.
    Caller must ensure duplicate structure nodes are unconstrained otherwise.
    """
    pairs=list(pairs)  # Validate and then apply even when caller supplies a generator.
    seen=set()
    for raft_node,structure_node in pairs:
        if raft_node not in raft['nodes'] or structure_node in seen or structure_node in raft['nodes']:
            raise ValueError('unique external structure node paired to a raft node required')
        if len(ops.nodeDisp(structure_node))!=6:raise ValueError('structure node must have6DOF')
        if np.linalg.norm(np.array(ops.nodeCoord(raft_node))-ops.nodeCoord(structure_node))>1e-8:
            raise ValueError('coincident nodes required; no nearest-node snapping')
        seen.add(structure_node)
    for a,b in pairs:ops.equalDOF(a,b,*(range(1,4) if translations_only else range(1,7)))
    return dict(pairs=[list(p) for p in pairs],handler='Plain',dofs=list(range(1,4) if translations_only else range(1,7)))


def displacement_path(g,raft,control_node,axis,targets_m,*,pattern_tag=800040,steps_per_segment=20):
    """EXPERIMENTAL quasi-static translation-controlled diagnostic.
    Assembled large-slip reversal has NOT passed; do not use for design capacity.
    targets are increments from the displacement at entry. A single-control
    equilibrium path, NOT prescribed multi-support earthquake motion.
    """
    from .solver import configure_static,force_balance
    from capmat.analysis import advance_displacement
    if g['stage'] not in ('gravity_equilibrated','static_loaded'):raise ValueError('equilibrated static model required')
    if control_node not in raft['nodes'] or axis not in (1,2,3):raise ValueError('raft translation control required')
    targets=np.asarray(targets_m,dtype=float)
    if targets.ndim!=1 or not len(targets) or not np.all(np.isfinite(targets)) or type(steps_per_segment) is not int or steps_per_segment<1:
        raise ValueError('finite targets and positive integer steps')
    ops.timeSeries('Linear',pattern_tag);ops.pattern('Plain',pattern_tag,pattern_tag)
    area=sum(raft['weights_m2'].values());reference={}
    for t,a in raft['weights_m2'].items():
        f=np.zeros(6);f[axis-1]=a/area;reference[t]=f;ops.load(t,*f)
    configure_static(max(.001,g['mass_kg']*g['config']['gravity_m_s2']*1e-10))
    offset=ops.nodeDisp(control_node,axis);current=offset;rows=[]
    # Recursive subdivision uses the same convergence tolerance and committed
    # state check, just as the existing RC displacement controller.
    for target in targets:
        end=offset+target
        for value in np.linspace(current,end,steps_per_segment+1)[1:]:
            if not advance_displacement(ops,control_node,axis,float(value)):
                raise RuntimeError('controlled displacement failed after subdivision')
            rows.append(dict(control_increment_m=float(ops.nodeDisp(control_node,axis)-offset),
                             total_reference_force_N=ops.getLoadFactor(pattern_tag),contact=contact_summary(raft)))
        current=end
    factor=ops.getLoadFactor(pattern_tag);ops.loadConst('-time',0.)
    for t,f in reference.items():
        old=g['applied_loads'].get(t,np.zeros(6));v=np.zeros(6);v[:len(old)]+=old;v+=factor*f;g['applied_loads'][t]=v
    balance=force_balance(g);g['stage']='static_loaded';ops.wipeAnalysis()
    return dict(axis=axis,control_node=control_node,rows=rows,force_balance=balance,
                scope='single controlled translation; simultaneous imposed supports not implemented')
