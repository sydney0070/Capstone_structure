"""Unilateral vertical foundation + optional elastic tension ties.
Lateral springs represent INDEPENDENT SHEAR GUIDES, not Coulomb friction.
Elastic small-displacement shell; no anchor failure simulation or buckling.
"""
from pathlib import Path
import sys,json,argparse
import numpy as np
import openseespy.opensees as ops
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'structure_stage2'))
from model import default as structural_default,mesh,surface_loads,get_case,GRAVITY
BASE=Path(__file__).resolve().parent

def default():
    c=structural_default();c.update(support='contact_anchor',horizontal_force_kN=0.,
        pressure_path_kPa=[0.,.025,.05,.075,.1,.2,.4,.6,.8,1.,.6,.2,0.],
        max_pressure_step_kPa=.05,
        anchor=dict(stiffness_per_m=50000.,allowable_tension_per_m=5.),
        scope='UNILATERAL_VERTICAL_CONTACT_WITH_INDEPENDENT_LATERAL_GUIDES')
    c.pop('pressure_kPa');return c

def admissible_contact(xyz,boundary,loads):
    """Necessary global equilibrium condition for nonnegative vertical reactions.
    All support points lie in z=0. Lateral guides do not resist overturning.
    """
    F=loads.sum(axis=0);M=np.cross(xyz,loads).sum(axis=0);R=-F[2]
    if R<=1e-9:return False,'nonpositive_required_compression'
    point=np.array([M[1]/R,-M[0]/R])
    for i,j in zip(boundary,boundary[1:]+boundary[:1]):
        e=xyz[j,:2]-xyz[i,:2];d=point-xyz[i,:2]
        if e[0]*d[1]-e[1]*d[0]<-1e-8:return False,'resultant_outside_support_polygon'
    return True,'necessary_condition_only'

def run(c):
    if c['units']!='kN,m,s,tonne,kPa':raise ValueError('unit mismatch')
    if c['support'] not in ('contact_only','contact_anchor'):raise ValueError('unknown support')
    if c['shape']!='dome':raise ValueError('stage 3 pressure protocol requires a dome cap')
    if c.get('material_config'):raise ValueError('stage 3 currently uses preset + material_override; import calibrated E/nu/rho explicitly')
    path=np.array(c['pressure_path_kPa'],float)
    if len(path)<2 or not np.isfinite(path).all() or path[0]!=0 or min(path)<0:raise ValueError('nonnegative pressure path must start at zero')
    if c['max_pressure_step_kPa']<=0:raise ValueError('positive pressure increment required')
    p=get_case(c['material'])['params'];p.update(c.get('material_override',{}))
    if 'E' not in p or min(p['E'],p['rho'],c['thickness'])<=0 or not -1<p['nu']<.5:raise ValueError('isotropic positive material required')
    if min(c['foundation'].values())<=0 or min(c['anchor'].values())<=0:raise ValueError('positive support values required')
    if c['gravity_factor']<=0:raise ValueError('positive gravity required to seat contact initially')
    xyz,quads,boundary,lengths=mesh(c);area,pvec=surface_loads(xyz,quads)
    mass=area*p['rho']*c['thickness'];weight=mass*GRAVITY[c['planet']]*c['gravity_factor']
    gravity=np.zeros_like(xyz);gravity[:,2]=-weight
    lateral=np.zeros_like(xyz);lateral[:,0]=mass/mass.sum()*c['horizontal_force_kN']
    L=np.array([lengths[i] for i in boundary]);kc=L*c['foundation']['k_vertical_per_m']
    ka=L*c['anchor']['stiffness_per_m'] if c['support']=='contact_anchor' else np.zeros_like(L)
    capacity=L*c['anchor']['allowable_tension_per_m']
    history=[];frames=[];pressure=0.;lateral_factor=0.;stopped=None;attempt=None
    ops.wipe()
    try:
        ops.model('basic','-ndm',3,'-ndf',6)
        for i,pos in enumerate(xyz,1):ops.node(i,*pos)
        ops.section('ElasticMembranePlateSection',1,p['E'],p['nu'],c['thickness'],0.)
        for e,nodes in enumerate(quads,1):ops.element('ShellMITC4',e,*map(int,nodes+1),1)
        ground=[];contact_elements=[];anchor_elements=[]
        for j,i in enumerate(boundary):
            nd=len(xyz)+j+1;ground.append(nd);ops.node(nd,*xyz[i]);ops.fix(nd,1,1,1,1,1,1)
            tag=10*j+1;kt=L[j]*c['foundation']['k_horizontal_per_m']
            ops.uniaxialMaterial('Elastic',tag,kt)
            ops.uniaxialMaterial('ENT',tag+1,float(kc[j]))
            e=len(quads)+2*j+1;contact_elements.append(e)
            ops.element('zeroLength',e,nd,i+1,'-mat',tag,tag,tag+1,'-dir',1,2,3)
            if c['support']=='contact_anchor':
                ops.uniaxialMaterial('Elastic',tag+2,float(ka[j]),0.,0.)
                ops.element('zeroLength',e+1,nd,i+1,'-mat',tag+2,'-dir',3)
                anchor_elements.append(e+1)
        def add_loads(tag,loads):
            ops.timeSeries('Linear',tag);ops.pattern('Plain',tag,tag)
            for i,f in enumerate(loads,1):ops.load(i,*map(float,f),0.,0.,0.)
        add_loads(1,gravity)
        ops.constraints('Transformation');ops.numberer('RCM');ops.system('UmfPack')
        ops.test('NormUnbalance',1e-7,60);ops.algorithm('Newton')
        ops.integrator('LoadControl',.1);ops.analysis('Static')
        if ops.analyze(10)!=0:raise RuntimeError('gravity seating failed')
        ops.loadConst('-time',0.)
        def record(phase,pval,hfactor):
            u=np.array([ops.nodeDisp(i+1)[:3] for i in range(len(xyz))]);z=u[boundary,2]
            ops.reactions();reactions=np.array([ops.nodeReaction(n)[:3] for n in ground])
            contact=-np.array([ops.eleResponse(e,'material',3,'stress')[0] for e in contact_elements])
            anchor=np.array([ops.eleResponse(e,'material',1,'stress')[0] for e in anchor_elements]) if anchor_elements else np.zeros_like(contact)
            external=gravity+lateral*hfactor+pvec*pval
            residual=reactions.sum(axis=0)+external.sum(axis=0)
            moment=np.cross(xyz[boundary],reactions).sum(axis=0)+np.cross(xyz,external).sum(axis=0)
            law_error=max(max(abs(contact-kc*np.maximum(-z,0.))),max(abs(anchor-ka*np.maximum(z,0.))))
            row=dict(phase=phase,pressure_kPa=float(pval),horizontal_load_factor=float(hfactor),
                max_u_m=float(max(np.linalg.norm(u,axis=1))),max_gap_m=float(max(np.maximum(z,0.))),
                contact_length_fraction=float(sum(L[contact>1e-7])/sum(L)),
                contact_compression_kN=float(sum(contact)),anchor_tension_kN=float(sum(anchor)),
                max_anchor_utilization=float(max(anchor/capacity)),
                relative_force_residual=float(np.linalg.norm(residual)/max(np.linalg.norm(external,axis=1).sum(),1.)),
                moment_residual_kNm=float(np.linalg.norm(moment)),constitutive_error_kN=float(law_error),
                max_contact_tension_kN=float(max(np.maximum(-contact,0.))),
                max_anchor_compression_kN=float(max(np.maximum(-anchor,0.))))
            if not np.isfinite(list(v for v in row.values() if isinstance(v,float))).all():raise RuntimeError('nonfinite response')
            history.append(row);frames.append(dict(u=u,contact=contact,anchor=anchor,reactions=reactions))
        record('gravity',0.,0.)
        add_loads(2,lateral)
        for target in np.linspace(0,1,11)[1:]:
            ok,reason=admissible_contact(xyz,boundary,gravity+lateral*target)
            if c['support']=='contact_only' and not ok:
                stopped='NO_STATIC_CONTACT_EQUILIBRIUM';attempt=dict(phase='lateral',factor=float(target),reason=reason);break
            ops.integrator('LoadControl',float(target-lateral_factor))
            if ops.analyze(1)!=0:
                stopped='NUMERICAL_NONCONVERGENCE';attempt=dict(phase='lateral',factor=float(target));break
            lateral_factor=float(target);record('lateral',0.,lateral_factor)
            if history[-1]['max_anchor_utilization']>1.+1e-9:
                stopped='ANCHOR_ALLOWABLE_EXCEEDED';attempt=dict(phase='lateral',factor=lateral_factor);break
        if stopped is None:
            ops.loadConst('-time',0.);add_loads(3,pvec)
            for target in path[1:]:
                start=pressure;count=max(1,int(np.ceil(abs(target-start)/c['max_pressure_step_kPa'])))
                for value in np.linspace(start,target,count+1)[1:]:
                    ok,reason=admissible_contact(xyz,boundary,gravity+lateral+pvec*value)
                    if c['support']=='contact_only' and not ok:
                        stopped='NO_STATIC_CONTACT_EQUILIBRIUM';attempt=dict(phase='pressure',pressure_kPa=float(value),reason=reason);break
                    ops.integrator('LoadControl',float(value-pressure))
                    if ops.analyze(1)!=0:
                        stopped='NUMERICAL_NONCONVERGENCE';attempt=dict(phase='pressure',pressure_kPa=float(value));break
                    pressure=float(value);record('pressure',pressure,1.)
                    if history[-1]['max_anchor_utilization']>1.+1e-9:
                        stopped='ANCHOR_ALLOWABLE_EXCEEDED';attempt=dict(phase='pressure',pressure_kPa=pressure);break
                if stopped:break
        summary=dict(status=stopped or 'COMPLETED_REQUESTED_PATH',attempted_stop=attempt,
            interpretation='NECESSARY_EQUILIBRIUM_CHECK_AND_ELASTIC_TIE_DEMAND_NOT_COLLAPSE_PREDICTION',
            last_converged_pressure_kPa=history[-1]['pressure_kPa'],planet=c['planet'],support=c['support'],
            horizontal_force_kN=c['horizontal_force_kN'],mass_tonne=float(sum(mass)),weight_kN=float(sum(weight)),
            total_anchor_allowable_kN=float(sum(capacity)) if anchor_elements else 0.,
            max_gap_m=max(r['max_gap_m'] for r in history),
            peak_anchor_tension_kN=max(r['anchor_tension_kN'] for r in history),
            peak_anchor_utilization=max(r['max_anchor_utilization'] for r in history),
            max_force_residual=max(r['relative_force_residual'] for r in history),
            max_moment_residual_kNm=max(r['moment_residual_kNm'] for r in history),
            opensees=ops.version(),steps=len(history),mesh_n=c['mesh']['n'])
        frame=frames[-1]
        return summary,history,dict(xyz=xyz,quads=quads,boundary=np.array(boundary),
            boundary_lengths=L,**frame)
    finally:ops.wipe()

def save(c,out):
    s,h,a=run(c);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for name,value in (('config',c),('summary',s),('history',h)):
        (out/f'{name}.json').write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(out/'last_converged.npz',**a)
    return s

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path)
    ap.add_argument('--out',type=Path,default=BASE/'single_run');ap.add_argument('--write-default',type=Path)
    args=ap.parse_args()
    if args.write_default:args.write_default.write_text(json.dumps(default(),indent=2)+'\n')
    else:print(json.dumps(save(json.loads(args.config.read_text()) if args.config else default(),args.out),indent=2))
