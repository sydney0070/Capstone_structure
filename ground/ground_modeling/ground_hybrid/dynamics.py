"""Explicit local viscous far boundaries; fixed-base and incident-wave modes.
These are reference implementations, not a free-field boundary or a validated
site earthquake model. Impedances use record G/K, not evolving tangent moduli.
"""
import math
from collections import defaultdict
import numpy as np
import openseespy.opensees as ops
from capmat.ssi import register_dashpot


def boundary_coefficients(g):
    coeff=defaultdict(lambda:np.zeros(3));bottom=defaultdict(lambda:np.zeros(3))
    for f in g['faces']:
        if f['side']=='top':continue
        m=g['materials'][f['material']];rho=m['density_kg_m3']
        Vs=math.sqrt(m['G_Pa']/rho);Vp=math.sqrt((m['K_Pa']+4*m['G_Pa']/3)/rho)
        axis=2 if f['side']=='bottom' else (0 if f['side'].startswith('x') else 1)
        # Vertical sides and flat bottom only. Tangential impedance added once
        # per incident face; edges/corners intentionally receive both faces.
        if abs(abs(f['normal'][axis])-1)>1e-6:raise ValueError('axis-aligned far faces required')
        for node,area in zip(f['nodes'],f['weights_m2']):
            C=np.full(3,rho*Vs*area);C[axis]=rho*Vp*area;coeff[node]+=C
            if f['side']=='bottom':bottom[node]+=C
    return coeff,bottom


def activate_absorbing_boundaries(g,*,tag_start=600000,pattern_tag=800010):
    """Replace SP supports by their equilibrated forces and viscous elements.
    Must be called at a converged static state. No static solve after release:
    free rigid-body modes are intentional, transient inertia supplies stability.
    """
    if g['stage'] not in ('gravity_equilibrated','static_loaded') or g['boundary']!='static_box':
        raise ValueError('requires equilibrated static box exactly once')
    coeff,bottom=boundary_coefficients(g)
    if set(range(tag_start,tag_start+len(coeff)))&set(ops.getNodeTags()):raise ValueError('boundary node tag collision')
    if set(range(tag_start,tag_start+3*len(coeff)))&set(ops.getEleTags()):raise ValueError('boundary element tag collision')
    ops.loadConst('-time',0.);ops.reactions()
    balancing={}
    for t,fixed in g['fixities'].items():
        r=ops.nodeReaction(t);load=np.zeros(3)
        for d,yes in enumerate(fixed,1):
            if yes:load[d-1]=r[d-1]
        if any(fixed):balancing[t]=load
    ops.timeSeries('Constant',pattern_tag);ops.pattern('Plain',pattern_tag,pattern_tag)
    for t,f in balancing.items():ops.load(t,*f)
    for t,fixed in g['fixities'].items():
        for d,yes in enumerate(fixed,1):
            if yes:ops.remove('sp',t,d)
    ops.model('basic','-ndm',3,'-ndf',3)
    dashpots=[]
    for i,(t,C) in enumerate(sorted(coeff.items())):
        anchor=tag_start+i;ops.node(anchor,*ops.nodeCoord(t));ops.fix(anchor,1,1,1)
        for d,c in enumerate(C,1):
            eid=tag_start+3*i+d-1
            dashpots.append(register_dashpot(ops,eid,eid,t,anchor,float(c),d))
    g['boundary_force_loads']=balancing;g['boundary']='local_viscous'
    g['stage']='dynamic_ready';g['dynamic_boundary']=dict(dashpots=dashpots,bottom_coefficients=dict(bottom),coefficients=dict(coeff))
    return dict(kind='local_viscous',dashpots=len(dashpots),
                released_reaction_force_N=np.sum(list(balancing.values()),axis=0).tolist(),
                scope='plane-wave approximation with fixed record impedances; not free-field, no automatic Rayleigh damping')


def mesh_frequency_guidance(g,frequency_Hz):
    if not math.isfinite(frequency_Hz) or frequency_Hz<=0:raise ValueError('positive target frequency')
    worst=math.inf;max_dt=math.inf
    for e in g['elements'].values():
        xyz=np.array([g['nodes'][n] for n in e['nodes']]);m=g['materials'][e['material']]
        edges=[np.linalg.norm(xyz[a]-xyz[b]) for a,b in [(0,1),(0,3),(0,4),(2,6),(4,5),(4,7)]]
        Vs=math.sqrt(m['G_Pa']/m['density_kg_m3']);Vp=math.sqrt((m['K_Pa']+4*m['G_Pa']/3)/m['density_kg_m3'])
        worst=min(worst,Vs/(frequency_Hz*max(edges)))
        max_dt=min(max_dt,min(edges)/(10*Vp))
    return dict(target_frequency_Hz=frequency_Hz,min_elements_per_shear_wavelength=worst,
                suggested_dt_s=min(max_dt,1/(20*frequency_Hz)),mesh_target_pass=worst>=10,
                scope='10 elements/wavelength and temporal resolution heuristics; implicit stability does not imply accuracy; uses record stiffness')


def add_incident_base_velocity(g,velocity_m_s,dt_s,*,axis=1,tag=800020):
    """Upgoing incident velocity at viscous base: traction=2*rho*c*v_inc.
    Input is NOT outcrop acceleration. No automatic factor-of-two conversion
    from arbitrary records, no simultaneous UniformExcitation on released base.
    """
    if g['boundary']!='local_viscous':raise ValueError('requires released viscous base')
    v=np.asarray(velocity_m_s,float)
    if axis not in (1,2,3) or not math.isfinite(dt_s) or dt_s<=0 or v.ndim!=1 or len(v)<2 or not np.all(np.isfinite(v)) or abs(v[0])>1e-12:
        raise ValueError('axis1/2/3, positive dt, finite velocity beginning at0 required')
    ops.timeSeries('Path',tag,'-dt',dt_s,'-values',*v.tolist())
    ops.pattern('Plain',tag,tag)
    for t,C in g['dynamic_boundary']['bottom_coefficients'].items():
        f=np.zeros(3);f[axis-1]=2*C[axis-1];ops.load(t,*f)
    return dict(input='incident velocity m/s',axis=axis,dt_s=dt_s,samples=len(v),traction_factor=2)


def run_transient(g,*,dt_s,steps,record_nodes,rayleigh=None,
                  visualization_nodes=None,max_visualization_frames=151):
    if g['stage'] not in ('gravity_equilibrated','static_loaded','dynamic_ready','dynamic_ran'):
        raise ValueError('equilibrated initial state required')
    if not math.isfinite(dt_s) or dt_s<=0 or type(steps) is not int or steps<1:raise ValueError('positive dt/steps')
    record_nodes=list(record_nodes)
    if not record_nodes or not set(record_nodes)<=set(ops.getNodeTags()):raise ValueError('existing record nodes required')
    visual=None;sample_steps=set()
    if visualization_nodes is not None:
        visual_nodes=sorted(set(visualization_nodes))
        if not visual_nodes or not set(visual_nodes)<=set(ops.getNodeTags()):raise ValueError('existing visualization nodes required')
        if type(max_visualization_frames) is not int or max_visualization_frames<2:raise ValueError('at least2 visualization frames')
        # Store actual committed solutions, including the pre-transient state.
        # Subsampling affects display only, never integration dt or probe rows.
        sample_steps=set(np.linspace(0,steps,min(steps+1,max_visualization_frames),dtype=int))
        initial=[list(ops.nodeDisp(t)[:3]) for t in visual_nodes]
        visual=dict(node_tags=visual_nodes,time_s=[ops.getTime()],displacement_m=[initial],
                    reference_displacement_m=initial,integration_steps=steps,
                    scope='sampled committed nodal translations; sampled maxima may miss between-frame peaks')
    if rayleigh is not None:
        zeta,f1,f2=rayleigh
        if not 0<=zeta<.2 or not 0<f1<f2:raise ValueError('Rayleigh(zeta,f1,f2) invalid')
        w1,w2=2*math.pi*f1,2*math.pi*f2
        # Initial stiffness proportional term avoids evolving negative tangent
        # contribution, but can overdamp high frequencies. Explicit opt-in only.
        ops.rayleigh(2*zeta*w1*w2/(w1+w2),0.,2*zeta/(w1+w2),0.)
    ops.wipeAnalysis();ops.constraints('Plain');ops.numberer('RCM');ops.system('UmfPack')
    ops.test('NormUnbalance',max(.001,g['mass_kg']*g['config']['gravity_m_s2']*1e-10),80,0)
    ops.algorithm('Newton');ops.integrator('Newmark',.5,.25);ops.analysis('Transient')
    rows=[];attempts=[]
    def advance(dt,level=0):
        before=ops.getTime();ok=ops.analyze(1,dt);after=ops.getTime()
        attempts.append(dict(before=before,dt=dt,after=after,code=ok))
        if ok:
            if abs(after-before)>1e-9:raise RuntimeError('failed transient step changed committed time')
            if level>=8:raise RuntimeError('transient failed after8 subdivisions')
            advance(dt/2,level+1);advance(dt/2,level+1)
    for k in range(steps):
        advance(dt_s)
        data={t:dict(u=list(ops.nodeDisp(t)),v=list(ops.nodeVel(t)),a=list(ops.nodeAccel(t))) for t in record_nodes}
        if not all(np.all(np.isfinite([*d['u'],*d['v'],*d['a']])) for d in data.values()):raise RuntimeError('nonfinite transient')
        rows.append(dict(time_s=ops.getTime(),nodes=data))
        if visual is not None and k+1 in sample_steps:
            u=[list(ops.nodeDisp(t)[:3]) for t in visual_nodes]
            if not np.all(np.isfinite(u)):raise RuntimeError('nonfinite visualization displacement')
            visual['time_s'].append(ops.getTime());visual['displacement_m'].append(u)
    g['stage']='dynamic_ran';ops.wipeAnalysis()
    return dict(dt_s=dt_s,steps=steps,rows=rows,attempts=attempts,boundary=g['boundary'],
                damping='none' if rayleigh is None else list(rayleigh),visualization=visual)
