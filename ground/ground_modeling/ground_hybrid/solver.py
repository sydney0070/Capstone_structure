"""Checked load-factor continuation; no tolerance relaxation on failure."""
import numpy as np
import openseespy.opensees as ops


def configure_static(tolerance_N):
    ops.wipeAnalysis()
    # Identity equalDOF + homogeneous SP: required by tested contact assembly.
    ops.constraints('Plain');ops.numberer('RCM');ops.system('UmfPack')
    ops.test('NormUnbalance',tolerance_N,80,0);ops.algorithm('NewtonLineSearch')
    ops.integrator('LoadControl',0.);ops.analysis('Static')


def continue_load(pattern_tag, *, steps=20, min_increment=1e-6, progress=None, algorithms=("KrylovNewton","NewtonLineSearch","Newton"), attempt_observer=None, initial_increment=None, max_attempts=10000, first_step_algorithm=None):
    if type(steps) is not int or steps<1:raise ValueError('positive integer steps')
    if type(max_attempts) is not int or max_attempts<1:raise ValueError('positive integer attempt budget')
    if not algorithms or any(a not in ('KrylovNewton','NewtonLineSearch','Newton','NewtonInitial') for a in algorithms):raise ValueError('supported nonempty algorithms required')
    if first_step_algorithm not in (None,'NewtonInitial'):raise ValueError('only explicit NewtonInitial first-step option supported')
    if initial_increment is not None and (not np.isfinite(initial_increment) or not 0<initial_increment<=1./steps):raise ValueError('initial increment must be positive and <=1/steps')
    inc=initial_increment if initial_increment is not None else 1./steps;attempts=[]
    while ops.getLoadFactor(pattern_tag)<1.-1e-10:
        before=ops.getLoadFactor(pattern_tag);step=min(inc,1.-before)
        ops.integrator('LoadControl',step)
        active_algorithms=(first_step_algorithm,)+tuple(a for a in algorithms if a!=first_step_algorithm) if before==0 and first_step_algorithm else algorithms
        for algorithm in active_algorithms:
            if len(attempts)>=max_attempts:raise RuntimeError('load continuation attempt budget exceeded')
            if algorithm=='NewtonInitial':ops.algorithm('Newton','-initial')
            else:ops.algorithm(algorithm)
            if attempt_observer is not None:
                attempt_observer(dict(event='attempt_started',pattern=pattern_tag,before=before,increment=step,algorithm=algorithm))
            ok=ops.analyze(1);after=ops.getLoadFactor(pattern_tag)
            norms=list(ops.testNorm());iterations=ops.testIter()
            final_norm=float(norms[min(iterations,len(norms))-1]) if norms and iterations>0 else None
            if final_norm is not None and (not np.isfinite(final_norm) or (ok and final_norm==0)):final_norm=None
            attempts.append(dict(before=before,increment=step,after=after,code=ok,algorithm=algorithm,iterations=iterations,unbalance_norm_N=final_norm))
            if attempt_observer is not None:attempt_observer(dict(event='attempt_finished',pattern=pattern_tag,**attempts[-1]))
            if ok==0:break
            if abs(before-after)>1e-8:raise RuntimeError('failed step modified committed load factor')
        if ok:
            if abs(before-after)>1e-8:raise RuntimeError('failed step modified committed load factor')
            inc=step/2
            if inc<min_increment:raise RuntimeError(f'load continuation failed at {before:.6g}; no physics/tolerance altered')
        elif abs(after-before-step)>1e-8:raise RuntimeError('unexpected committed load factor')
        elif ops.testIter() < 12:
            inc=min(1./steps,inc*1.25)
        if progress is not None:
            progress(attempts[-1])
    if abs(ops.getLoadFactor(pattern_tag)-1.)>1e-8:raise RuntimeError('load factor overshoot')
    return attempts


def _add_loads(loads):
    for t,f in loads.items():
        ndf=len(ops.nodeDisp(t));v=list(f)
        if len(v)>ndf or len(v)<3 or not np.all(np.isfinite(v)):raise ValueError('finite load vector with3..ndf entries required')
        ops.load(t,*(v+[0.]*(ndf-len(v))))


def force_balance(g):
    ops.reactions()
    # Auxiliary rotational SP constraints may carry bearing artifacts. Only
    # translational support reactions enter this force equilibrium gate.
    support=np.sum([ops.nodeReaction(t)[:3] for t in ops.getNodeTags()],axis=0)
    loads=list(g.get('applied_loads',{}).values())
    applied=np.sum([np.array(v[:3]) for v in loads],axis=0) if loads else np.zeros(3)
    residual=support+applied
    scale=max(sum(np.linalg.norm(v[:3]) for v in loads),1.)
    err=float(np.linalg.norm(residual)/scale)
    disp=np.array([ops.nodeDisp(t)[:3] for t in ops.getNodeTags()])
    if not np.all(np.isfinite(disp)) or err>1e-5:
        raise RuntimeError(f'force balance/finite gate failed: {err}')
    return dict(applied_N=applied.tolist(),reaction_N=support.tolist(),residual_N=residual.tolist(),relative_error=err)


def initialize_gravity(g, *, extra_loads=None, pattern_tag=800001, steps=20):
    """Build soil, raft and structure BEFORE this call; all selfweights once.
    No construction-stage displacement reset. Reference is captured after all
    assembled selfweights; later settlement is incremental to this state.
    """
    if g['stage']!='built':raise ValueError('gravity can be initialized only once')
    pdmy=[m['tag'] for m in g['materials'].values() if m['active_law']=='PressureDependMultiYield02']
    for tag in pdmy:ops.updateMaterialStage('-material',tag,'-stage',0)
    loads={t:np.array(v,dtype=float) for t,v in g['gravity_loads'].items()}
    for t,v in (extra_loads or {}).items():
        if t in loads:
            a=np.zeros(max(len(loads[t]),len(v)));a[:len(loads[t])]+=loads[t];a[:len(v)]+=v;loads[t]=a
        else:loads[t]=np.array(v,dtype=float)
    ops.timeSeries('Linear',pattern_tag);ops.pattern('Plain',pattern_tag,pattern_tag);_add_loads(loads)
    tolerance=max(.001,sum(np.linalg.norm(v[:3]) for v in loads.values())*1e-10)
    configure_static(tolerance)
    stage_attempts=[]
    # A full elastic gravity ramp can leave local stresses outside PDMY yield
    # surfaces beneath loaded frame corners. Use a small elastic seating load,
    # then integrate the actual nonlinear gravity path to the full load factor.
    # This changes the explicitly documented initialization path, NOT strength.
    seating=g['config'].get('pdmy_elastic_seating_fraction',.001)
    if pdmy:
        if not 0<seating<=1:raise ValueError('PDMY seating fraction must be in(0,1]')
        ops.integrator('LoadControl',seating)
        if ops.analyze(1):raise RuntimeError('elastic seating failed')
        for tag in pdmy:ops.updateMaterialStage('-material',tag,'-stage',1)
        ops.integrator('LoadControl',0.)
        for algorithm in ('NewtonLineSearch','KrylovNewton','Newton'):
            ops.algorithm(algorithm)
            ok=ops.analyze(1)
            stage_attempts.append(dict(algorithm=algorithm,code=ok,factor=seating))
            if ok==0:break
        if ok:raise RuntimeError('post-stage-switch equilibrium failed; tolerance and material unchanged')
    attempts=continue_load(pattern_tag,steps=steps)
    ops.loadConst('-time',0.)
    ops.integrator('LoadControl',0.)
    if ops.analyze(1):raise RuntimeError('final gravity re-equilibration failed')
    g['applied_loads']=loads;balance=force_balance(g)
    g['displacement_reference']={t:np.array(ops.nodeDisp(t)) for t in ops.getNodeTags()}
    g['stage']='gravity_equilibrated'
    g['gravity_result']=dict(final_factor=ops.getLoadFactor(pattern_tag),attempts=attempts,
                             post_stage_equilibrated=True,stage_attempts=stage_attempts,force_balance=balance,
                             construction='all assembled selfweights ramped together; PDMY small elastic seating then nonlinear gravity; no stress-free staged installation',
                             pdmy_elastic_seating_fraction=seating if pdmy else None,
                             max_ground_displacement_m=max(np.linalg.norm(ops.nodeDisp(t)) for t in g['nodes']))
    ops.wipeAnalysis()
    return g['gravity_result']


def apply_static_loads(g,loads,*,pattern_tag=800002,steps=20):
    if g['stage'] not in ('gravity_equilibrated','static_loaded'):raise ValueError('initialize gravity first; static solve invalid after absorbing transition')
    ops.timeSeries('Linear',pattern_tag);ops.pattern('Plain',pattern_tag,pattern_tag);_add_loads(loads)
    configure_static(max(.001,sum(np.linalg.norm(v[:3]) for v in g['applied_loads'].values())*1e-10))
    attempts=continue_load(pattern_tag,steps=steps);ops.loadConst('-time',0.)
    for t,v in loads.items():
        old=g['applied_loads'].get(t,np.zeros(len(v)));a=np.zeros(max(len(old),len(v)))
        a[:len(old)]+=old;a[:len(v)]+=v;g['applied_loads'][t]=a
    balance=force_balance(g);g['stage']='static_loaded'
    result=dict(final_factor=ops.getLoadFactor(pattern_tag),attempts=attempts,force_balance=balance,
                max_incremental_ground_displacement_m=max(np.linalg.norm(np.array(ops.nodeDisp(t))-g['displacement_reference'][t]) for t in g['nodes']))
    ops.wipeAnalysis();return result


def initial_stress_summary(g):
    """Point stresses are diagnostics, not a claim of imposed Jaky K0."""
    rows=[]
    for tag,e in g['elements'].items():
        s=np.asarray(ops.eleResponse(tag,'material',1,'stress'),float)
        if len(s)<6 or not np.all(np.isfinite(s)):raise RuntimeError('invalid solid stress response')
        rows.append(dict(element=tag,material=e['material'],depth_m=e['depth_m'],
                         stress_Pa=s[:6].tolist(),Kx=float(s[0]/s[2]) if abs(s[2])>1 else None,
                         Ky=float(s[1]/s[2]) if abs(s[2])>1 else None))
    return rows
