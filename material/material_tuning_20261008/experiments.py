"""Independent specimens. Calling simulate() wipes the OpenSees domain.
Stress/force is positive in tension; engineering shear strain is dimensionless.
"""
import numpy as np
import openseespy.opensees as ops
from materials import uniaxial, continuum
from presets import validate

def path(c):
    t = c['test']
    a = t['amplitude']
    if c['material'] == 'Concrete02':
        peaks = [0., .0001, 0., -.001, 0., -.002, 0., -.004, 0., -.008, 0.]
    elif c['kind'] == 'viscous':
        # x is prescribed strain/deformation RATE for this isolated material test.
        return .1*np.sin(np.linspace(0,6*np.pi,6*t['steps']+1))
    else:
        peaks = [0.] + [v for _ in range(t['cycles']) for v in (a,0.,-a,0.)]
    return np.r_[0., np.concatenate([np.linspace(lo,hi,t['steps']+1)[1:]
                  for lo,hi in zip(peaks[:-1],peaks[1:])])]

def static_analysis():
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr',1.e-10,60)
    ops.algorithm('Newton')
    ops.integrator('LoadControl',1.)
    ops.analysis('Static')

def checked(n=1,dt=None):
    result = ops.analyze(n) if dt is None else ops.analyze(n,dt)
    if result != 0:
        raise RuntimeError(f'OpenSees did not converge: code={result}, time={ops.getTime()}')

def imposed_series(tag, values, dt=1., start=0.):
    times = start + np.arange(len(values))*dt
    ops.timeSeries('Path',tag,'-time',*times,'-values',*values,'-useLast')
    ops.pattern('Plain',tag,tag)

def simulate(c, x=None):
    validate(c)
    x = path(c) if x is None else np.asarray(x,dtype=float)
    if len(x)<2 or not np.isfinite(x).all() or abs(x[0]) > 1e-12:
        raise ValueError('history requires >=2 finite samples, starting at zero')
    # Subdivide along HISTORY, never sort strain: sorting destroys hysteresis.
    sub = int(c['test']['substeps'])
    dense = np.r_[x[0],np.concatenate([np.linspace(a,b,sub+1)[1:]
                     for a,b in zip(x[:-1],x[1:])])]
    ops.wipe()
    dispatch = dict(uniaxial=uniaxial_test,spring=uniaxial_test,soil3d=soil3d_test,
                    soil2dcv=constant_volume_test,contact=contact_test,viscous=viscous_test,
                    shell=shell_test)
    try:
        result = dispatch[c['kind']](c,dense)
        result = {k: np.asarray(v)[::sub] for k,v in result.items()}
        if any(not np.isfinite(v).all() for v in result.values()):
            raise RuntimeError('nonfinite response')
        return result
    finally:
        ops.wipe()

def uniaxial_test(c,x):
    ops.model('basic','-ndm',1,'-ndf',1)
    ops.node(1,0.); ops.node(2,1. if c['kind']=='uniaxial' else 0.)
    ops.node(3,1.)
    ops.fix(1,1)
    uniaxial(1,c)
    if c['kind']=='uniaxial':
        ops.element('truss',1,1,2,1.,1)  # L=A=1: force/A equals stress
    else:
        ops.element('zeroLength',1,1,2,'-mat',1,'-dir',1)
    # Separate unloaded elastic DOF lets the solver handle a fully imposed coupon.
    # It is not connected to node 2, so no artificial specimen stiffness is added.
    ops.uniaxialMaterial('Elastic',99,1.)
    ops.element('truss',99,1,3,1.,99)
    imposed_series(1,x)
    ops.sp(2,1,1.)
    static_analysis()
    y=[0.]; reaction=[0.]
    for _ in x[1:]:
        checked()
        ops.reactions()
        y.append(ops.eleResponse(1,'material',1,'stress')[0]
                 if c['kind']=='spring' else ops.eleResponse(1,'material','stress')[0])
        reaction.append(ops.nodeReaction(2,1))
    return dict(x=x,y=y,reaction=reaction)

def soil3d_test(c,x):
    ops.model('basic','-ndm',3,'-ndf',3)
    xyz=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),
         (0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    for i,coord in enumerate(xyz,1):
        ops.node(i,*coord)
        ops.fix(i,*( (1,1,1) if i<=4 else (0,1,0) ))
    for i in (6,7,8):
        ops.equalDOF(5,i,1,3)
    continuum(1,c,3)
    ops.element('SSPbrick',1,*range(1,9),1,0.,0.,0.)
    nonlinear = c['material'] in ('PressureDependMultiYield02','PressureIndependMultiYield')
    if nonlinear:
        ops.updateMaterialStage('-material',1,'-stage',0)
    ops.timeSeries('Linear',1); ops.pattern('Plain',1,1)
    ops.load(5,0.,0.,-c['test']['sigma_v'])
    static_analysis(); ops.integrator('LoadControl',.05)
    checked(20)
    ops.loadConst('-time',0.)
    if nonlinear:
        ops.updateMaterialStage('-material',1,'-stage',1)
        ops.integrator('LoadControl',0.); checked()
    initial=np.array(ops.eleResponse(1,'strain'))
    ops.timeSeries('Linear',2); ops.pattern('Plain',2,2); ops.load(5,1.,0.,0.)
    rows=[]
    def record(target):
        s=np.array(ops.eleResponse(1,'stress'))
        e=np.array(ops.eleResponse(1,'strain'))
        ops.reactions()
        rows.append([target,s[5],np.sum(e[:3]-initial[:3]),s[2],
                     -sum(ops.nodeReaction(i,1) for i in (1,2,3,4)),
                     ops.nodeDisp(5,1),-np.mean(s[:3])])
    record(0.)
    for target in x[1:]:
        dx=float(target)-ops.nodeDisp(5,1)
        if abs(dx)>1e-15:
            ops.integrator('DisplacementControl',5,1,dx)
            checked()
        record(target)
    a=np.array(rows)
    return dict(zip(('x','y','eps_v','sigma_vertical','reaction','actual_x','p_mean'),a.T))

def constant_volume_test(c,x):
    """Plane-strain, constant-volume material proxy. No pore fluid DOFs.
    ru_equivalent is effective vertical stress loss, not simulated fluid pressure.
    """
    t=c['test']
    ops.model('basic','-ndm',2,'-ndf',2)
    for n,coord in enumerate(((0,0),(1,0),(1,1),(0,1)),1):
        ops.node(n,*coord)
        ops.fix(n,1,1 if n<3 else 0)
    continuum(1,c,2)
    ops.element('SSPquad',1,1,2,3,4,1,'PlaneStrain',1.,0.,0.)
    ops.updateMaterialStage('-material',1,'-stage',0)
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
    for n in (3,4):ops.load(n,0.,-t['sigma_v']/2.)
    static_analysis();ops.integrator('LoadControl',.05);checked(20)
    ops.loadConst('-time',0.)
    ops.timeSeries('Constant',2);ops.pattern('Plain',2,2)
    for n in (3,4):ops.sp(n,2,ops.nodeDisp(n,2))
    ops.updateMaterialStage('-material',1,'-stage',1)
    ops.setParameter('-val',0,'-ele',1,'FirstCall','1')
    for n in (3,4):ops.remove('sp',n,1)
    # Independent unloaded DOF, no artificial stiffness at specimen nodes.
    ops.node(99,0.,0.);ops.fix(99,0,1)
    ops.uniaxialMaterial('Elastic',99,1.)
    ops.element('zeroLength',99,1,99,'-mat',99,'-dir',1)
    imposed_series(3,x)
    for n in (3,4):ops.sp(n,1,1.)
    ops.wipeAnalysis();static_analysis()
    ops.integrator('LoadControl',0.);checked()
    sv0=-ops.eleResponse(1,'stress')[1]
    if sv0<=0:raise RuntimeError('initial vertical stress must be compressive')
    ops.integrator('LoadControl',1.)
    rows=[]
    def record(v):
        s=ops.eleResponse(1,'stress')
        rows.append([v,s[2],1.+s[1]/sv0,s[1],ops.nodeDisp(3,1)])
    record(0.)
    for v in x[1:]:checked();record(v)
    return dict(zip(('x','y','ru_equivalent','sigma_vertical','actual_x'),np.array(rows).T))


def contact_test(c,x):
    p,t=c['params'],c['test']
    ops.model('basic','-ndm',2,'-ndf',2)
    ops.node(1,0.,0.); ops.node(2,0.,0.); ops.node(3,0.,0.)
    ops.fix(1,1,1); ops.fix(3,0,1)
    ops.uniaxialMaterial('Elastic',99,1.)
    ops.element('zeroLength',99,1,3,'-mat',99,'-dir',1)
    ops.element('zeroLengthContactASDimplex',1,1,2,p['Kn'],p['Kt'],p['mu'],
                '-orient',0.,1.,0.,'-intType',0)
    ops.timeSeries('Constant',1); ops.pattern('Plain',1,1)
    ops.load(2,0.,-t['normal_force'])
    imposed_series(2,x); ops.sp(2,1,1.)
    static_analysis()
    ops.integrator('LoadControl',0.); checked()
    ops.integrator('LoadControl',1.)
    rows=[]
    def record(v):
        ops.reactions()
        rows.append([v,ops.nodeReaction(2,1),ops.nodeReaction(1,2),ops.nodeDisp(2,2)])
    record(0.)
    for v in x[1:]: checked(); record(v)
    return dict(zip(('x','y','normal_reaction','normal_displacement'),np.array(rows).T))

def viscous_test(c,x):
    # Direct rate-controlled material test, not a displacement-only static test.
    uniaxial(1,c); ops.testUniaxialMaterial(1)
    y=[]
    for velocity in x:
        ops.setStrain(0.,float(velocity))
        y.append(ops.getStress())
    return dict(x=x,y=y)

def shell_test(c,x):
    """Unit square shell coupon; x is axial strain, y is membrane stress kPa.
    Transverse edges are free in plane, all out-of-plane DOFs restrained.
    """
    from materials import shell_section
    ops.model('basic','-ndm',3,'-ndf',6)
    for n,xy in enumerate(((0,0),(1,0),(1,1),(0,1)),1):
        ops.node(n,*xy,0.)
        ops.fix(n,1 if n in (1,4) else 0,1 if n in (1,2) else 0,1,1,1,1)
    shell_section(1,c,mat_tag=1,plate_tag=2)
    ops.element('ShellMITC4',1,1,2,3,4,1)
    imposed_series(1,x)
    for n in (2,3): ops.sp(n,1,1.)
    static_analysis()
    rows=[]
    def record(v):
        ops.reactions()
        force=sum(ops.nodeReaction(n,1) for n in (2,3))
        rows.append([v,force/c['params']['thickness'],ops.nodeDisp(3,2)])
    record(0.)
    for v in x[1:]: checked(); record(v)
    return dict(zip(('x','y','transverse_strain'),np.array(rows).T))
