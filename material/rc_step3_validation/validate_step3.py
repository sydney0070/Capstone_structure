"""RC Step 3 follow-up: axial equilibrium, fiber tracking, cyclic step sensitivity.
Requires openseespy numpy matplotlib. Prototype-only C30/SD400 material parameters.
Run from this folder. Does not overwrite the original Step 3 results.
"""
from pathlib import Path
import csv, json, math, traceback
import numpy as np
import openseespy.opensees as ops
from rc_section import create_rc_fiber_section, SectionGeometry
import test_rc_step3 as base

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results_rc_validation'
OUT.mkdir(exist_ok=True)
GEOM=SectionGeometry()


def add_fiber_recorders(mode, folder, mat_tags=(1,2)):
    """Concrete fibers near outer faces and rebar at top/bottom, at integration point 1 for beams."""
    folder.mkdir(parents=True,exist_ok=True)
    y_conc=GEOM.depth/2 - GEOM.depth/(2*GEOM.ny)
    y_steel=GEOM.depth/2 - GEOM.cover - GEOM.bar_diameter/2
    for name,y,mt in [('concrete_top',y_conc,1),('concrete_bottom',-y_conc,1),
                      ('steel_top',y_steel,2),('steel_bottom',-y_steel,2)]:
        args=(['section'] if mode=='section' else ['section',1]) + ['fiber',float(y),0.0,int(mt),'stressStrain']
        ops.recorder('Element','-file',str(folder/(name+'.txt')),'-ele',1,*args)


def read_fibers(folder):
    result={}
    for p in folder.glob('*.txt'):
        if p.stat().st_size==0:
            result[p.stem]={'status':'EMPTY'}; continue
        try:
            a=np.loadtxt(p,ndmin=2)
            if a.shape[1]<2: raise ValueError('expected stress and strain columns')
            result[p.stem]={'status':'RECORDED','points':len(a),'min_stress_MPa':float(a[:,0].min()/1e6),
                            'max_stress_MPa':float(a[:,0].max()/1e6),
                            'min_strain':float(a[:,1].min()),'max_strain':float(a[:,1].max())}
        except Exception as exc:result[p.stem]={'status':'UNREADABLE','error':str(exc)}
    return result


def axial_equilibrium():
    ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
    create_rc_fiber_section(ops)
    ops.node(1,0,0);ops.node(2,0,0)
    ops.fix(1,1,1,1);ops.fix(2,0,1,0)
    ops.element('zeroLengthSection',1,1,2,10)
    fdir=OUT/'axial_fibers';add_fiber_recorders('section',fdir)
    target=-360000.0
    ops.timeSeries('Constant',1);ops.pattern('Plain',1,1)
    ops.load(2,target,0,0)
    base.solver()
    ops.integrator('LoadControl',1.0)
    ok=ops.analyze(1)
    if ok!=0:raise RuntimeError('axial preload failed')
    ops.reactions()
    support_rx=float(ops.nodeReaction(1,1))
    balance=abs(support_rx+target)
    ops.remove('recorders')
    return {'load_N':target,'support_reaction_N':support_rx,
            'absolute_force_imbalance_N':balance,
            'relative_force_imbalance':balance/abs(target),
            'equilibrium_pass_1e-4':balance/abs(target)<1e-4,
            'fiber_results':read_fibers(fdir)}


def beam_with_fiber_recorders(substeps=60, n_ip=5, name='baseline'):
    ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
    create_rc_fiber_section(ops)
    ops.node(1,0,0);ops.node(2,0,3.0);ops.fix(1,1,1,1)
    ops.geomTransf('PDelta',1)
    ops.beamIntegration('Lobatto',1,10,n_ip)
    ops.element('forceBeamColumn',1,1,2,1,1)
    fdir=OUT/name/'fibers'; add_fiber_recorders('beam',fdir)
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
    ops.load(2,1,0,0)
    base.solver()
    targets=(.003,-.003,.006,-.006,.012,-.012,.018,-.018,0.)
    rows=[(0.,0.)]; u=0.
    for target in targets:
        du=(target-u)/substeps
        for j in range(substeps):
            ops.integrator('DisplacementControl',2,1,du)
            if ops.analyze(1)!=0:
                base.write_csv(OUT/name/'partial.csv', ['tip_displacement_m','applied_horizontal_force_kN'],rows)
                raise RuntimeError(f'{name}: nonconvergence target={target}, step={j}')
            u=float(ops.nodeDisp(2,1)); rows.append((u,float(ops.getLoadFactor(1))/1000.))
    ops.reactions()
    reaction_x=float(ops.nodeReaction(1,1))
    reaction_m=float(ops.nodeReaction(1,3))
    ops.remove('recorders')
    base.write_csv(OUT/name/'beam_cyclic.csv',['tip_displacement_m','applied_horizontal_force_kN'],rows)
    a=np.array(rows)
    integrated=float(np.trapezoid(a[:,1]*1000,a[:,0])) if hasattr(np,'trapezoid') else float(np.trapz(a[:,1]*1000,a[:,0]))
    return {'points':len(rows),'max_abs_force_kN':float(np.max(np.abs(a[:,1]))),
            'endpoint_force_kN':float(a[-1,1]), 'endpoint_displacement_m':float(a[-1,0]),
            'support_reaction_x_N':reaction_x,'support_reaction_moment_Nm':reaction_m,
            'horizontal_equilibrium_rel':abs(reaction_x+a[-1,1]*1000)/max(1,abs(a[-1,1]*1000)),
            'path_integral_J':integrated,'fiber_results':read_fibers(fdir)}


def compare_histories(x, y):
    """Compare force at matching targets, use interpolation over progression fraction (not on nonmonotonic displacement)."""
    n=min(len(x),len(y))
    a=np.asarray(x,float);b=np.asarray(y,float)
    ta=np.linspace(0,1,len(a));tb=np.linspace(0,1,len(b))
    q=np.linspace(0,1,500)
    fa=np.interp(q,ta,a[:,1]);fb=np.interp(q,tb,b[:,1])
    return {'rms_force_difference_kN':float(np.sqrt(np.mean((fa-fb)**2))),
            'max_force_difference_kN':float(np.max(abs(fa-fb)))}


def main():
    result={'status':'PROTOTYPE_NUMERICAL_VERIFICATION_NOT_CALIBRATION','tests':{},'warnings':[]}
    try: result['tests']['axial_equilibrium']=axial_equilibrium()
    except Exception as exc:result['tests']['axial_equilibrium']={'status':'FAILED','error':str(exc),'traceback':traceback.format_exc()}
    for n,sub in [('increment_60',60),('increment_120',120)]:
        try:
            result['tests'][n]={'status':'COMPLETED',**beam_with_fiber_recorders(substeps=sub,name=n)}
        except Exception as exc:result['tests'][n]={'status':'FAILED','error':str(exc),'traceback':traceback.format_exc()}
    pa=OUT/'increment_60'/'beam_cyclic.csv';pb=OUT/'increment_120'/'beam_cyclic.csv'
    if pa.exists() and pb.exists():
        a=np.genfromtxt(pa,delimiter=',',skip_header=1);b=np.genfromtxt(pb,delimiter=',',skip_header=1)
        result['tests']['increment_sensitivity']=compare_histories(a,b)
        result['tests']['increment_sensitivity']['peak_force_rel_difference']=abs(np.max(np.abs(a[:,1]))-np.max(np.abs(b[:,1])))/max(1,np.max(np.abs(b[:,1])))
    result['warnings'] += [
      'Applied load at tip is a reference unit load; reported load factor has units N.',
      'At end of displacement-controlled cycles, u=0 does not mean F=0. Residual displacement requires load-controlled unloading to zero force.',
      'Fiber recorder covers selected concrete/steel fibers only, not all fibers.',
      'Directly comparing path-integral values for open paths does not prove dissipated energy.',
      'Increment sensitivity on 60 vs 120 subdivisions does not establish mesh/objectivity or experimental calibration.']
    (OUT/'validation_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result['tests'],ensure_ascii=False,indent=2))
    print('Saved:',OUT)

if __name__=='__main__':main()
