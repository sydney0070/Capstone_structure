"""RC Step 3: zeroLengthSection moment-curvature + one forceBeamColumn cantilever.
These are illustrative numerical verification, NOT calibrated performance tests.
"""
import csv, json, math
from dataclasses import asdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openseespy.opensees as ops
from rc_section import create_rc_fiber_section, SectionGeometry
from rc_materials import ConcreteParams, SteelParams

OUT=Path(__file__).resolve().parent/'results_rc_step3'
OUT.mkdir(exist_ok=True)
GEOM=SectionGeometry()


def solver():
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormUnbalance',1e-5,60)
    ops.algorithm('Newton')
    ops.analysis('Static')


def step_safe(disp_node, disp_dof, inc):
    # Preserve path direction; retry with a smaller displacement increment.
    for div in (1,2,5,10):
        ops.integrator('DisplacementControl',disp_node,disp_dof,inc/div)
        ok=True
        for _ in range(div):
            if ops.analyze(1)!=0:
                ok=False
                break
        if ok:
            return True
        # A partial increment may have committed; abort rather than silently corrupting the path.
        raise RuntimeError('Nonconvergence; partial increment possible. Review log / reduce global step.')
    return False


def write_csv(path, header, rows):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)


def plot_xy(path, rows, xlabel, ylabel, title):
    a=np.asarray(rows,float)
    plt.figure(figsize=(7,4.5)); plt.plot(a[:,0],a[:,1],lw=1.5)
    plt.xlabel(xlabel);plt.ylabel(ylabel);plt.title(title);plt.grid(True,alpha=.3)
    plt.tight_layout();plt.savefig(path,dpi=160);plt.close()


def moment_curvature(axial_load=-0.05*30e6*GEOM.width*GEOM.depth, max_k=0.025, steps=150):
    ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
    geom=create_rc_fiber_section(ops)
    ops.node(1,0,0);ops.node(2,0,0)
    ops.fix(1,1,1,1);ops.fix(2,0,1,0)
    ops.element('zeroLengthSection',1,1,2,10)
    ops.timeSeries('Constant',1);ops.pattern('Plain',1,1)
    ops.load(2,axial_load,0,0)
    solver();ops.integrator('LoadControl',0.0)
    if ops.analyze(1)!=0: raise RuntimeError('Axial preload failed')
    ops.loadConst('-time',0.0)
    ops.timeSeries('Linear',2);ops.pattern('Plain',2,2)
    ops.load(2,0,0,1.0)
    rows=[]
    for i in range(steps):
        ops.integrator('DisplacementControl',2,3,max_k/steps)
        if ops.analyze(1)!=0: raise RuntimeError(f'M-K nonconvergence at step {i}')
        curvature=ops.nodeDisp(2,3)
        moment=ops.getLoadFactor(2) # reference 1 N*m moment: multiplier is moment in N*m
        rows.append((curvature,moment/1e3))
    write_csv(OUT/'moment_curvature.csv',['curvature_1_per_m','moment_kNm'],rows)
    plot_xy(OUT/'moment_curvature.png',rows,'Curvature (1/m)','Moment (kN m)','RC section moment-curvature')
    return {'completed_steps':len(rows),'max_abs_moment_kNm':max(abs(v) for _,v in rows),'axial_load_N':axial_load}


def beam_cyclic(height=3.0, amplitude=0.018, cycles=(0.003,-0.003,0.006,-0.006,0.012,-0.012,0.018,-0.018,0.0),substeps=30):
    ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
    create_rc_fiber_section(ops)
    ops.node(1,0,0);ops.node(2,0,height)
    ops.fix(1,1,1,1)
    ops.geomTransf('PDelta',1)
    ops.beamIntegration('Lobatto',1,10,5)
    ops.element('forceBeamColumn',1,1,2,1,1)
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
    ops.load(2,1.0,0,0)  # reference horizontal load for displacement control
    solver()
    u=0.0;rows=[(0.0,0.0)]
    for target in cycles:
        du=(target-u)/substeps
        for j in range(substeps):
            ops.integrator('DisplacementControl',2,1,du)
            if ops.analyze(1)!=0:
                raise RuntimeError(f'Beam nonconvergence: target={target}, substep={j}. Prior points saved only in memory.')
            u=float(ops.nodeDisp(2,1))
            shear=float(ops.getLoadFactor(1))
            rows.append((u,shear/1e3))
    write_csv(OUT/'beam_cyclic.csv',['tip_displacement_m','applied_horizontal_force_kN'],rows)
    plot_xy(OUT/'beam_cyclic.png',rows,'Tip displacement (m)','Lateral force (kN)','RC cantilever cyclic lateral response')
    return {'completed_points':len(rows),'max_abs_force_kN':max(abs(f) for _,f in rows),'max_abs_disp_m':max(abs(v) for v,_ in rows),'height_m':height}


def main():
    results={'status':'PROTOTYPE_NOT_EXPERIMENTALLY_CALIBRATED','units':'N,m,Pa','geometry_provisional':asdict(GEOM),'concrete_provisional':asdict(ConcreteParams()),'steel_provisional':asdict(SteelParams()),'documentation':['https://openseespydoc.readthedocs.io/en/latest/src/MomentCurvature.html','https://opensees.github.io/OpenSeesDocumentation/user/manual/model/beamIntegrations/Lobatto.html']}
    for label,func in [('moment_curvature',moment_curvature),('beam_cyclic',beam_cyclic)]:
        try:
            results[label]={'status':'PASS_NUMERICAL_EXECUTION',**func()}
        except Exception as e:
            results[label]={'status':'FAILED','error':str(e)}
        print(label,results[label])
    (OUT/'summary.json').write_text(json.dumps(results,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Saved:',OUT)
    if any(results[x]['status']=='FAILED' for x in ('moment_curvature','beam_cyclic')):
        raise SystemExit('One or more tests failed; inspect summary.json')

if __name__=='__main__': main()
