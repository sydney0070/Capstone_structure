"""Fixed physical hinge-length formulation; lp is an explicit integration input.
No experimental hinge-length calibration. Compare fiber/time discretization
at fixed lp; varying lp is physical-model sensitivity, NOT convergence.
"""
import sys,json,math
from pathlib import Path
import openseespy.opensees as ops
R=Path(__file__).resolve().parents[1];sys.path[:0]=[str(R),str(R/'verification')]
from capmat.rc import register_rc_section
import rc_convergence as c
OUT=R/'research'/'rc_hinge_results';OUT.mkdir(exist_ok=True);c.OUT=OUT
BARS=[(y,z,.02) for y in [-.245,.245] for z in [-.145,0,.145]]
LP=.3

def model(ny,nip,form='force',ne=1):
 ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
 meta=register_rc_section(ops,10,1,2,width_m=.4,depth_m=.6,bars=BARS,ny=ny,nz=2)
 EA=30e9*meta['concrete_area_m2']+200e9*meta['steel_area_m2']
 EI=30e9*meta['concrete_Iz_m4']+200e9*sum(math.pi*d*d/4*y*y for y,z,d in BARS)
 ops.section('Elastic',20,30e9,EA/30e9,EI/30e9)
 ops.node(1,0.,0.);ops.node(2,0.,3.);ops.fix(1,1,1,1);ops.geomTransf('PDelta',1)
 ops.beamIntegration('HingeRadau',1,10,LP,10,LP,20)
 ops.element('forceBeamColumn',1,1,2,1,1,'-iter',40,1e-12)
 ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(2,1.,0.,0.)
 ops.constraints('Plain');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',.05,100,0);ops.algorithm('Newton');ops.integrator('LoadControl',0.);ops.analysis('Static')
 return 2
c.model=model
runs=[];summary=dict(scope='physical hinge length explicit, no-overlap section, unconfined C30/SD400 assumptions, no axial preload',cases=[],comparisons=[])
for ny,steps,lp in [(24,60,.3),(36,60,.3),(48,60,.3),(36,120,.3),(36,60,.15),(36,60,.35)]:
 LP=lp;c.OUT=OUT/f'lp{lp}';c.OUT.mkdir(exist_ok=True)
 r=c.run(ny,6,steps);runs.append(r);r[0]['lp_m']=lp;summary['cases'].append(r[0])
for i,j in [(0,1),(1,2),(1,3)]:summary['comparisons'].append(c.compare(runs[i],runs[j]))
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
