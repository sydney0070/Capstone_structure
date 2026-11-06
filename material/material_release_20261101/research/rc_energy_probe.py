"""Exploratory energy regularization, not adopted unless full-path checks pass.
Reference Gft450J/m² and compression softening work10800J/m² are derived
from ORIGINAL Concrete02 and reference end quadrature length0.15m.
They are NOT measured fracture energies. No numeric fitted constants.
"""
import sys,math,json
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
from numpy.polynomial.legendre import Legendre
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'verification'));import rc_convergence as c
OUT=R/'research'/'rc_energy_results';OUT.mkdir(exist_ok=True);c.OUT=OUT

def model(ny,nip,form='force',ne=1):
 ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
 ops.uniaxialMaterial('Steel02',2,400e6,200e9,.01,18,.925,.15)
 poly=Legendre.basis(nip-1);x=np.r_[-1.,poly.deriv().roots(),1.];weights=1/(nip*(nip-1)*poly(x)**2) # [0,1]
 sections=[]
 for i,w in enumerate(weights):
  lch=3/ne*w;tag=100+i
  ops.uniaxialMaterial('Concrete02',tag,-30e6,-.002,-6e6,-.002-2*10800/(36e6*lch),.1,3e6,9e12*lch/(2*450))
  ops.section('Fiber',tag);ops.patch('rect',tag,ny,4,-.3,-.2,.3,.2)
  for y in [-.245,.245]:ops.layer('straight',2,3,math.pi*.02**2/4,y,-.145,y,.145)
  sections.append(tag)
 for j in range(ne+1):ops.node(j+1,0.,3*j/ne)
 ops.fix(1,1,1,1);ops.geomTransf('PDelta',1);ops.beamIntegration('Lobatto',1,nip,*sections)
 for j in range(ne):ops.element('forceBeamColumn',j+1,j+1,j+2,1,1,'-iter',30,1e-12)
 top=ne+1;ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(top,1.,0.,0.)
 ops.constraints('Plain');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',.05,100,0);ops.algorithm('Newton');ops.integrator('LoadControl',0.);ops.analysis('Static')
 return top
c.model=model
runs=[];summary=dict(hypothesis='energy per integration length regularization; reference energies derived from old baseline, not calibrated',cases=[],comparisons=[])
for args in [(36,5,60),(48,5,60),(48,7,60),(48,9,60),(48,10,60),(48,5,120),(48,7,120)]:
 r=c.run(*args);runs.append(r);summary['cases'].append(r[0])
for i,j in [(0,1),(1,2),(2,3),(3,4),(1,5),(2,6)]:summary['comparisons'].append(c.compare(runs[i],runs[j]))
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
