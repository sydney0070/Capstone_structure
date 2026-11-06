"""Geometry, elastic beam and axial response of the production no-overlap section."""
import json,math,sys,traceback
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.rc import rectangle_fibers,register_rc_section
BARS=[(y,z,.02) for y in [-.245,.245] for z in [-.145,0,.145]]
def check(ny):
 fibers,meta=rectangle_fibers(.4,.6,BARS,ny=ny,nz=12)
 Ac=.24-len(BARS)*math.pi*.02**2/4;Es=200e9;Ec=30e9
 Iz_gross=.4*.6**3/12;Iz_holes=sum(math.pi*.02**4/64+math.pi*.02**2/4*y*y for y,z,d in BARS);Iz_steel=sum(math.pi*.02**2/4*y*y for y,z,d in BARS)
 EI=Ec*(Iz_gross-Iz_holes)+Es*Iz_steel;EA=Ec*Ac+Es*(.24-Ac)
 ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
 data=register_rc_section(ops,10,1,2,width_m=.4,depth_m=.6,bars=BARS,ny=ny,nz=12,elastic=True)
 ops.node(1,0.,0.);ops.node(2,3.,0.);ops.fix(1,1,1,1);ops.geomTransf('Linear',1);ops.beamIntegration('Lobatto',1,10,5);ops.element('forceBeamColumn',1,1,2,1,1)
 ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(2,10000.,10000.,0.)
 ops.constraints('Plain');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',1e-5,40);ops.algorithm('Newton');ops.integrator('LoadControl',1.);ops.analysis('Static');ok=ops.analyze(1)
 ux,uy,rz=ops.nodeDisp(2);ex=10000*3/EA;ey=10000*3**3/(3*EI)
 e1=abs(ux-ex)/ex;e2=abs(uy-ey)/ey
 return dict(ny=ny,metadata=data,axial_relative_error=e1,bending_relative_error=e2,expected_deflections_m=[ex,ey],computed_deflections_m=[ux,uy],**{'pass':bool(ok==0 and max(e1,e2)<1e-5)})
if __name__=='__main__':
 r=dict(runtime=ops.version(),cases=[])
 for n in [24,36,48]:
  try:d=check(n)
  except Exception as e:d=dict(ny=n,error=str(e),traceback=traceback.format_exc(),**{'pass':False})
  r['cases'].append(d);print(d,flush=True)
 r['all_pass']=all(d['pass'] for d in r['cases']);(ROOT/'results'/'rc_section_patch.json').write_text(json.dumps(r,indent=2))
