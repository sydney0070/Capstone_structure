"""Force-controlled3DOF-soil/6DOF-foundation interface and free physical rotation."""
import json,sys
from pathlib import Path
import openseespy.opensees as o
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat.ssi import register_contact
from capmat.mechanics import contact_parameters
p=contact_parameters(1,30e6,12e6,1,1e5,.0001,.34)
o.wipe();o.model('basic','-ndm',3,'-ndf',3);o.node(1,0.,0.,0.);o.fix(1,1,1,1)
o.model('basic','-ndm',3,'-ndf',6);o.node(2,0.,0.,0.);o.fix(2,0,1,0,1,0,1);o.node(99,0.,0.,0.);o.fix(99,1,1,1,1,1,1)
o.uniaxialMaterial('Elastic',99,1e4);o.element('zeroLength',99,99,2,'-mat',99,'-dir',5)
register_contact(o,1,2,1,p,auxiliary_nodes=(101,102),normal_material_tag=10,rotation_material_tag=11,friction_tag=1,restore_ndf=6)
o.timeSeries('Linear',1);o.pattern('Plain',1,1);o.load(2,0.,0.,-1e5,0.,0.,0.)
o.constraints('Plain');o.numberer('RCM');o.system('BandGeneral');o.test('NormUnbalance',1e-6,100);o.algorithm('Newton');o.integrator('LoadControl',.1);o.analysis('Static');ok=o.analyze(10)
o.loadConst('-time',0.);o.timeSeries('Linear',2);o.pattern('Plain',2,2);o.load(2,10000.,0.,0.,0.,100.,0.);o.integrator('LoadControl',.1);ok+=o.analyze(10)
u=o.nodeDisp(2);expected=[10000/p['Kt_N_m'],-1e5/p['Kn_N_m'],.01]
errors=[abs(u[i]-v)/abs(v) for i,v in zip([0,2,4],expected)]
r=dict(runtime=o.version(),constraint_handler='Plain',pass_=ok==0,relative_errors=errors,expected_ux_uz_ry=expected,actual_ux_uz_ry=[u[i] for i in [0,2,4]],physical_rotation_unrestrained_by_contact=True,**{'pass':bool(ok==0 and max(errors)<1e-6)})
(ROOT/'results'/'ssi_transfer.json').write_text(json.dumps(r,indent=2));print(r)
