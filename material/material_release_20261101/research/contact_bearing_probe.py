import sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'verification'));import ssi_contact as t
import openseespy.opensees as ops
src=Path(t.__file__).read_text();start=src.index('def case(');end=src.index('def viscous()',start);code=src[start:end]
code=code.replace("'-ndf',3", "'-ndf',6").replace('ops.fix(1,1,1,1);ops.fix(99,0,1,1)', 'ops.fix(1,1,1,1,1,1,1);ops.fix(99,0,1,1,1,1,1)')
code=code.replace("ops.fix(2,*[int(d not in [axis,tangent]) for d in [1,2,3]])", "ops.fix(2,*[int(d not in [axis,tangent]) for d in [1,2,3]],1,1,1)")
code=code.replace("register_contact(ops,1,2,1,p,normal_axis=axis)", """ops.frictionModel('Coulomb',1,p['mu'])
 ops.uniaxialMaterial('Elastic',10,p['Kn_N_m']);ops.uniaxialMaterial('Elastic',11,0.)
 normal=[0.,0.,0.];normal[axis-1]=1.;direction=[0.,0.,0.];direction[tangent-1]=1.
 ops.element('flatSliderBearing',1,1,2,1,p['Kt_N_m'],'-P',10,'-T',11,'-My',11,'-Mz',11,'-orient',*normal,*direction,'-iter',100,1e-10)""")
code=code.replace('np.linalg.norm(r+g)', 'np.linalg.norm((r+g)[:3])')
t.OUT=R/'research'/'contact_bearing_results';t.OUT.mkdir(exist_ok=True)
exec(code,t.__dict__)
for a in [1,2,3]:print(t.case(1,0,a))
