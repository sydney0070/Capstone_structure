"""Analytical patch checks, units, mass metadata, and principal-axis laminate.
Finite reference tests, not source-data calibration.
"""
import json,sys,math,traceback
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from capmat import material,register_uniaxial,register_plate,laminate_abd,shielding,contact_parameters,radiation_dashpots
from capmat.laminate import register_crossply
results=[]
def save():
 (ROOT/'results'/'material_patch.json').write_text(json.dumps(dict(runtime=ops.version(),cases=results,all_pass=all(r['pass'] for r in results)),indent=2))
def add(name,fn):
 try:
  r=fn();results.append(dict(name=name,**r));print(name,r,flush=True)
 except Exception as e:results.append(dict(name=name,pass_=False,error=str(e),traceback=traceback.format_exc()));results[-1]['pass']=False
 save()
def axial(name):
 ops.wipe();meta=register_uniaxial(ops,name,1);ops.testUniaxialMaterial(1)
 m=material(name);E=m.get('E_Pa',m.get('params',{}).get('E'));eps=1e-6
 ops.setStrain(eps);s=ops.getStress();err=abs(s-E*eps)/(E*eps)
 info=dict(relative_initial_tangent_error=err,stress_Pa=s,metadata=meta)
 if m['law']=='J2PerfectPlastic':
  ops.setStrain(.01);stress=ops.getStress();info['stress_at_one_percent_Pa']=stress;info['yield_relative_error']=abs(stress-m['fy_Pa'])/m['fy_Pa'];err=max(err,info['yield_relative_error'])
 return {**info,'pass':err<1e-3}
def shell(name,elastic=False,angles=None):
 ops.wipe();ops.model('basic','-ndm',3,'-ndf',6)
 if angles:
  meta=register_crossply(ops,10,100,angles);lam=laminate_abd(angles);A=lam['A_N_m'];eff=A[0,0]-A[0,1]**2/A[1,1];nu=A[0,1]/A[1,1]
 else:
  meta=register_plate(ops,name,100,101,10,.01,elastic=elastic);m=material(name);eff=m['E_Pa']*.01;nu=m['nu']
 for i,(x,y) in enumerate([(0,0),(1,0),(1,1),(0,1)],1):ops.node(i,x,y,0);ops.fix(i,int(x==0),int(y==0),1,1,1,1)
 ops.element('ShellMITC4',1,1,2,3,4,10);ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
 eps=1e-5
 for n in [2,3]:ops.sp(n,1,eps)
 ops.constraints('Transformation');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',1e-7,50);ops.algorithm('Newton');ops.integrator('LoadControl',1);ops.analysis('Static');ok=ops.analyze(1)
 ops.reactions();force=-sum(ops.nodeReaction(n,1) for n in [1,4]);uy=ops.nodeDisp(3,2)
 ef=abs(force-eff*eps)/max(abs(eff*eps),1e-9);en=abs(uy+nu*eps)/max(abs(nu*eps),1e-12)
 return dict(pass_=ok==0,**{'pass':bool(ok==0 and max(ef,en)<1e-5)},force_N=force,expected_N=eff*eps,force_relative_error=ef,poisson_relative_error=en,metadata=meta)
def mechanics():
 d=shielding(2,10,1.62,density_kg_m3=1700);p=contact_parameters(2,30e6,12e6,1,1e5,.0001,.3);r=radiation_dashpots(1800,30e6,12e6,2)
 sym=laminate_abd([0,90,90,0]);balanced=laminate_abd([45,-45,-45,45]);iso={'E1_Pa':100.,'E2_Pa':100.,'nu12':.25,'G12_Pa':40.}
 from capmat.mechanics import qbar,reduced_Q
 err=float(np.max(abs(qbar(reduced_Q(iso),37)-reduced_Q(iso))))
 passed=d['mass_kg']==34000 and d['weight_N']==55080 and p['estimated_penetration_m']<=.0001 and np.max(abs(sym['B_N']))<1e-8 and np.max(abs(balanced['A_N_m'][[0,1],2]))<1e-6 and err<1e-10
 return {'pass':bool(passed),'shield':d,'contact':p,'dashpots':r,'isotropic_rotation_absolute_error':err,'symmetric_laminate_B_max':float(np.max(abs(sym['B_N'])))}
if __name__=='__main__':
 for n in ['steel_s355_coupon','steel_s355_nominal','al6061_t6_extrusion','rebar_sd400_prototype','sintered_jsc1a_high_porosity']:add('uniaxial_'+n,lambda n=n:axial(n))
 for n in ['steel_s355_nominal','al6061_t6_extrusion','sintered_jsc1a_high_porosity']:add('shell_'+n,lambda n=n:shell(n))
 for a in [[0,0,0,0],[90,90,90,90],[0,90,90,0]]:add('crossply_'+str(a),lambda a=a:shell('cfrp',angles=a))
 add('independent_conversions_and_CLT',mechanics)
