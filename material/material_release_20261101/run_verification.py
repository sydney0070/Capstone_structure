"""Reproduce current release evidence. RC cyclic convergence remains a reported gate.
No installs, original-source edits, plots, or parameter refits.
"""
import hashlib,json,platform,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
TESTS=['material_patch','soil_material','ground_adapter','steel_coupon','dynamic_mass','rc_section_patch','ssi_contact','ssi_transfer','rc_convergence','rc_hinge','rc_3d','rc_continuum','laminate_bending','shell_mass','metal_shell_plastic']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.time();runs=[]
 for name in TESTS:
  p=ROOT/'verification'/f'{name}.py'
  with (ROOT/'results'/f'{name}.log').open('w') as f:job=subprocess.run([sys.executable,str(p)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
  runs.append(dict(test=name,exit_code=job.returncode,script_sha256=sha(p)))
  print(name,'exit',job.returncode,flush=True)
 evidence={}
 for name in ['material_patch','dynamic_mass','rc_section_patch','ssi_transfer','rc_3d','rc_continuum','laminate_bending','shell_mass','metal_shell_plastic']:
  d=json.loads((ROOT/'results'/f'{name}.json').read_text());evidence[name]=d.get('all_pass',d.get('pass',False))
 for name in ['soil_material','steel_coupon','ssi_contact','rc_hinge']:
  d=json.loads((ROOT/'results'/name/'summary.json').read_text());evidence[name]=d['all_pass']
 d=json.loads((ROOT/'results'/'ground_adapter.json').read_text());evidence['ground_adapter']=len(d['cases'])==9 and all(c['pass_gravity'] for c in d['cases'])
 d=json.loads((ROOT/'results'/'rc_convergence'/'summary.json').read_text());rejected_diagnostic={'arbitrary_Lobatto_RC_all_comparisons_pass':d['all_comparisons_pass'],'decision':'Rejected; replaced by fixed explicit physical HingeRadau formulation tested in rc_hinge. Not a request to change softening material constants until old arbitrary integration passes.'}
 manifest=dict(python=sys.version,executable=sys.executable,platform=platform.platform(),elapsed_s=time.time()-started,runs=runs,gates=evidence,rejected_formulations=rejected_diagnostic,all_required_gates_pass=all(evidence.values()),
  status='PRE_INTEGRATION_PACKAGE_IN_PROGRESS' if not all(evidence.values()) else 'NUMERICALLY_VERIFIED_PREINTEGRATION_WITH_EXPLICIT_LIMITS',
  code_hashes={str(p.relative_to(ROOT)):sha(p) for folder in ['capmat','verification','configs'] for p in sorted((ROOT/folder).rglob('*')) if p.is_file() and p.suffix in ['.py','.json']})
 (ROOT/'results'/'verification_manifest.json').write_text(json.dumps(manifest,indent=2))
 print(json.dumps(manifest['gates'],indent=2));return 0 if manifest['all_required_gates_pass'] and all(x['exit_code']==0 for x in runs) else 1
if __name__=='__main__':sys.exit(main())
