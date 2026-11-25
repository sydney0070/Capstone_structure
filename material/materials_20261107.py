"""2026-11-07 catalogue audit; no refitting and no edits to released coefficients."""
from pathlib import Path
import copy,hashlib,json,math
import numpy as np
ROOT=Path(__file__).resolve().parent
RELEASE=ROOT/'material_release_20261101'/'configs'/'catalog.json'
GROUND=ROOT.parent/'ground'/'ground_modeling'/'ground_hybrid'/'configs'/'materials.json'

def validate_record(name,m):
    def positive(k):
        if k not in m or not math.isfinite(float(m[k])) or float(m[k])<=0:raise ValueError(f'{name}: positive finite {k} required')
    positive('density_kg_m3')
    law=m['law']
    if law not in {'Steel02','Concrete02','DruckerPrager','ElasticIsotropic','ElasticOrthotropic','J2PerfectPlastic','MassOnly','PressureDependMultiYield02'}:
        raise ValueError(f'{name}: unsupported law {law}')
    required={'Steel02':{'E','Fy','R0','b','cR1','cR2'},'Concrete02':{'fpc','epsc0','fpcu','epsU','lam','ft','Ets'}}
    if law in required and (not isinstance(m.get('params'),dict) or required[law]-set(m['params'])):
        raise ValueError(f'{name}: missing required {law} parameters')
    if 'nu' in m and not (-1<float(m['nu'])<.5):raise ValueError(f'{name}: invalid nu')
    if law=='Steel02':
        p=m['params']
        if not all(math.isfinite(float(v)) for v in p.values()) or min(p['E'],p['Fy'],p['R0'])<=0 or not 0<=p['b']<1:raise ValueError(f'{name}: invalid Steel02 parameters')
    if law=='J2PerfectPlastic':positive('fy_Pa')
    if law=='ElasticOrthotropic':
        for key in ['E1_Pa','E2_Pa','G12_Pa']:positive(key)
        E1,E2=m['E1_Pa'],m['E2_Pa'];n12,n23=m['nu12'],m['nu23']
        S=np.array([[1/E1,-n12/E1,-n12/E1],[-n12/E1,1/E2,-n23/E2],[-n12/E1,-n23/E2,1/E2]])
        if not np.isfinite(S).all() or min(np.linalg.eigvalsh(S))<=0:raise ValueError(f'{name}: non-positive orthotropic compliance')
    if law=='PressureDependMultiYield02':
        for key in ['G_Pa','K_Pa','p_ref_Pa','peak_shear_strain']:positive(key)
        for key in ['contrac1','contrac2','contrac3','dilat1','dilat2','dilat3','numerical_c_Pa']:
            if not math.isfinite(m[key]) or m[key]<0:raise ValueError(f'{name}: invalid {key}')
        if not 0<m['phi_deg']<60 or m['n_surfaces']<2:raise ValueError(f'{name}: invalid PDMY yield surfaces/angle')
    if law=='Concrete02':
        p=m['params']
        if not all(math.isfinite(float(v)) for v in p.values()) or not p['fpc']<p['fpcu']<0 or not p['epsU']<p['epsc0']<0 or min(p['ft'],p['Ets'])<=0:raise ValueError(f'{name}: invalid Concrete02')
    if 'E_Pa' in m:
        positive('E_Pa')
        if 'G_Pa' in m and 'K_Pa' in m:
            positive('G_Pa');positive('K_Pa')
            G,K=m['G_Pa'],m['K_Pa'];E=9*K*G/(3*K+G);nu=(3*K-2*G)/(2*(3*K+G))
            if not math.isclose(E,m['E_Pa'],rel_tol=1e-8) or not math.isclose(nu,m['nu'],abs_tol=1e-8):raise ValueError(f'{name}: inconsistent E/nu/G/K')
    if law=='DruckerPrager':
        if not all(math.isfinite(float(m[k])) for k in ['phi_deg','c_Pa','psi_deg']):raise ValueError(f'{name}: nonfinite DP parameters')
        if not 0<m['phi_deg']<89 or m['c_Pa']<0 or not 0<=m['psi_deg']<=m['phi_deg']:raise ValueError(f'{name}: invalid DP strength/flow')
    return copy.deepcopy(m)

def compare_parameters(name, a, b, path=''):
    """Compare shared constitutive fields recursively, retaining metadata freedom."""
    if isinstance(a, dict) and isinstance(b, dict):
        keys = set(a) | set(b) if path == 'params' or path.startswith('params.') else set(a) & set(b)
        for key in keys:
            field = f'{path}.{key}' if path else key
            if key not in a or key not in b:
                raise ValueError(f'Conflicting duplicate catalogue parameter: {name}/{field}')
            compare_parameters(name, a[key], b[key], field)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10):
            raise ValueError(f'Conflicting duplicate catalogue value: {name}/{path}')
    elif path == 'law' or path == 'params' or path.startswith('params.'):
        if a != b:
            raise ValueError(f'Conflicting duplicate catalogue parameter: {name}/{path}')


def audit_catalogs():
    merged={};sources=[];duplicates=[];all_records=[]
    for path in [RELEASE,GROUND]:
        raw=path.read_bytes();d=json.loads(raw)
        sources.append(dict(path=str(path),sha256=hashlib.sha256(raw).hexdigest(),records=len(d['materials'])))
        for n,m in d['materials'].items():
            validate_record(n,m)
            all_records.append(m)
            if n in merged:
                compare_parameters(n,merged[n],m)
                duplicates.append(n)
            merged[n]=m
    source_keys=set()
    for path in [RELEASE,GROUND]:source_keys.update(json.loads(path.read_text()).get('sources',{}))
    missing={key for m in all_records for key in m.get('sources',[]) if key not in source_keys}
    if missing:raise ValueError(f'Unresolved source identifiers: {missing}')
    return dict(sources=sources,checked_records=len(merged),duplicate_names=duplicates,
        units='m,N,s,kg,Pa',restrictions={n:dict(status=m.get('status'),assumptions=m.get('assumptions',[])) for n,m in merged.items()},
        conclusion='Numerical range/consistency checks only; not experimental calibration or design certification')

if __name__=='__main__':print(json.dumps(audit_catalogs(),indent=2))
