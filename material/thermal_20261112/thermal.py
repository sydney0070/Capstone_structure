"""Thermal strain conversion for the 19-record overlay.

Material-point histories consume total strain and temperature and pass
mechanical strain = total strain - thermal strain to the original native law.
This is NOT an implicit 3D OpenSees material plugin or a heat-transfer solver.
"""
import copy, hashlib, json, math, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'material_release_20261101'
sys.path.insert(0,str(BASE))
from capmat.materials import register_uniaxial, register_nd

def catalog():
    d=json.loads((HERE/'thermal_catalog.json').read_text())
    if hashlib.sha256(Path(d['mechanical_base_path']).read_bytes()).hexdigest()!=d['mechanical_base_sha256']:
        raise RuntimeError('Mechanical baseline changed; repeat validation before using thermal overlay')
    return d
def record(name):return copy.deepcopy(catalog()['materials'][name])
def _temperature(T):
    if not math.isfinite(T) or T<=0:raise ValueError('positive finite absolute temperature in K required')
def _length(T,rec):
    _temperature(T)
    if not 4<=T<=300:raise ValueError('NIST 6061 expansion fit valid only4–300K; no silent extrapolation')
    f=rec['length_fit']
    y=f['low_constant'] if T<18 else sum(c*T**i for i,c in enumerate(f['coefficients']))
    return 1+y*f['scale']
def alpha(name,T):
    """Instantaneous linear CTE in local material axes, not volume CTE."""
    _temperature(T);r=record(name)
    if r['model']=='nist_6061_length_curve':
        length=_length(T,r);f=r['length_fit']
        v=0. if T<18 else sum(i*c*T**(i-1) for i,c in enumerate(f['coefficients']) if i)*f['scale']/length
        return np.array([v,v,v])
    if r['model']=='constant_orthotropic':return np.array(r['alpha_principal_per_K'],float)
    return np.full(3,r['alpha_linear_per_K'])
def thermal_strain(name,T,T_ref=None,*,angle_deg=0.,allow_research_assumptions=False):
    r=record(name);_temperature(T);T_ref=r['reference_temperature_K'] if T_ref is None else T_ref;_temperature(T_ref)
    if r['model']!='nist_6061_length_curve' and not allow_research_assumptions:
        raise ValueError('Mars-temperature use is a declared research assumption; pass allow_research_assumptions=True explicitly')
    lo,hi=r['research_temperature_range_K']
    if not (lo<=T<=hi and lo<=T_ref<=hi) and r['model']!='nist_6061_length_curve':
        raise ValueError(f'Configured research envelope{lo}–{hi}K exceeded; supply new evidence/configuration')
    if not math.isfinite(angle_deg):raise ValueError('finite ply angle required')
    a=alpha(name,T)*(T-T_ref) if r['model']!='nist_6061_length_curve' else np.full(3,_length(T,r)/_length(T_ref,r)-1.)
    c=math.cos(math.radians(angle_deg));s=math.sin(math.radians(angle_deg))
    return np.array([a[0]*c*c+a[1]*s*s,a[0]*s*s+a[1]*c*c,a[2],2*(a[0]-a[1])*c*s,0.,0.])
def mechanical_strain(name,total_strain,T,T_ref=None,**options):
    e=np.asarray(total_strain,dtype=float)
    if e.shape!=(6,) or not np.isfinite(e).all():raise ValueError('six finite total engineering strains required')
    return e-thermal_strain(name,T,T_ref,**options)
def mass_only_state(name,T,*,mass_kg,volume_ref_m3,T_ref=None,allow_research_assumptions=False):
    if not record(name).get('mechanical_thermal_stress_applicable') is False:raise ValueError('MassOnly record required')
    if not math.isfinite(mass_kg) or not math.isfinite(volume_ref_m3) or min(mass_kg,volume_ref_m3)<=0:raise ValueError('positive mass and reference volume required')
    e=thermal_strain(name,T,T_ref,allow_research_assumptions=allow_research_assumptions)
    V=volume_ref_m3*np.prod(1+e[:3])
    return dict(mass_kg=mass_kg,assumed_volume_m3=float(V),assumed_density_kg_m3=float(mass_kg/V),thermal_force_N=0.,structural_stiffness_added=False,bulk_volume_change_is_assumption=True)

def register_uniaxial_fixed_temperature(ops,name,tag,base_tag,T,*,T_ref=None,allow_research_assumptions=False):
    """Fixed-temperature native wrapper, only for the supported uniaxial laws.

    InitStrainMaterial commits its initial state on creation/update in engine3.4.
    Hence plastic thermal cycles/free-expansion paths must use the explicit
    material-point driver until an equilibrium-consistent FE adapter is validated.
    """
    if tag==base_tag:raise ValueError('distinct material tags required')
    m=json.loads((BASE/'configs/catalog.json').read_text())['materials'][name]
    if m['law'] not in ('Steel02','Concrete02','J2PerfectPlastic','ElasticIsotropic'):raise ValueError('Unsupported uniaxial law')
    eps=thermal_strain(name,T,T_ref,allow_research_assumptions=allow_research_assumptions)[0]
    # Constructor precommit may activate plasticity before free expansion. Refuse
    # nonzero thermal offsets on nonlinear laws rather than quietly preload them.
    if abs(eps)>1e-15 and m['law']!='ElasticIsotropic':
        raise NotImplementedError('Engine3.4 InitStrain precommits nonlinear offsets; use material-point driver or validated modern FE coupling')
    register_uniaxial(ops,name,base_tag)
    ops.uniaxialMaterial('InitStrainMaterial',tag,base_tag,-float(eps))
    return dict(material_tag=tag,thermal_strain=float(eps),scope='fixed-temperature uniaxial wrapper only')

def require_native_3d_thermal_support(ops):
    """Fail explicitly on installed3.4; never claim metadata enables FE coupling."""
    version=tuple(int(x) for x in ops.version().split('.')[:2])
    if version<(3,7):raise NotImplementedError('Installed OpenSees '+ops.version()+' lacks the validated3D InitStrain adapter; full SSI thermal coupling is not implemented. Material-point history is available.')
    raise NotImplementedError('Newer engine requires project benchmarks and adapter validation before full SSI thermal coupling')

def laminate_thermal_resultants(angles_deg,thicknesses_m,temperatures_K,T_ref=None):
    from capmat.mechanics import qbar,reduced_Q
    m=json.loads((BASE/'configs/catalog.json').read_text())['materials']['cfrp_as4_8552_ud_rtd']
    a=np.asarray(angles_deg,float);h=np.asarray(thicknesses_m,float);T=np.asarray(temperatures_K,float)
    if not(len(a)==len(h)==len(T)) or not len(a) or not np.isfinite(h).all() or min(h)<=0:raise ValueError('matching finite positive ply thicknesses required')
    z=np.r_[-h.sum()/2,-h.sum()/2+np.cumsum(h)];N=np.zeros(3);M=np.zeros(3)
    for i,angle in enumerate(a):
        et=thermal_strain('cfrp_as4_8552_ud_rtd',float(T[i]),T_ref,angle_deg=float(angle),allow_research_assumptions=True)[[0,1,3]]
        sigma=qbar(reduced_Q(m),float(angle))@et
        N+=sigma*h[i];M+=sigma*(z[i+1]**2-z[i]**2)/2
    return dict(N_thermal_N_m=N,M_thermal_N=M,scope='CLT thermal resultants, not automatically applied to native shell elements')
