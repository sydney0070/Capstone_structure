"""kN-m-s-tonne-kPa. Values are demonstration seeds, NOT calibrated properties."""
from copy import deepcopy

CASES = {
    'elastic_soil': dict(kind='soil3d', material='ElasticIsotropic',
        params=dict(E=180000., nu=.3, rho=1.9),
        bounds=dict(E=[90000., 300000.])),
    'sand_pdmy02': dict(kind='soil3d', material='PressureDependMultiYield02',
        params=dict(rho=1.9, G=100000., K=233000., phi=33.5, gamma_peak=.1,
            p_ref=101., d=.5, phi_pt=25.5, contraction1=.045, contraction3=.15,
            dilation1=.06, dilation3=.15, surfaces=20, contraction2=5.,
            dilation2=3., liquefaction1=1., liquefaction2=0., e=.7,
            cs1=.9, cs2=.02, cs3=.7, pa=101.),
        bounds=dict(G=[60000., 140000.], phi=[30., 38.], contraction1=[.005,.10])),
    'clay_pimy': dict(kind='soil3d', material='PressureIndependMultiYield',
        params=dict(rho=1.5, G=60000., K=300000., cohesion=37., gamma_peak=.1),
        bounds=dict(G=[20000.,100000.], cohesion=[15.,70.])),
    'sand_pm4': dict(kind='soil2dcv', material='PM4Sand',
        params=dict(Dr=.35, G0=476., hpo=.53, rho=1.9, pa=101.3,
            emax=.8, emin=.5, nu=.333333333333),
        bounds=dict(G0=[300.,700.], hpo=[.2,1.])),
    'silt_pm4': dict(kind='soil2dcv', material='PM4Silt',
        params=dict(Su=25., Su_ratio=.25, G0=500., hpo=.5, rho=1.8,
            Su_factor=1., pa=101.3, nu=.333333333333, e=.9),
        bounds=dict(Su=[15.,40.], G0=[250.,800.], hpo=[.2,1.])),
    'concrete_cover': dict(kind='uniaxial', material='Concrete02',
        params=dict(fpc=-30000., epsc0=-.002, fpcu=-6000., epsU=-.006,
            lam=.1, ft=2500., Ets=1500000.),
        bounds=dict(fpc=[-40000.,-20000.], epsc0=[-.003,-.0015],
            lam=[.05,.3])),
    'concrete_core': dict(kind='uniaxial', material='Concrete02',
        params=dict(fpc=-36000., epsc0=-.0025, fpcu=-10000., epsU=-.015,
            lam=.1, ft=2500., Ets=1500000.),
        bounds=dict(fpc=[-50000.,-30000.], epsc0=[-.004,-.002],
            epsU=[-.025,-.008])),
    'rebar_steel02': dict(kind='uniaxial', material='Steel02',
        params=dict(Fy=400000., E=200000000., b=.01, R0=18., cR1=.925,
            cR2=.15, a1=0., a2=1., a3=0., a4=1.),
        bounds=dict(Fy=[300000.,500000.], b=[.002,.03], R0=[10.,25.])),
    'pile_py': dict(kind='spring', material='PySimple1',
        params=dict(soil_type=2, pult=200., y50=.01, Cd=.1, c=0.),
        bounds=dict(pult=[100.,350.], y50=[.003,.025], Cd=[0.,.3])),
    'pile_tz': dict(kind='spring', material='TzSimple1',
        params=dict(soil_type=2, tult=100., z50=.005, c=0.),
        bounds=dict(tult=[50.,200.], z50=[.001,.02])),
    'pile_qz': dict(kind='spring', material='QzSimple1',
        params=dict(qz_type=2, qult=1000., z50=.01, suction=0., c=0.),
        bounds=dict(qult=[500.,2000.], z50=[.003,.03])),
    'normal_gap': dict(kind='spring', material='ElasticPPGap',
        params=dict(K=100000., Fy=-500., gap=-.001, eta=0.),
        bounds=dict(K=[50000.,200000.], Fy=[-900.,-200.], gap=[-.003,-.0001])),
    'interface_friction': dict(kind='contact', material='zeroLengthContactASDimplex',
        params=dict(Kn=1000000., Kt=200000., mu=.5),
        bounds=dict(Kt=[100000.,400000.], mu=[.2,.8])),
    'radiation_dashpot': dict(kind='viscous', material='Viscous',
        params=dict(C=400., alpha=1.), bounds=dict(C=[100.,800.])),
}

# Candidate shell materials. All except the explicitly identified geopolymer E
# are chosen demo seeds. Equal names never imply certified grades or recipes.
for name,E,nu,rho in [
    ('shell_concrete',30.e6,.2,2.4),
    ('shell_sintered_regolith',20.e6,.2,2.0),
    ('shell_geopolymer',8.13e6,.2,1.75),
    ('shell_sulfur_concrete',15.e6,.2,2.0),
]:
    CASES[name]=dict(kind='shell',material='ElasticIsotropic',
        params=dict(E=E,nu=nu,rho=rho,thickness=.02),
        bounds=dict(E=[.5*E,1.5*E]))
for name,E,nu,rho,fy in [('shell_steel',200.e6,.3,7.85,355000.),
                         ('shell_aluminum',69.e6,.33,2.7,250000.)]:
    CASES[name]=dict(kind='shell',material='J2Plasticity',
        params=dict(E=E,nu=nu,rho=rho,sig0=fy,sigInf=fy*1.2,
                    delta=15.,H=E*.005,thickness=.02),
        bounds=dict(sig0=[fy*.7,fy*1.3],H=[E*.001,E*.02]))
CASES['shell_cfrp_elastic']=dict(kind='shell',material='ElasticOrthotropic',
    params=dict(Ex=130.e6,Ey=9.e6,Ez=9.e6,nu_xy=.28,nu_yz=.35,nu_zx=.01938461538,
                Gxy=5.e6,Gyz=3.3e6,Gzx=5.e6,rho=1.6,thickness=.02),
    bounds=dict(Ex=[90.e6,170.e6]))
for planet,rho in [('moon',1.5),('mars',1.65)]:
    CASES[f'{planet}_regolith_dp']=dict(kind='soil3d',material='DruckerPrager',
        params=dict(rho=rho,G=20000.,K=45000.,cohesion=.5,phi=35.,psi=5.,
                    H=100.,theta=1.,pa=101.),
        bounds=dict(G=[5000.,50000.],cohesion=[.05,2.],phi=[28.,45.]))

def get_case(name):
    c = deepcopy(CASES[name])
    c.update(name=name, status='DEMONSTRATION_NOT_CALIBRATED',
             units='kN,m,s,tonne,kPa',
             test=dict(sigma_v=100., steps=40, amplitude=.005, cycles=3,
                       normal_force=100., substeps=1))
    if c['kind'] == 'spring':
        c['test']['amplitude'] = .05
    if name == 'elastic_soil':
        c['test']['amplitude'] = .0001
    if c['kind']=='shell' and c['material']!='J2Plasticity':
        c['test']['amplitude']=.0001
    if name.startswith(('moon_regolith','mars_regolith')):
        planet=name.split('_')[0]
        g={'moon':1.62,'mars':3.71}[planet]
        c['test']['sigma_v']=c['params']['rho']*g*1.0
        c['test']['substeps']=10
        c['planet']=dict(name=planet,gravity=g,depth=1.,surcharge=0.,
                         condition='dry_demo_not_site_calibrated')
    return c

def validate(c):
    p, t, m = c['params'], c['test'], c['material']
    for key in ('G','K','E','rho','G0','hpo','Su','p_ref','pa','pult','tult',
                'qult','y50','z50','Kn','Kt','C'):
        if key in p and p[key] <= 0:
            raise ValueError(f'{key} must be positive')
    if 'nu' in p and not (-1 < p['nu'] < .5):
        raise ValueError('-1 < nu < 0.5 required')
    if m == 'Concrete02' and not (p['fpc'] < p['fpcu'] <= 0 and
            p['epsU'] < p['epsc0'] < 0 and 0 < p['lam'] < 1 and
            p['ft'] >= 0 and p['Ets'] > 0):
        raise ValueError('invalid concrete compression/envelope parameters')
    if m == 'PressureDependMultiYield02' and not (0 < p['phi_pt'] <= p['phi'] < 60):
        raise ValueError('require 0 < phi_pt <= phi < 60')
    if m == 'PM4Sand' and not (0 < p['Dr'] < 1 and p['emax'] > p['emin'] > 0):
        raise ValueError('invalid relative density or void ratio')
    if m=='DruckerPrager' and not (0 < p['phi'] < 60 and 0 <= p['psi'] <= p['phi'] and p['cohesion']>0):
        raise ValueError('DP requires positive cohesion and 0 <= psi <= phi < 60')
    if c['kind']=='shell' and p['thickness']<=0:
        raise ValueError('shell thickness must be positive')
    if m=='J2Plasticity' and not (p['sigInf']>=p['sig0']>0 and p['H']>=0):
        raise ValueError('J2 requires sigInf >= sig0 > 0, H >= 0')
    if m == 'QzSimple1' and not 0 <= p['suction'] <= .1:
        raise ValueError('0 <= suction <= 0.1 required')
    if m == 'ElasticPPGap' and not (p['Fy'] < 0 and p['gap'] <= 0):
        raise ValueError('compression-only gap requires negative Fy and gap')
    if m == 'zeroLengthContactASDimplex' and p['mu'] < 0:
        raise ValueError('mu must be nonnegative')
    if c['kind'] == 'spring' and p.get('c',0.) != 0:
        raise ValueError('static spring driver cannot calibrate c; use dynamic SSI tests')
    if c['kind'] == 'viscous' and p['alpha'] <= 0:
        raise ValueError('alpha must be positive')
    if t['steps'] < 2 or t['cycles'] < 1 or t['substeps'] < 1:
        raise ValueError('steps >= 2, cycles >= 1, substeps >= 1 required')
    if t['sigma_v'] <= 0 or t['amplitude'] <= 0 or t['normal_force'] <= 0:
        raise ValueError('positive pressure, normal force and amplitude required')
