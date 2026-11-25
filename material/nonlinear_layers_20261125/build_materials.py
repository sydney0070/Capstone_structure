"""Literature-traceable, all-elastoplastic Mars-layer trial records.

DP simple-shear matching is distinct from triaxial-compression matching.
Rock-mass DP records approximate the published HB equations over a stated range;
GSI and analogue transfer are uncertainty scenarios, not InSight measurements.
"""
import copy
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=ROOT/'material/insight_layer_specimens_20261119/materials.json'


def shear_matching(c,phi):
    t=math.tan(math.radians(phi))
    sine=3*t/(2*math.sqrt(3)+t)
    angle=math.asin(sine)
    return c*(3-sine)/(2*math.sqrt(3)*math.cos(angle)),math.degrees(angle)


def hb_constants(gsi,mi=22.4,D=0.):
    # Canonical Hoek&Brown2019 Eqs4–6 use MINUS exp(-20/3).
    # Villeneuve2021 Eq8 prints a plus; the canonical equation is used here.
    return dict(mb=mi*math.exp((gsi-100)/(28-14*D)),s=math.exp((gsi-100)/(9-3*D)),
                a=.5+(math.exp(-gsi/15)-math.exp(-20/3))/6)


def hb_sigma1(pressure,constants,sigci):
    return pressure+sigci*(constants['mb']*pressure/sigci+constants['s'])**constants['a']


def rock_fit(gsi):
    porosity=.04
    sigci=33.32e6*porosity**(-.408)
    constants=hb_constants(gsi)
    x=np.linspace(0,500000,101);y=hb_sigma1(x,constants,sigci)
    # Minimize maximum relative error, rather than fitting only high-stress points.
    A=np.vstack([np.c_[np.ones(len(x)),x,-y],np.c_[-np.ones(len(x)),-x,-y]])
    fit=linprog([0.,0.,1.],A_ub=A,b_ub=np.r_[y,-y],bounds=[(0,None),(1,None),(0,None)],method='highs')
    if not fit.success:raise RuntimeError(fit.message)
    intercept,slope,error=fit.x
    c=intercept/(2*math.sqrt(slope));phi=math.degrees(math.asin((slope-1)/(slope+1)))
    return c,phi,dict(GSI=gsi,mi=22.4,D=0.,micro_porosity_analogue=.04,sigma_ci_Pa=sigci,
                      sigma_ci_source='Villeneuve2021 Table3 basalt empirical UCS=33.32*porosity^(-.408)MPa; SEE59.5MPa for full range',
                      mi_source='Villeneuve2021 Table2 Etna basalt porosity.04 mi22.4',GSI_source='55 recommended average volcanic rock mass;40/70 sensitivity assumptions, not Mars measurements',
                      constants=constants,fit_sigma3_range_Pa=[0.,500000.],fit_points=101,
                      fit_max_relative_error=float(error),fit_kind='minimax relative linear MC fit, then DP triaxial mapping',
                      independent_validation=False,
                      formula_note='Canonical a=.5+(exp(-GSI/15)-exp(-20/3))/6 from HoekBrown2019; plus sign in Villeneuve Eq8 not used.')


def build():
    old={r['id']:r for r in json.loads(OLD.read_text())['materials']}
    materials=[]
    def add(name,base,layer,c,phi,mode,target_source,extra=None):
        r=copy.deepcopy(old[base]);r.update(id=name,group='nonlinear_20261125',layer=layer,family='DruckerPrager',law='DruckerPrager',
                                          status='NONLINEAR_LITERATURE_CONSTRAINED_TRIAL_NOT_SITE_VALIDATED',mars_calibrated=False)
        p=r['parameters'];target_c,target_phi=c,phi
        if mode=='simple_shear':c,phi=shear_matching(c,phi)
        p.update(c_Pa=float(c),phi_deg=float(phi),psi_deg=0.,hardening=0.)
        r['calibration']=dict(target_mode=mode,target_MC_c_Pa=target_c,target_MC_phi_deg=target_phi,target_source=target_source,independent_validation=False)
        if extra:r['calibration'].update(extra)
        for key in ['c_Pa','phi_deg']:
            r['evidence'][key]=dict(source=target_source,locator='see calibration record',kind='derived_calibration' if mode=='simple_shear' or extra else 'literature_range_scenario',
                                   note='Numerical DP input; for simple-shear fitting it is not the physical MC friction/cohesion value.')
        r['evidence']['psi_deg']=dict(source='modelling_assumption',locator='nonassociated flow',kind='assumed',note='Zero dilation; requires volumetric-strain calibration')
        r['evidence']['hardening']=dict(source='modelling_assumption',locator='perfect plasticity',kind='assumed',note='No measured hardening/softening curve')
        r['notes']=['Actual DP c/phi and target physical MC c/phi are stored separately.',
                    'Elastic moduli remain source-based initial proxies; no unmeasured strain-softening law is claimed.',
                    'Elastoplastic law may remain on its elastic branch at service loads.']
        materials.append(r)
        return r
    for mode in ['triaxial','simple_shear']:
        add('crust_'+('tc' if mode=='triaxial' else 'ss'),'duricrust_dp','duricrust',5800.,30.,mode,'Spohn2022 slope inference at assumed phi30/rho1200; transferred to rho1025')
    for code,A,k,phi,rmin in [('mgs1',4.96,2.97,35.5,1.32),('mgs1s',5.28,2.93,34.,1.31)]:
        target=A*math.exp(k*1.30)
        for mode in ['triaxial','simple_shear']:
            add(code+'_'+('tc' if mode=='triaxial' else 'ss'),'sand_'+code+'_dp_1300','sand',target,phi,mode,'Dotson2024 Eq2/Table4; Fig6 nearest-density friction',
                dict(A_Pa=A,k=k,target_density_g_cm3=1.30,measured_density_min_g_cm3=rmin,density_extrapolated=True,
                     uncertainty_note='Slight density extrapolation; phi held at nearby density, fit covariance unavailable'))
        add(code+'_ss_nearest','sand_'+code+'_dp_1300','sand',A*math.exp(k*rmin),phi,'simple_shear','Dotson2024 Table4 fit at nearest tested density, transferred to1300',
            dict(A_Pa=A,k=k,reference_density_g_cm3=rmin,target_density_g_cm3=1.30,density_transferred=True))
    for phi in [20.,30.,35.,40.]:
        add(f'blocky_phi{int(phi)}','blocky_ejecta_dp_trial','blocky_ejecta',0.,phi,'triaxial','Hughes2022 Section6.2 broken/blasted cohesionless rock phi20–40deg; transferred from Gale/terrestrial range to InSight scenario')
    for c in [1000.,5000.]:
        add(f'blocky_c{int(c)}_phi30','blocky_ejecta_dp_trial','blocky_ejecta',c,30.,'triaxial','Hughes2022 granular regolith c1–5kPa; bonded-mixture sensitivity, not measured InSight blocky ejecta')
    for phi in [30.,35.]:
        add(f'blocky_ss_phi{int(phi)}','blocky_ejecta_dp_trial','blocky_ejecta',0.,phi,'simple_shear','Hughes2022 phi20–40deg range; simple-shear-matched zero-cohesion scenario')
    for gsi in [40,55,70]:
        c,phi,extra=rock_fit(gsi)
        add(f'basalt_hb_gsi{gsi}','fractured_basalt_elastic','fractured_basalt',c,phi,'triaxial','VilleneuveHeap2021 + HoekBrown2019 rock-mass analogue',extra)
    r=add('basalt_lab_reference','fractured_basalt_elastic','laboratory_reference_only',46.5e6,46.7,'triaxial','VilleneuveHeap2021 Table1 Etna basalt lava porosity.04; c46.5±3.95MPa phi46.7±.95deg',dict(c_uncertainty_Pa=3.95e6,phi_uncertainty_deg=.95))
    r['notes'].append('Strength-reference fixture only; E13GPa/rho2760 inherited to isolate yield-envelope implementation, not measured Etna elastic data.')
    payload=dict(date='2026-11-25',units='N m kg s Pa',source_specimens=str(OLD),source_sha256=hashlib.sha256(OLD.read_bytes()).hexdigest(),materials=materials,
                 scope=f'All{len(materials)} records use nD DruckerPrager elastoplasticity. Calibration/analogue/sensitivity are labelled; no independent InSight validation claimed.')
    (HERE/'materials.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    with (HERE/'parameters.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['specimen','layer','E_MPa','density_kg_m3','input_DP_c_Pa','input_DP_phi_deg','target_MC_c_Pa','target_MC_phi_deg','mode','source','HB_fit_error_pct'])
        for r in materials:
            p=r['parameters'];c=r['calibration'];w.writerow([r['id'],r['layer'],p['E_Pa']/1e6,p['density_kg_m3'],p['c_Pa'],p['phi_deg'],c['target_MC_c_Pa'],c['target_MC_phi_deg'],c['target_mode'],c['target_source'],100*c['fit_max_relative_error'] if 'fit_max_relative_error' in c else ''])
    print('Created',len(materials),'nonlinear records')


if __name__=='__main__':build()
