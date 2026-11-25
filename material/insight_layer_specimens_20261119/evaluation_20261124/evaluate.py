"""Audit literature inputs, independent formulas, and recorded specimen behaviour.

No base input/result is changed. Null reference/error means not identifiable,
not zero error. All project tolerances are separated from published bounds.
"""
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
PACK=HERE.parent
ROOT=PACK.parents[1]
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/capstone_specimen_mpl')


def save_csv(name,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (HERE/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def pct(value,reference):
    return 100*(value/reference-1) if reference else None


def write_json(path,data):
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def main():
    records=json.loads((PACK/'materials.json').read_text())['materials']
    original=json.loads((PACK/'results/summary.json').read_text())
    results={r['id']:r for r in original['cases']}
    sources=json.loads((PACK/'source_manifest.json').read_text())
    for key,file,url in [('QzSimple1_command','D15_QzSimple1_20261119.html','QzSimple1'),('PySimple1_command','D13_PySimple1_20261119.html','PySimple1'),('TzSimple1_command','D14_TzSimple1_20261119.html','TzSimple1')]:
        p=ROOT/'docs/references/05_official_model_docs'/file
        sources.append(dict(id=key,file=str(p),url=f'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/{url}.html',exists=p.is_file(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    extra=json.loads((HERE/'additional_tests/summary.json').read_text())
    criteria=[
        dict(id='N_E',metric='initial E',formula='abs(Efit/Ein-1)',tolerance_relative=1e-5,tolerance_pct=.001,scope='numerical_formula',basis='Existing project tolerance; not a paper validation standard'),
        dict(id='N_Y',metric='DP triaxial yield',formula='abs(qFE/qanalytic-1)',tolerance_relative=1e-4,tolerance_pct=.01,scope='numerical_formula',basis='Existing project tolerance; MC-matched triaxial compression only'),
        dict(id='N_SHEAR',metric='DP fully developed simple-shear yield',formula='tau=2sqrt(3)*cos(phi)/(3-sin(phi)) * (c+p*tan(phi))',tolerance_relative=1e-4,tolerance_pct=.01,scope='added_numerical_formula',basis='Derived from the implemented DP yield surface at equal normal stresses and zero dilation; not a paper strength fit'),
        dict(id='N_M',metric='cube mass',formula='abs((spring_k/eigenvalue)/(rho*1m3)-1)',tolerance_relative=1e-7,tolerance_pct=.00001,scope='numerical_formula',basis='Existing project tolerance; density implementation, not a density experiment'),
        dict(id='N_L',metric='lateral pressure',formula='max(abs(sigma_x+p),abs(sigma_y+p))/p',tolerance_relative=1e-5,tolerance_pct=.001,scope='numerical_boundary',basis='Existing tolerance'),
        dict(id='N_F',metric='force balance',formula='max(abs(F+R))/max(max(abs(F)),1N); triaxial scale=max(qpeak*A,p*A,1N)',tolerance_relative=1e-5,tolerance_pct=.001,scope='numerical_equilibrium',basis='Existing tolerance'),
        dict(id='N_V',metric='constant vertical pressure',formula='max(abs(sigma_z+p))/p',tolerance_relative=1e-5,tolerance_pct=.001,scope='numerical_boundary',basis='Existing tolerance'),
        dict(id='N_U',metric='DP elastic unloading to q0',formula='abs(qfinal)/qyield',tolerance_relative=1e-4,tolerance_pct=.01,scope='numerical_formula',basis='Existing tolerance; residual strain imposed at 2*qyield/E'),
        dict(id='N_C',metric='elastic curve',formula='abs(q-E*eps) <= .01Pa + 1e-6*abs(E*eps)',tolerance_relative=1e-6,tolerance_pct=.0001,scope='numerical_formula',basis='Existing mixed absolute/relative tolerance'),
        dict(id='N_0',metric='elastic return to origin',formula='abs(qfinal) < .01Pa',tolerance_absolute=.01,unit='Pa',scope='numerical_formula',basis='Existing absolute tolerance; relative error undefined at zero'),
        dict(id='N_CAP',metric='spring capacity bound',formula='max(abs(F))/Fult <= 1.001',tolerance_relative=.001,tolerance_pct=.1,scope='constitutive_bound',basis='Existing upper bound; reaching Fult is not required at finite displacement'),
        dict(id='N_GAP',metric='ElasticPPGap capacity',formula='abs(Fmin/Fyield-1)',tolerance_relative=1e-6,tolerance_pct=.0001,scope='synthetic_mechanics',basis='Existing tolerance; not a measured contact force'),
        dict(id='N_HYS',metric='Hysteretic initial stiffness',formula='abs(KFE/Kinput-1)',tolerance_relative=1e-5,tolerance_pct=.001,scope='synthetic_mechanics',basis='Existing tolerance'),
        dict(id='L_RANGE',metric='reported interval membership',formula='reported_low <= trial_value <= reported_high',scope='literature_plausibility',basis='Published range; no assumed confidence level or layer-specific validation'),
        dict(id='L_ENV',metric='approximate MC envelope difference',formula='100*(tau_FE/tau_reference-1)',tolerance_pct=5.,scope='provisional_literature_target',basis='NEW project screening target, not a paper/ASTM tolerance or statistical acceptance criterion. Density/fixture mismatch remains.'),
        dict(id='L_CURVE',metric='experimental cyclic/force-displacement curve RMSE',formula='RMSE=sqrt(mean((prediction-observation)^2))',scope='not_evaluable',basis='Matched raw experimental curves unavailable. Not calculated as zero.'),
        dict(id='L_50',metric='nominal spring half-capacity point',formula='100*(abs(F(z50))/(.5*Fult)-1)',scope='nominal_diagnostic',basis='Nominal shape diagnostic, not a universal pass tolerance. Simple1 docs do not establish exact Simple2 shape equivalence.'),
    ]
    numerical=[];behaviour=[]
    def metric(name,case,quantity,ref,val,tol=None,unit='',kind='numerical_formula',absolute_tol=None):
        delta=val-ref if ref is not None else None
        rel=abs(delta/ref) if ref else None
        passed=(rel<=tol if tol is not None and rel is not None else abs(delta)<=absolute_tol if absolute_tol is not None and delta is not None else None)
        utilization=None
        if kind=='finite_amplitude_utilization_not_error':
            utilization=100*val/ref;rel=None;passed=val<=ref*1.001;tol=.001
        numerical.append(dict(specimen=name,case=case,quantity=quantity,unit=unit,reference=ref,computed=val,absolute_difference=delta,relative_error=rel,error_pct=100*rel if rel is not None else None,tolerance_relative=tol,tolerance_absolute=absolute_tol,passed=passed,comparison_kind=kind,utilization_pct=utilization))
    for r in records:
        name=r['id'];p=r['parameters'];res=results[name]
        if 'mass_test' in res:
            m=res['mass_test'];metric(name,'1m3_rigid_translation','mass',m['expected_kg'],m['recovered_kg'],1e-7,'kg')
        for t in res['tests']:
            pressure=t.get('pressure_Pa');case=f'p{int(pressure)}' if pressure else 'spring_history'
            row=dict(specimen=name,family=r['family'],case=case,pressure_Pa=pressure,numerical_pass=t['passed'])
            if t['test']=='triaxial_load_unload':
                a=np.loadtxt(PACK/'results'/name/f'triaxial_p{int(pressure)}.csv',delimiter=',',skiprows=1)
                Efit=float(np.polyfit(a[:20,1],a[:20,2],1)[0]);qmax=float(np.max(a[:,2]));scale=max(qmax,pressure,1.)
                metric(name,case,'E_initial',p['E_Pa'],Efit,1e-5,'Pa')
                metric(name,case,'lateral_pressure_max_error',1.,1.+float(t['confinement_max_error_Pa'])/pressure,1e-5,'normalized')
                metric(name,case,'force_balance_max_error',1.,1.+float(t['force_balance_max_error_N'])/scale,1e-5,'normalized')
                row.update(initial_E_MPa=Efit/1e6,peak_q_kPa=qmax/1000,final_q_Pa=float(a[-1,2]),final_axial_strain=float(a[-1,1]))
                if r['family']=='DruckerPrager':
                    sine=math.sin(math.radians(p['phi_deg']))
                    qref=(2*pressure*sine+2*p['c_Pa']*math.cos(math.radians(p['phi_deg'])))/(1-sine)
                    metric(name,case,'q_yield',qref,float(a[150,2]),1e-4,'Pa')
                    Eu=float(np.polyfit(a[151:,1],a[151:,2],1)[0]);metric(name,case,'E_unloading',p['E_Pa'],Eu,1e-5,'Pa',kind='added_formula_diagnostic')
                    metric(name,case,'unloaded_q_normalized',1.,1.+abs(float(a[-1,2]))/qref,1e-4,'normalized')
                    metric(name,case,'residual_axial_strain',2*qref/p['E_Pa'],float(a[-1,1]),1e-5,'strain',kind='imposed_protocol_consistency_not_experiment')
                    row.update(analytic_yield_q_kPa=qref/1000,unloading_E_MPa=Eu/1e6)
                else:
                    abs_curve=float(np.max(abs(a[:,2]-p['E_Pa']*a[:,1])))
                    metric(name,case,'elastic_curve_max_absolute_error',0.,abs_curve,unit='Pa',kind='numerical_formula',absolute_tol=.01+1e-6*float(np.max(abs(p['E_Pa']*a[:,1]))))
                    metric(name,case,'return_to_origin_q',0.,float(a[-1,2]),unit='Pa',absolute_tol=.01)
            elif t['test'].startswith('drained_'):
                a=np.loadtxt(PACK/'results'/name/f'cyclic_shear_p{int(pressure)}.csv',delimiter=',',skiprows=1)
                scale=max(float(np.max(abs(a[:,3]))),1.)
                eq=float(np.max(abs(a[:,3]+a[:,4]))/scale);vert=float(np.max(abs(a[:,2]+pressure))/pressure)
                metric(name,case,'force_balance_max_error',1.,1.+eq,1e-5,'normalized')
                metric(name,case,'vertical_stress_max_error',1.,1.+vert,1e-5,'normalized')
                initG=float(np.polyfit(a[:5,0],a[:5,1],1)[0]);taupos=float(a[419,1])
                row.update(initial_G_MPa=initG/1e6,peak_abs_tau_kPa=float(np.max(abs(a[:,1])))/1000,
                           tau_at_zero_final_kPa=float(a[-1,1])/1000,secant_G_at_first_positive_0p5pct_MPa=taupos/.005/1e6,
                           final_vertical_displacement_mm=float(a[-1,5])*1000,initial_mean_pressure_kPa=t['initial_mean_pressure_Pa']/1000,
                           bisections=t['step_bisections'],initial_G_literature_error=None,cyclic_curve_RMSE=None)
            else:
                a=np.loadtxt(PACK/'results'/name/'cyclic_spring.csv',delimiter=',',skiprows=1)
                row.update(positive_force_kN=float(np.max(a[:,1]))/1000,negative_force_kN=float(np.min(a[:,1]))/1000,final_force_kN=float(a[-1,1])/1000)
                if 'ultimate_force_N' in p:
                    F=p['ultimate_force_N'];observed=float(np.max(abs(a[:,1])))
                    row.update(capacity_kN=F/1000,mobilized_capacity_pct=100*observed/F,capacity_shortfall_pct=100*(1-observed/F),nominal_half_capacity_difference_pct=pct(abs(float(a[199,1])),F/2),force_at_virgin_negative_z50_kN=float(a[199,1])/1000)
                    metric(name,case,'capacity_upper_bound',F,observed,unit='N',kind='finite_amplitude_utilization_not_error')
                if r['family']=='ElasticPPGap':metric(name,case,'compression_yield_force',p['yield_force_N'],float(np.min(a[:,1])),1e-6,'N',kind='synthetic_mechanics')
                if r['family']=='Hysteretic':
                    K=p['positive_backbone'][0][0]/p['positive_backbone'][0][1]
                    metric(name,case,'initial_spring_stiffness',K,float(a[0,1]/a[0,0]),1e-5,'N/m',kind='synthetic_mechanics')
            behaviour.append(row)

    # Numeric source references are explicit; assumptions have no fabricated target.
    refs={}
    def put(name,key,reference,source,locator,kind='source_input',lo=None,hi=None,note='',raw=None):
        refs[(name,key)]=dict(reference_value=reference,source=source,locator=locator,comparison_kind=kind,reported_low=lo,reported_high=hi,note=note,raw_reference_value=raw)
    for r in records:
        name=r['id'];p=r['parameters'];family=r['family']
        if r['group']=='first_comparison':
            if name.startswith('surface_'):
                put(name,'E_Pa',7.5e6,'Kenda2020','Table1, PDF4','effective_reference_transfer');put(name,'nu',.22,'Kenda2020','Table1, PDF4','effective_reference_transfer');put(name,'density_kg_m3',1200.,'Spohn2022','Section6.7','inferred_reference')
            elif name.startswith('duricrust_'):
                put(name,'E_Pa',11.48e6,'Brinkman2022','Table4, PDF13','effective_shallow_proxy',9.25e6,17.39e6,note='Not isolated duricrust')
                put(name,'nu',.28,'Brinkman2022','Table4, PDF13','effective_shallow_proxy',.229,.40)
                put(name,'density_kg_m3',None,'Spohn2022','Section6.7','inferred_range',950.,1100.,'Midpoint chosen; no point measurement')
                if family=='DruckerPrager':
                    put(name,'c_Pa',5800.,'Spohn2022','Section6.4','conditional_inference',2000.,15000.,'Conditional on phi30 and rho1200; applied rho1025')
                    put(name,'phi_deg',30.,'Spohn2022','Section6.4','source_assumption_not_measurement')
            elif name.startswith('sand_'):
                # Recompute from published predicted waves, at the chosen target density.
                rho=1300.;G=rho*94.6**2;K=rho*(157.9**2-4*94.6**2/3);E=9*K*G/(3*K+G);nu=(3*K-2*G)/(2*(3*K+G))
                for key,value in [('E_Pa',E),('nu',nu),('G_Pa',G),('K_Pa',K)]:put(name,key,value,'Murdoch2018','Table2, PDF15–16','derived_at_target_density',note='Predicted 1m state, not measured static modulus')
                put(name,'density_kg_m3',None,'Spohn2022','Section6.7','inferred_range',1300.,1500.)
                if family=='DruckerPrager':
                    put(name,'phi_deg',35.5 if 'mgs1_dp' in name else 34.,'Dotson2024','Figure6','approximate_graph_transfer')
                    put(name,'c_Pa',280.,'Dotson2024','Figure5, earlier graph reading','approximate_graph_transfer',note='Approximate graph reading; Table4 fit checked separately below')
            elif name.startswith('local_gravel_'):
                put(name,'density_kg_m3',1600.,'Spohn2022','Section6.7','mixture_estimate')
            elif name.startswith('blocky_'):
                put(name,'E_Pa',520e6,'Kenda2020','Table1, PDF4','expected_effective_reference');put(name,'nu',.24,'Kenda2020','Table1, PDF4','expected_effective_reference')
            elif name.startswith('fractured_'):
                put(name,'E_Pa',13e9,'Kenda2020','Table1, PDF4','expected_effective_reference');put(name,'nu',.28,'Kenda2020','Table1, PDF4','expected_effective_reference');put(name,'density_kg_m3',2760.,'Morgan2018','PDF39','terrestrial_analogue')
            elif name.startswith('sediment_'):
                put(name,'Vp_m_s',750.,'Hobiger2021','Table2, PDF7','prior_not_inversion');put(name,'Vs_m_s',417.,'Hobiger2021','Table2, PDF7','prior_not_inversion')
        if r['group']=='followup_reference':
            dr=int(name.split('_dr')[1].split('_')[0]);i=[30,40,50,60,75].index(dr)
            raw={'density_kg_m3':[1700.,1800.,1900.,2000.,2100.][i], 'G_Pa':[60e6,90e6,100e6,110e6,130e6][i], 'K_Pa':[160e6,220e6,233e6,240e6,260e6][i],
                 'phi_deg':[31.,32.,33.5,35.,36.5][i], 'PTAngle_deg':[31.,26.,25.5,26.,26.][i], 'contrac1':[.087,.067,.045,.028,.013][i], 'contrac3':[.18,.23,.15,.05,0.][i], 'dilat1':[0.,.06,.06,.1,.3][i], 'dilat3':[0.,.27,.15,.05,0.][i], 'e':[.85,.77,.70,.65,.55][i],
                 'p_ref_Pa':101000.,'peak_shear_strain':.1,'pressure_exponent':.5,'n_surfaces':20,'contrac2':5.,'dilat2':3.,'liquefac1':1.,'liquefac2':0.,'cs1':.9,'cs2':.02,'cs3':.7,'pa_Pa':101000.,'numerical_c_Pa':300.}
            for key,value in raw.items():
                converted=value*(101/80)**.5 if key in ('G_Pa','K_Pa') else value
                put(name,key,converted,'PDMY02','Suggested parameters table / defaults','converted_earth_reference' if key in ('G_Pa','K_Pa') else 'earth_reference',note='No Mars cyclic calibration',raw=value)
        if r['family'] in ['QzSimple1','QzSimple2','PySimple2','TzSimple2']:
            for key,value in p.items():
                if key=='damping':put(name,key,.05,'Shallow2008','PDF17','static_override',note='Original Tz example has extra optional argument; only documented one damping argument used in fixture')
                else:put(name,key,value,'Shallow2008','PDF17','same_example_input' if family!='QzSimple1' else 'transfer_between_backbones')
    inputs=[]
    for r in records:
        for key,value in r['parameters'].items():
            ev=r['evidence'][key]
            ref=refs.get((r['id'],key),dict(reference_value=None,source=ev['source'],locator=ev['locator'],comparison_kind=ev['kind'],reported_low=None,reported_high=None,note=ev['note'],raw_reference_value=None))
            rv=ref['reference_value'];numeric=isinstance(value,(int,float)) and not isinstance(value,bool)
            within=(ref['reported_low']<=value<=ref['reported_high']) if numeric and ref['reported_low'] is not None and ref['reported_high'] is not None else None
            inputs.append(dict(specimen=r['id'],parameter=key,input_value=json.dumps(value,ensure_ascii=False),**ref,signed_difference=value-rv if numeric and rv is not None else None,signed_difference_pct=pct(value,rv) if numeric and rv is not None else None,within_reported_interval=within,experimental_validation=False))
    # Compare inferred wave speeds only with the appropriate source state; retain
    # the same-source equality as a consistency check, never an independent test.
    waves=[]
    def wave_row(r,quantity,value,ref,source,lo=None,hi=None,context=''):
        waves.append(dict(specimen=r['id'],quantity=quantity,unit='m/s',predicted=value,reference=ref,signed_difference_pct=pct(value,ref),reported_low=lo,reported_high=hi,within_reported_interval=lo<=value<=hi if lo is not None else None,source=source,context=context,independent_validation=False))
    for r in records:
        p=r['parameters']
        if r['group']!='first_comparison':continue
        vs=math.sqrt(p['G_Pa']/p['density_kg_m3']);vp=math.sqrt((p['K_Pa']+4*p['G_Pa']/3)/p['density_kg_m3'])
        if r['id'].startswith(('surface_','duricrust_')):
            wave_row(r,'Vs',vs,63.,'Brinkman2022',56.,74.,'Effective shallow-volume plausibility; not isolated-layer measurements')
            wave_row(r,'Vp',vp,119.,'Brinkman2022',98.,164.,'Effective shallow-volume plausibility; not isolated-layer measurements')
        if r['id'].startswith('sand_'):
            wave_row(r,'Vs',vs,94.6,'Murdoch2018',context='Predicted1m source used for input; rounding consistency only')
            wave_row(r,'Vp',vp,157.9,'Murdoch2018',context='Predicted1m source used for input; rounding consistency only')
        if r['id']=='fractured_basalt_elastic':
            wave_row(r,'Vs',vs,1700.,'Morgan2018/Hobiger2021',context='Different model reference: compare cross-source scenario, not a fit requirement')
            wave_row(r,'Vp',vp,3000.,'Morgan2018/Hobiger2021',context='Different model reference: compare cross-source scenario, not a fit requirement')
        if r['id']=='sediment_elastic_prior':
            wave_row(r,'Vs',vs,417.,'Hobiger2021',context='Same prior used to calculate inputs: algebraic identity')
            wave_row(r,'Vp',vp,750.,'Hobiger2021',context='Same prior used to calculate inputs: algebraic identity')
    # Table4 coefficients read visually from the supplied v1 PDF. rho in g/cm3;
    # A in Pa. Below measured density bounds: mark extrapolation explicitly.
    coefficients={'mgs1':dict(A_Pa=4.96,A_reported_uncertainty=.34,k=2.97,k_reported_uncertainty=.22,rho_min=1.32,rho_max=1.86,phi_deg=35.5),
                  'mgs1s':dict(A_Pa=5.28,A_reported_uncertainty=1.14,k=2.93,k_reported_uncertainty=.22,rho_min=1.31,rho_max=1.83,phi_deg=34.)}
    write_json(HERE/'dotson_table4_extracted.json',dict(source='Dotson2024 archived v1 PDF12 Table4; Eq2 PDF10; Fig6 PDF13',units='cohesion Pa, density g/cm3, phi degrees',materials=coefficients,note='Fit covariance/raw point residuals not supplied here; no confidence interval propagated. phi graph/text value at nearest measured density, not a fit at1300.'))
    strength=[];strength_stats=[]
    for t in extra['cases']:
        if not t.get('completed'):continue
        key='mgs1s' if 'mgs1s' in t['id'] else 'mgs1';c=coefficients[key];pressure=t['pressure_Pa'];phi=math.radians(c['phi_deg'])
        cref=c['A_Pa']*math.exp(c['k']*1.3);cnear=c['A_Pa']*math.exp(c['k']*c['rho_min'])
        actual=t['final_tau_Pa'];old=280.+pressure*math.tan(phi);target=cref+pressure*math.tan(phi);near=cnear+pressure*math.tan(phi)
        dp_theory=None;dp_theory_difference=None
        if '_dp_' in t['id']:
            factor=2*math.sqrt(3)*math.cos(phi)/(3-math.sin(phi))
            dp_theory=factor*old;dp_theory_difference=pct(actual,dp_theory)
            metric(t['id'],f'additional_p{int(pressure)}','DP_steady_simple_shear',dp_theory,actual,1e-4,'Pa',kind='added_numerical_formula')
        strength.append(dict(specimen=t['id'],normal_stress_Pa=pressure,gamma=t['final_gamma'],computed_tau_Pa=actual,approx_graph_MC_tau_Pa=old,graph_envelope_signed_difference_pct=pct(actual,old),
                             dp_steady_shear_theory_Pa=dp_theory,dp_formula_difference_pct=dp_theory_difference,
                             table4_nearest_density_MC_tau_Pa=near,table4_nearest_signed_difference_pct=pct(actual,near),
                             table4_1300_extrapolated_MC_tau_Pa=target,table4_1300_signed_difference_pct=pct(actual,target),
                             c_input_Pa=280. if '_dp_' in t['id'] else None,c_table4_1300_Pa=cref,c_table4_nearest_Pa=cnear,
                             extrapolated_density=True,reference_density_min_kg_m3=1000*c['rho_min'],target_density_kg_m3=1300,
                             preliminary_5pct_target_met=abs(pct(actual,old))<=5.,numerical_pass=t['numerical_pass'],tail_variation_pct=100*t['tail_range_relative'],
                             interpretation='MC-fit envelope comparison only; DP near plateau / PDMY endpoint at0.2% may still harden; not raw experimental curve RMSE'))
    for name in sorted(set(x['specimen'] for x in strength)):
        rows=[x for x in strength if x['specimen']==name];a=np.array([x['computed_tau_Pa'] for x in rows]);b=np.array([x['approx_graph_MC_tau_Pa'] for x in rows]);b2=np.array([x['table4_1300_extrapolated_MC_tau_Pa'] for x in rows])
        strength_stats.append(dict(specimen=name,n=len(rows),graph_reference_MAPE_pct=float(np.mean(abs((a/b-1)*100))),graph_reference_max_error_pct=float(np.max(abs((a/b-1)*100))),graph_reference_RMSE_Pa=float(np.sqrt(np.mean((a-b)**2))),
                                   graph_reference_NRMSE_pct=float(np.sqrt(np.mean((a-b)**2))/np.mean(b)*100),table4_extrapolated_MAPE_pct=float(np.mean(abs((a/b2-1)*100))),
                                   note='RMSE across four normal-stress envelope points; not a stress-strain time-series RMSE or independent validation'))
    assessments=[]
    for r in records:
        name=r['id'];family=r['family']
        if family=='PressureDependMultiYield02':
            decision='추가 보정 필요';reason='반복 곡선의 독립 문헌 오차 산정 불가. 지구 기준 계수/화성 이전 가정을 분리. 같은 평균압 기준도 추가 확인 필요.'
        elif family=='DruckerPrager':
            decision='정적 비교용 초안 · 강도 보정 필요';reason='삼축압축 구성식 일치. 직접전단/단순전단 강도는 같은 c·phi라도 달라질 수 있음.'
        elif family=='ElasticIsotropic':
            decision='탄성 기준 시편';reason='입력·구성식 구현 확인. 파속/원자료의 공간 범위·상태 차이와 대용값 가정을 유지.'
        else:
            decision='기능 비교용 · 기초 적용 보정 필요';reason='문헌 예제/합성 곡선. 화성 링 기초의 분담면적·하중–침하 곡선은 미보정.'
        if name=='surface_elastic' or name=='surface_dp_trial':reason+=' 선택 물성의 Vs50.61m/s와 Vp84.47m/s는 Brinkman 천부 유효 보고범위보다 낮음. 얇은 먼지층 직접값으로 확정 불가.'
        if name=='duricrust_dp':reason+=' c5.8kPa는 phi30·rho1200 조건 추정값을 rho1025에 이전. 경화층 자체 E 실측 아님.'
        if name=='fractured_basalt_elastic':reason+=' E13GPa·rho2760의 조합은 Vs1356m/s이며 다른 참조의1700m/s와 약20.21% 차이.'
        if name=='sediment_elastic_prior':reason+=' 파속 원값과0% 차이는 같은 사전모델로 역산한 항등관계.'
        if name in [x['specimen'] for x in strength_stats]:
            st=next(x for x in strength_stats if x['specimen']==name);reason+=f' 추가 전단의 기존 그림 포락선 대비 평균 절대차 {st["graph_reference_MAPE_pct"]:.2f}%; 잠정5% 목표 초과.'
        assessments.append(dict(specimen=name,family=family,base_numerical_pass=results[name]['passed'],independent_Mars_validation=False,experimental_curve_RMSE=None,decision=decision,reason=reason))

    save_csv('criteria.csv',criteria);save_csv('numerical_errors.csv',numerical);save_csv('literature_input_comparison.csv',inputs);save_csv('wave_speed_comparison.csv',waves)
    save_csv('strength_envelope_comparison.csv',strength);save_csv('strength_statistics.csv',strength_stats);save_csv('behaviour_metrics.csv',behaviour);save_csv('specimen_assessments.csv',assessments)
    worst={}
    for q in ['E_initial','q_yield','mass','force_balance_max_error','vertical_stress_max_error','E_unloading','DP_steady_simple_shear']:
        rr=[x for x in numerical if x['quantity']==q and x['error_pct'] is not None];w=max(rr,key=lambda x:x['error_pct']);worst[q]=dict(specimen=w['specimen'],case=w['case'],max_error_pct=w['error_pct'])
    outcome=dict(evaluation_date='2026-11-24',base_specimens=28,base_loading_histories=50,base_mass_checks=22,base_numerical_pass_count=sum(r['passed'] for r in results.values()),
                 additional_shear_tests=len(extra['cases']),additional_numerical_pass=sum(r.get('numerical_pass',False) for r in extra['cases']),independently_validated_Mars_specimens=0,
                 interpretation='Base passes verify formulas/equilibrium only. Input agreement and MC envelope comparison are not independent experimental validation.',worst_numerical_errors=worst,
                 strength_statistics=strength_stats,assessments=assessments)
    (HERE/'결과요약.txt').write_text('시편 평가 결과 — 2026-11-24\n\n기존 28개 수치 검사, 추가 문헌 압력 전단16개 모두 수치 검사 통과.\n문헌 강도 일치는 별도이며, 독립 화성 실험곡선 검증 완료 시편은0개.\n\n평가기준(프로젝트 수치 허용오차)\n초기 E:0.001%, DP 삼축압축 항복:0.01%, 질량:0.00001%, 힘 평형/압력:0.001%.\n문헌 전단 포락선 잠정 목표 ±5%는 이번 평가용 목표이며 논문/ASTM 허용기준이 아님.\n\n기존 그림에서 읽은 MC 포락선 대비 평균 절대차(MAPE)\nMGS-1 DP:16.5701%, MGS-1S DP:17.6607%\nMGS-1 PDMY02:18.2608%, MGS-1S PDMY02:17.3384%\n97/250/500/680Pa 4개 정상응력, PDMY02는 γ0.2% 시점(극한강도 아님).\n\nDP 차이의 주원인: 삼축압축 기준 변환의 단순전단 강도 계수1.165702/1.176607.\nDP 단순전단 해석 자체는 그 이론식과 최대 약0.000105% 차이.\n\nDotson 표4/식2를 ρ1300에 적용하면 c235.664/238.156Pa.\n기존 c280Pa와도 차이가 있으나 원 시험 최소밀도1320/1310보다 낮은 작은 외삽이므로 독립 검증값으로 확정하지 않음.\n\n나머지 층/반복/스프링에서 대응하는 독립 실험곡선이 없으면 오차는 산정 불가.\n0%로 채우지 않았음. 각 후보의 수치와 근거는 report.html 및8개 CSV에 있음.\n')
    write_json(HERE/'evaluation.json',outcome)
    plot(strength,waves,behaviour)
    report(outcome,criteria,inputs,numerical,waves,strength,strength_stats,behaviour,assessments,sources)
    inputs_manifest=[]
    for p in [PACK/'materials.json',PACK/'specimens.py',PACK/'results/summary.json',HERE/'additional_tests/summary.json']+[PACK/'results'/r['id']/'summary.json' for r in records]:
        inputs_manifest.append(dict(file=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    write_json(HERE/'input_manifest.json',inputs_manifest)
    write_json(HERE/'source_manifest.json',sources)
    manifest=[]
    for p in sorted(HERE.rglob('*')):
        if p.is_file() and p.name!='file_manifest.json' and '__pycache__' not in p.parts:manifest.append(dict(file=str(p.relative_to(HERE)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
    write_json(HERE/'file_manifest.json',manifest)
    print(json.dumps(dict(worst_errors=worst,strength_statistics=strength_stats,base_pass=outcome['base_numerical_pass_count'],extra_pass=outcome['additional_numerical_pass']),indent=2))


def plot(strength,waves,behaviour):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    for ax,key,title in zip(axes,['mgs1','mgs1s'],['MGS-1','MGS-1S']):
        dp=sorted([x for x in strength if x['specimen']==f'sand_{key}_dp_1300'],key=lambda x:x['normal_stress_Pa'])
        pd=sorted([x for x in strength if x['specimen']==f'sand_{key}_pdmy02_1300'],key=lambda x:x['normal_stress_Pa'])
        x=np.array([r['normal_stress_Pa'] for r in dp])/1000
        for name,rows,k,style in [('DP',dp,'computed_tau_Pa','o-'),('PDMY02 at 0.2% strain',pd,'computed_tau_Pa','s-'),('Earlier graph MC envelope',dp,'approx_graph_MC_tau_Pa','--'),('Table4 fit at1300 (extrapolation)',dp,'table4_1300_extrapolated_MC_tau_Pa',':')]:
            ax.plot(x,np.array([r[k] for r in rows])/1000,style,label=name,lw=1.5)
        ax.set(title=title,xlabel='Normal stress (kPa)',ylabel='Shear stress (kPa)');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.subplots_adjust(bottom=.22,wspace=.25);fig.text(.5,.035,'Idealized shear vs literature-derived envelopes; not raw experimental curve validation.',ha='center',fontsize=9)
    fig.savefig(HERE/'strength_comparison.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5.4));names=['surface_elastic','duricrust_elastic','sand_elastic_1300']
    titles=['Surface trial','Duricrust proxy','1m sand reference']
    for i,name in enumerate(names):
        vs=next(x['predicted'] for x in waves if x['specimen']==name and x['quantity']=='Vs')
        ax.barh(i,vs,color=['#9e6e3e','#347e91','#4769ad'][i]);ax.text(vs+1,i,f'{vs:.2f}',va='center')
    ax.axvspan(56,74,alpha=.12,color='#339966',label='Brinkman shallow effective reported range56–74')
    ax.axvline(63,color='#339966',ls='--',label='Brinkman central63 m/s');ax.axvline(94.6,color='#5665aa',ls=':',label='Murdoch1m prediction94.6 m/s')
    ax.set(yticks=range(3),yticklabels=titles,xlabel='Vs (m/s)',xlim=(0,120));ax.grid(axis='x',alpha=.2)
    handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,fontsize=7.5,loc='lower center',bbox_to_anchor=(.55,.07))
    fig.subplots_adjust(left=.23,bottom=.32);fig.text(.5,.02,'Different depth/measurement volumes: comparison is a source-state check, not a common test.',ha='center',fontsize=8)
    fig.savefig(HERE/'wave_speed_comparison.png',dpi=160);plt.close(fig)


def report(outcome,criteria,inputs,numerical,waves,strength,stats,behaviour,assessments,sources):
    esc=lambda x:html.escape(str(x))
    def fmt(x):
        if x is None:return '산정 불가/해당 없음'
        if isinstance(x,bool):return '예' if x else '아니오'
        if isinstance(x,float):return f'{x:.6g}'
        return esc(x)
    def table(rows,fields):
        h='<div class="scroll"><table><thead><tr>'+''.join(f'<th>{esc(label)}</th>' for key,label in fields)+'</tr></thead><tbody>'
        for row in rows:h+='<tr>'+''.join('<td>'+fmt(row.get(key))+'</td>' for key,label in fields)+'</tr>'
        return h+'</tbody></table></div>'
    out=['''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>시편 문헌 비교·오차 평가 2026-11-24</title><style>body{font:16px/1.65 -apple-system,BlinkMacSystemFont,"Apple SD Gothic Neo",sans-serif;max-width:1400px;padding:30px;margin:auto;color:#223042}h1{font-size:30px}h2{margin-top:42px}a{color:#126594}p{max-width:1150px}.note{padding:18px 22px;background:#fff1dc;border-left:5px solid #b47720}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:10px;border-bottom:1px solid #dae3e9;text-align:left;vertical-align:top}th{background:#edf3f7}td:first-child{min-width:140px;overflow-wrap:anywhere}nav a{display:inline-block;margin:5px 18px 5px 0}img{max-width:1050px;width:100%;height:auto}.muted{color:#617283;font-size:13px}pre{background:#f3f6f8;padding:15px;white-space:pre-wrap}summary{cursor:pointer}details{margin-top:20px}.stats{display:flex;gap:20px;flex-wrap:wrap}.stats p{background:#edf3f7;padding:15px 25px}.stats b{font-size:27px}</style><h1>시편 문헌 비교와 오차 평가</h1><p>2026-11-24 · 기존 28개 시편 / 추가 문헌 압력 전단16개 · N–m–kg–s–Pa</p>
<p class="note"><strong>이전 28/28 통과는 수치 구현 검사였습니다.</strong> 이번 평가는 입력 출처 일치, 구성식·평형 오차, 문헌 강도 포락선 차이를 분리했습니다. 독립된 화성 시편의 실험곡선을 재현·검증한 재료는 아직 없습니다. 실험 기준값이 없는 항목은 0%가 아니라 <strong>산정 불가</strong>로 표시했습니다.</p>''']
    out.append(f'<div class="stats"><p><b>28/28</b><br>기존 수치 검사 통과</p><p><b>{outcome["additional_numerical_pass"]}/16</b><br>추가 저압 전단 수치 검사</p><p><b>0</b><br>독립 화성 실험 검증 완료</p></div>')
    out.append('''<nav><a href="#criteria">평가기준</a><a href="#numerical">수치 오차</a><a href="#strength">문헌 강도 비교</a><a href="#wave">물성·파속</a><a href="#behaviour">거동 결과</a><a href="#assessment">후보별 판단</a></nav>
<p><a href="criteria.csv">기준 CSV</a> · <a href="numerical_errors.csv">수치 오차 전체</a> · <a href="literature_input_comparison.csv">입력–문헌 대조</a> · <a href="strength_envelope_comparison.csv">강도 비교</a> · <a href="strength_statistics.csv">강도 오차 통계</a> · <a href="wave_speed_comparison.csv">파속 비교</a> · <a href="behaviour_metrics.csv">거동 수치</a> · <a href="specimen_assessments.csv">후보별 판단</a> · <a href="evaluation.json">종합 JSON</a></p>
<h2 id="criteria">평가기준을 어떻게 정했는가</h2><p><strong>수치 오차 기준은 프로젝트에서 정한 허용오차</strong>이며 논문에서 정한 합격기준이나 ASTM 재료 허용오차가 아닙니다. 문헌의 보고범위는 별도로 사용합니다. 강도 포락선의 ±5%는 이번 비교를 위한 잠정 목표이며 실험 불확도·현장 안전 기준이 아닙니다.</p>
<p>상대오차 = |해석값−기준값|/|기준값|×100. 부호 있는 차이는 (해석값/기준값−1)×100이며 +는 과대입니다. 기준값이0이면 상대오차를 정의하지 않고 Pa/N 단위 절대오차를 사용했습니다. 힘 평형은 이력 최대력으로 정규화해0 부근 분모 문제를 피했습니다.</p>''')
    out.append(table(criteria,[('id','ID'),('metric','항목'),('formula','계산식'),('tolerance_pct','허용 상대오차 %'),('tolerance_absolute','허용 절대오차'),('scope','평가 범위'),('basis','근거/한계')]))
    out.append('<h2 id="numerical">구성식·경계·힘 평형 오차</h2><p>50개 기존 하중 이력과22개 질량 검사 결과를 다시 계산했습니다. 아래의 극히 작은 값은 구성식/입력과의 일치이며 실제 토양의 예측 오차가 아닙니다. DP 삼축압축 이론값 qᵧ = [2σ₃sinφ + 2c cosφ]/[1−sinφ]를 사용했습니다.</p>')
    wr=[dict(quantity=k,**v) for k,v in outcome['worst_numerical_errors'].items()]
    out.append(table(wr,[('quantity','항목'),('max_error_pct','최대 오차 %'),('specimen','해당 시편'),('case','조건')]))
    dp=[x for x in numerical if x['quantity']=='q_yield'];out.append('<details><summary>DP 6개 × 구속압2개: 항복 기준값·결과·오차</summary>'+table(dp,[('specimen','시편'),('case','조건'),('reference','이론 q Pa'),('computed','해석 q Pa'),('absolute_difference','차이 Pa'),('error_pct','오차 %'),('passed','기준 충족')])+'</details>')
    out.append('''<p class="muted">제하 시 잔류축변형 2qᵧ/E는 이번 변위 제어 이력에서 의도한 값입니다. 실제 실험의 잔류침하 검증이 아닙니다. PDMY02의 γ=0에서 남은 τ는 구속 상태의 잔류응력이며 자유 제하 잔류변형과 구분했습니다.</p>
<h2 id="strength">문헌 강도와의 차이: MGS-1 / MGS-1S</h2><p>같은 저장 재료를 문헌 정상응력 97/250/500/680Pa에서 추가 단순전단했습니다. DP는 포락선 부근 평탄 응답을 얻었습니다. PDMY02는 γ0.2% 시점이며 일부는 여전히 경화하므로 극한강도로 부르지 않습니다. 정상응력 범위를 맞췄어도 논문10.17×10.17cm 직접전단 상자와 이번1m 균질요소의 경계조건·응력경로는 다릅니다.</p>
<p><strong>① 기존 그림 기반 기준:</strong> τref = 280Pa + σₙtanφ, φ=35.5°/34°. 입력 선정에 사용한 근사값이므로 독립 검증 데이터는 아닙니다.<br><strong>② Table4·식2 재확인:</strong> c=A₁exp(kρ), ρ는g/cm³. MGS-1: A₁4.96Pa,k2.97; MGS-1S: A₁5.28Pa,k2.93. ρ1.30에서 c235.664/238.156Pa입니다. 기존280Pa는 각각 +18.81/+17.57% 높습니다. 다만 원 논문 최소밀도1.32/1.31보다 낮아 이 값은 <strong>작은 외삽</strong>입니다. 최소 시험밀도에서의 적합식 c250.086/245.237Pa와의 비교도 CSV에 병기했습니다. φ는 인접 밀도의 값으로 유지했습니다.</p>
<p>A₁·k의 공분산과 원 시험점 잔차가 없으므로 이 적합식에서 엄밀한 신뢰구간을 만들지 않았습니다. 그림 근사값–적합식 차이를 시험 오차로 단정하지 않습니다.</p><img src="strength_comparison.png" alt="MGS 전단 응답과 문헌 포락선 비교">''')
    out.append(table(stats,[('specimen','시편'),('graph_reference_MAPE_pct','그림 포락선 MAPE %'),('graph_reference_RMSE_Pa','그림 포락선 RMSE Pa'),('graph_reference_NRMSE_pct','NRMSE %'),('table4_extrapolated_MAPE_pct','표4 외삽 기준 MAPE %')]))
    out.append('<p class="muted">MAPE는 정상응력4개 지점의 평균 절대백분율 차이. NRMSE=RMSE/기준 전단강도 평균×100. 이 RMSE는 응력–변형률 원 시험곡선의 오차가 아닙니다. 모두 잠정 ±5% 목표를 넘으므로 현재 강도 일치를 주장하지 않습니다.</p>')
    out.append(table(strength,[('specimen','시편'),('normal_stress_Pa','σn Pa'),('computed_tau_Pa','τ 해석 Pa'),('approx_graph_MC_tau_Pa','그림 기준 Pa'),('graph_envelope_signed_difference_pct','그림 대비 %'),('table4_1300_extrapolated_MC_tau_Pa','표4 외삽 기준 Pa'),('table4_1300_signed_difference_pct','표4 대비 %'),('tail_variation_pct','마지막20점 변화 %')]))
    out.append('''<p><strong>DP 차이의 원인:</strong> 현재의 삼축압축 일치 DP 매핑에서, 무팽창·충분한 단순전단 뒤 세 정상응력이−p로 같아지면 항복식으로부터 τDP = [2√3 cosφ/(3−sinφ)]·(c+p tanφ)가 됩니다. 이 계수는 MGS-1에서1.165702, MGS-1S에서1.176607입니다. 즉 근사 MC 포락선 대비 <strong>+16.57/+17.66%는 수치 오차가 아니라 이 매핑과 응력경로가 만드는 구조적인 강도 차이</strong>입니다. 실제 추가 해석은 이 DP 단순전단 이론과 최대 약0.000105% 차이로 맞았습니다. 기존 c의 근사값 차이는 별도의 항목입니다.</p>
<p>따라서 직접전단 포락선을 목표로 DP를 보정하면 기존 삼축압축 매칭도 다시 점검해야 합니다. PDMY02의 수축·팽창/반복 계수를 유지한 채 초기 강성과φ만 바꾼 초안도 이 시험의 문헌 포락선과 일치하지 않습니다. 이번에는 재료값을 자동 수정하지 않았습니다.</p>
<h2 id="wave">탄성 물성과 파속 대조</h2><p>Vs=√(G/ρ), Vp=√[(K+4G/3)/ρ]. Brinkman 천부 유효 보고범위와 얇은 개별층의 대용값을 구분했습니다. 표면 대용값은 보고범위보다 느리고, 경화층 대용값은 보고범위 안입니다. 이것만으로 층별 실측과 일치한다고 판정할 수 없습니다. 모래의 1m 파속 일치는 원래 그 예측 파속으로 물성을 만들었기 때문에 반올림 일관성 검사입니다.</p><img src="wave_speed_comparison.png" alt="문헌의 깊이/범위를 구분한 Vs 비교">''')
    out.append(table(waves,[('specimen','시편'),('quantity','항목'),('predicted','계산 m/s'),('reference','문헌 중심/참조 m/s'),('signed_difference_pct','차이 %'),('reported_low','보고 하한'),('reported_high','보고 상한'),('within_reported_interval','범위 내'),('context','비교 조건')]))
    out.append('''<p>경화층 E11.48MPa·ν.28에서 계산한 G4.484MPa/K8.697MPa는 Brinkman 표의 독립 최빈값 G4.47/K7.79와 각각 약 +0.322/+11.643% 차이입니다. K는 보고범위6.24–9.39MPa 안입니다. 이 차이를 구성식 오차로 계산하지 않습니다. 밀도1025는950–1100 추정범위의 선택값입니다. 파쇄 현무암 E13GPa/ν.28/ρ2760은 다른 참조의 Vp3000/Vs1700과 동시 일치하지 않습니다.</p>
<details><summary>28개 모든 입력의 문헌값·차이·추정/가정 구분</summary>''')
    out.append(table(inputs,[('specimen','시편'),('parameter','입력'),('input_value','값 SI'),('reference_value','동일 기준 참조'),('signed_difference_pct','차이 %'),('comparison_kind','근거 구분'),('within_reported_interval','범위 내'),('source','출처'),('locator','위치'),('note','주의')])+'</details>')
    out.append('''<h2 id="behaviour">모든 시편의 거동 수치</h2><p>탄성은 압축 후 원점 복귀, DP는 항복·제하·잔류축변형, PDMY02는 최대 전단응력·최종 구속응력·수직변위, 스프링은 최대력과 용량 동원율을 정리했습니다. 최대력/Fult는 유한 변위의 동원율이며 (1−동원율)을 문헌 예측 오차로 부르지 않습니다. 전체 수치는 behaviour_metrics.csv에 있습니다.</p>''')
    out.append(table(behaviour,[('specimen','시편'),('case','조건'),('initial_E_MPa','초기 E MPa'),('peak_q_kPa','최대 q kPa'),('final_axial_strain','최종 축변형'),('initial_G_MPa','PDMY 초기 G MPa'),('peak_abs_tau_kPa','최대 |τ| kPa'),('tau_at_zero_final_kPa','γ0 최종 τ kPa'),('bisections','분할 횟수'),('positive_force_kN','양 최대력 kN'),('negative_force_kN','음 최대력 kN'),('mobilized_capacity_pct','용량 동원율 %'),('nominal_half_capacity_difference_pct','명목50%점 차이 %')]))
    out.append('''<p class="muted">QzSimple1의 z50는 문서상 단조재하50% 용량 기준이나, 저장 이력의 최초 압축 z=−z50에서 실제 값은 명목치보다8.93% 낮습니다. QzSimple2는11.83%, PySimple2는1.62%, TzSimple2는0.000405% 낮습니다. 이는 명목 배면곡선 진단이며 서로 다른 Simple1/2를 같은 입력만으로 동일 곡선이라 가정할 수 없음을 보여줍니다. 설계 접촉곡선 보정 전 확인 항목으로 남깁니다. 합성 Hysteretic/EPPGap에는 독립 실험 기준이 없어 문헌곡선 RMSE를 산정하지 않습니다.</p>
<h2 id="assessment">후보별 사용 판단</h2>''')
    out.append(table(assessments,[('specimen','시편'),('base_numerical_pass','수치 검사'),('decision','현 단계 판단'),('reason','근거/남은 검증'),('experimental_curve_RMSE','독립 실험곡선 RMSE')]))
    out.append('<h2>재현과 출처</h2><pre>cd '+str(HERE)+'\n/Users/sydney/anaconda3/envs/cap/bin/python run_additional_tests.py\n/Users/sydney/anaconda3/envs/cap/bin/python evaluate.py</pre><p><a href="../index.html">기존 시편 모음</a> · <a href="input_manifest.json">평가 입력 SHA-256</a> · <a href="source_manifest.json">보관 문헌 SHA-256</a> · <a href="dotson_table4_extracted.json">새로 대조한 Table4 계수</a> · <a href="additional_tests/summary.json">추가 시험 상세</a></p><ul>')
    for source in sources:out.append(f'<li><a href="{esc(os.path.relpath(source["file"],HERE))}">{source["id"]} 보관 원문</a> · <a href="{esc(source["url"])}">원 출처</a></li>')
    out.append('</ul><p class="muted">모든 비교 원문은 기존 docs/references에 보관되어 있습니다. 수치·물성 파일은 변경하지 않았습니다. 비교조건이 없는 자료에 임의 오차를 채우지 않았습니다.</p></html>')
    (HERE/'report.html').write_text(''.join(out))
    (HERE/'README.txt').write_text('2026-11-24 시편 평가\nreport.html을 여세요.\n기존 28개 시편의 수치 구현검사와 문헌 거동 일치는 별개입니다.\ncriteria.csv: 기준·계산식·허용오차\nnumerical_errors.csv: 구성식/질량/평형 수치 오차\nliterature_input_comparison.csv: 모든 입력의 문헌 근거·차이·가정\nstrength_envelope_comparison.csv: 문헌 정상응력의 추가 전단16개\nstrength_statistics.csv: 포락선 MAPE/RMSE (원 실험곡선 오차 아님)\nbehaviour_metrics.csv: 모든 시편의 거동 수치\nspecimen_assessments.csv: 28개 후보별 판단\n독립 실험곡선이 없으면 RMSE는 빈칸/산정 불가입니다.\n')


if __name__=='__main__':main()
