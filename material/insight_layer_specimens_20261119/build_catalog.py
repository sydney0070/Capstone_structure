"""Extracted source values -> explicitly labelled specimen parameter records (SI)."""
import copy
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LIBRARY = ROOT / 'docs/references'


def evidence(source, locator, kind, note=''):
    return dict(source=source, locator=locator, kind=kind, note=note)


def build():
    records = []
    sources = {
        'Spohn2022': ('01_ground/spohn_2022_insight_hp3.pdf', 'https://doi.org/10.1007/s11214-022-00941-z'),
        'Brinkman2022': ('01_ground/G07_Brinkman_et_al_2022_In_Situ_Regolith_Seismic_Velocity_Measurement_at_the_InSight_Lan.pdf', 'https://doi.org/10.1029/2022JE007229'),
        'Murdoch2018': ('01_ground/G14_Murdoch_2018_InSight_Lander_Flexible_Modes_AAM.pdf', 'https://doi.org/10.1007/s11214-018-0553-y'),
        'Dotson2024': ('01_ground/R13_Regolith_Cohesion_Shear_Strength_2024.pdf', 'https://arxiv.org/abs/2403.11029'),
        'Morgan2018': ('01_ground/G12_Morgan_2018_InSight_Regolith_Properties.pdf', 'https://elib.dlr.de/122371/1/Revised_InSight%20Landing%20Site%20Regolith%20Physical%20Properties%20Complete%202018-6-6.pdf'),
        'Hobiger2021': ('01_ground/G15_Hobiger_2021_InSight_Shallow_Structure.pdf', 'https://doi.org/10.1038/s41467-021-26957-7'),
        'Kenda2020': ('01_ground/G16_Kenda_2020_Compliance.pdf', 'https://doi.org/10.1029/2020JE006387'),
        'Shallow2008': ('04_SSI_and_numerics/R20_ShallowFoundationGen_2008.pdf', 'https://opensees.berkeley.edu/OpenSees/manuals/usermanual/documents/ShallowDocumentation.pdf'),
        'PDMY02': ('05_official_model_docs/PressureDependentMultiYield02.html', 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html'),
        'DP_command': ('05_official_model_docs/DruckerPrager.html', 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html'),
        'Hysteretic_command': ('05_official_model_docs/D31_Hysteretic_20261119.html', 'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Hysteretic.html'),
        'ElasticPPGap_command': ('05_official_model_docs/D16_ElasticPPGap.html', 'https://opensees.berkeley.edu/OpenSees/manuals/usermanual/170.htm'),
    }

    def solid(name, layer, E, nu, rho, E_src, rho_src, law='ElasticIsotropic', phi=None, c=None, strength=None):
        p = dict(E_Pa=E, nu=nu, G_Pa=E/(2*(1+nu)), K_Pa=E/(3*(1-2*nu)), density_kg_m3=rho)
        ev = dict(E_Pa=E_src, nu=E_src, density_kg_m3=rho_src)
        for k in ('G_Pa', 'K_Pa'):
            ev[k] = evidence('isotropic_elasticity', 'G=E/[2(1+nu)], K=E/[3(1-2nu)]', 'derived')
        if law == 'DruckerPrager':
            p.update(phi_deg=phi, c_Pa=c, psi_deg=0., hardening=0.)
            ev.update(phi_deg=strength, c_Pa=strength,
                      psi_deg=evidence('specimen_assumption', 'non-associated flow', 'assumed', 'Zero dilation; not measured'),
                      hardening=evidence('specimen_assumption', 'perfect plasticity', 'assumed', 'No measured hardening/softening curve'))
        r = dict(id=name, group='first_comparison', layer=layer, family=law, parameters=p,
                 evidence=ev, mars_calibrated=False, geometry=dict(type='cube', side_m=1., area_m2=1., volume_m3=1.),
                 status='SOURCE_BASED_INITIAL_TRIAL',
                 notes=['Unit cube is a material test, not the actual thickness of this geological unit.',
                        'No selfweight; confinement is an imposed laboratory test condition.'])
        records.append(r)
        return r

    surf = evidence('Kenda2020', 'Table 1, surface regolith', 'transferred', 'Expected effective surface properties, not dust-only measurements')
    sand = evidence('Murdoch2018', 'Table 2; existing 1300 kg/m3 candidate', 'transferred', 'Rounded E=28.4 MPa, nu=.22 from previous cross-source initial set')
    crust = evidence('Brinkman2022', 'Table 4, PDF p.13', 'transferred', 'Select E=11.48 MPa and nu=.28; effective shallow volume, not isolated crust. Derive G/K consistently; do not mix independent modal estimates.')
    block = evidence('Kenda2020', 'Table 1, coarse ejecta', 'published_reference')
    rock = evidence('Kenda2020', 'Table 1, fractured bedrock', 'published_reference')
    generic = evidence('specimen_assumption', 'unbonded granular sensitivity', 'assumed', 'c=0 and phi are trial envelope choices; no layer-specific shear test')
    solid('surface_elastic', 'local_surface_0_1cm', 7.5e6, .22, 1200., surf,
          evidence('Spohn2022', 'Section 6.7, Figure 20', 'inferred', 'Approximate loose surface mixture density'))
    r = copy.deepcopy(records[-1])
    solid('surface_dp_trial', r['layer'], 7.5e6, .22, 1200., surf, r['evidence']['density_kg_m3'], 'DruckerPrager', 30., 0., generic)
    rho_crust = evidence('Spohn2022', 'Section 6.7, range 950–1100 kg/m3', 'assumed_within_inferred_range', 'Midpoint 1025 chosen for this specimen, not directly measured')
    solid('duricrust_elastic', 'local_duricrust_1_19cm_approx', 11.48e6, .28, 1025., crust, rho_crust)
    solid('duricrust_dp', 'local_duricrust_1_19cm_approx', 11.48e6, .28, 1025., crust, rho_crust,
          'DruckerPrager', 30., 5800., evidence('Spohn2022', 'Section 6.4, slope analysis; phi=30 deg and rho1200 assumed', 'inferred_conditional', 'c=5.8 kPa is inferred conditional on phi and rho1200; transferred to rho1025 without new calibration. Not independent measured c and phi.'))
    rho_sand = evidence('Spohn2022', 'Section 6.7, sand 1300–1500 kg/m3', 'selected_inferred_reference', '1300 is also the previously saved target density')
    solid('sand_elastic_1300', 'local_sand_19_31cm_and_metre_scale_sandy_reference', 28.4e6, .22, 1300., sand, rho_sand)
    old = json.loads((ROOT/'ground/ground_modeling/ground_hybrid/configs/materials.json').read_text())['materials']
    for code, phi in [('mgs1',35.5), ('mgs1s',34.)]:
        saved = old[f'mars_insight_sand_{code}_dp_1300']
        r = solid(f'sand_{code}_dp_1300', 'sandy_regolith_1300', saved['E_Pa'], saved['nu'], saved['density_kg_m3'], sand, rho_sand,
                  'DruckerPrager', saved['phi_deg'], saved['c_Pa'], evidence('Dotson2024', 'Figures 5–6; previously digitized set', 'simulant_transfer_approximate', 'c≈280 Pa, phi near rho1310–1320 transferred to rho1300. Direct-shear normal stress .097–.68 kPa, not foundation pressure validation'))
        r['notes'].extend(saved['assumptions'])
    rho_gravel = evidence('Spohn2022', 'Section 6.7, gravel mixture estimate', 'inferred', 'Approximate 1600 from 15 vol-% rock and 1300 kg/m3 sand')
    matrix = evidence('Murdoch2018', 'sand reference reused as matrix proxy', 'proxy_assumption', 'No gravel-layer E measured; penetration resistance 4.9 MPa is NOT E')
    solid('local_gravel_elastic_proxy', 'local_gravel_below_31cm_penetrated_interval_only', 28.4e6, .22, 1600., matrix, rho_gravel)
    solid('local_gravel_dp_trial', 'local_gravel_below_31cm_penetrated_interval_only', 28.4e6, .22, 1600., matrix, rho_gravel, 'DruckerPrager', 35., 0., generic)
    rho_block = evidence('existing_project_reference', 'mars_insight_blocky_ejecta', 'transferred_model_density', '1700 carried over from published reference density scenario, not core measurement')
    solid('blocky_ejecta_elastic', 'metre_scale_blocky_ejecta_3_13m_reference', 520e6, .24, 1700., block, rho_block)
    solid('blocky_ejecta_dp_trial', 'metre_scale_blocky_ejecta_3_13m_reference', 520e6, .24, 1700., block, rho_block, 'DruckerPrager', 35., 0., generic)
    solid('fractured_basalt_elastic', 'fractured_basalt_13_30m_reference', 13e9, .28, 2760., rock,
          evidence('Morgan2018', 'Section 4.1, terrestrial fractured basalt analogy', 'analogue_transfer', 'Effective continuum; no joint failure strength inferred'))
    rho, vp, vs = 1700., 750., 417.
    G, K = rho*vs**2, rho*(vp**2-4*vs**2/3)
    E, nu = 9*K*G/(3*K+G), (3*K-2*G)/(2*(3*K+G))
    r = solid('sediment_elastic_prior', 'deep_low_velocity_sediment_hypothesis', E, nu, rho,
              evidence('Hobiger2021', 'Table 2 modified prior, sediment Vp750 Vs417 m/s', 'derived_from_prior', 'Prior values, not final inversion point or actual 30–75m measured stiffness; density assumed'),
              evidence('specimen_assumption', 'velocity-to-modulus conversion', 'assumed', '1700 kg/m3 requires density sensitivity'))
    r['parameters'].update(Vp_m_s=vp, Vs_m_s=vs)
    for key in ('Vp_m_s','Vs_m_s'):
        r['evidence'][key] = evidence('Hobiger2021', 'Table 2 modified prior', 'published_prior')

    # Official table: G/K are labelled at 80 kPa although refPress row says 101.
    # Resolve by scaling both to 101 using the documented exponent .5.
    table = [(30,1700,60e6,160e6,31,31,.087,.18,0,0,.85),
             (40,1800,90e6,220e6,32,26,.067,.23,.06,.27,.77),
             (50,1900,100e6,233e6,33.5,25.5,.045,.15,.06,.15,.70),
             (60,2000,110e6,240e6,35,26,.028,.05,.10,.05,.65),
             (75,2100,130e6,260e6,36.5,26,.013,0,.30,0,.55)]
    for dr,rho,G,K,phi,PT,c1,c3,d1,d3,e in table:
        factor=(101/80)**.5
        p=dict(density_kg_m3=rho,G_Pa=G*factor,K_Pa=K*factor,phi_deg=phi,peak_shear_strain=.1,
               p_ref_Pa=101000.,pressure_exponent=.5,PTAngle_deg=PT,contrac1=c1,contrac3=c3,dilat1=d1,dilat3=d3,
               n_surfaces=20,contrac2=5.,dilat2=3.,liquefac1=1.,liquefac2=0.,e=e,cs1=.9,cs2=.02,cs3=.7,pa_Pa=101000.,numerical_c_Pa=300.)
        ev={k:evidence('PDMY02','suggested parameters table / documented defaults','official_earth_reference') for k in p}
        for key in ('G_Pa','K_Pa'):
            ev[key]=evidence('PDMY02','suggested G/K at 80 kPa; scaled (101/80)^.5', 'derived', 'RefPress=101kPa matches gamma_peak reference; all SI units')
        records.append(dict(id=f'pdmy02_earth_dr{dr}_reference',group='followup_reference',layer='earth_sand_reference_only',family='PressureDependMultiYield02',
                            parameters=p,evidence=ev,mars_calibrated=False,status='OFFICIAL_EARTH_REFERENCE_NOT_MARS',
                            geometry=dict(type='cube',side_m=1.,area_m2=1.,volume_m3=1.),
                            notes=['Official rho is saturated bulk density. Dr is the official table label, not inferred Mars relative density.',
                                   'Dry/effective-stress mechanical fixture only; no pore pressure, u-p element or liquefaction claim.']))
    ref = copy.deepcopy(next(x for x in records if x['id']=='pdmy02_earth_dr50_reference'))
    for base_id in ('sand_mgs1_dp_1300','sand_mgs1s_dp_1300','local_gravel_dp_trial','blocky_ejecta_dp_trial'):
        base=next(x for x in records if x['id']==base_id)
        r=copy.deepcopy(ref)
        r.update(id=base_id.replace('_dp','_pdmy02'),group='followup_mars_trial',layer=base['layer'],status='UNCALIBRATED_MARS_SENSITIVITY_SEED')
        for key in ('G_Pa','K_Pa','density_kg_m3','phi_deg'):
            r['parameters'][key]=base['parameters'][key]
            r['evidence'][key]=copy.deepcopy(base['evidence'][key])
        # Source stiffness is imposed at a chosen 5kPa mean reference, NOT a measured pressure law.
        r['parameters']['p_ref_Pa']=5000.
        r['evidence']['p_ref_Pa']=evidence('specimen_assumption','5 kPa reference pressure','assumed','Anchor copied initial stiffness at 5 kPa; source did not identify this pressure dependence')
        copied=set(r['parameters'])-{'G_Pa','K_Pa','density_kg_m3','phi_deg','p_ref_Pa'}
        for key in copied:
            r['evidence'][key]=evidence('PDMY02','Dr50 defaults transferred to uncalibrated trial','proxy_assumption','Not identified by Mars observations. No claim of Mars Dr50 or matching physical c.')
        r['notes']=['Operational sensitivity seed only. All cyclic/critical-state parameters and e=.70 are transferred assumptions.',
                    'e=.70 is not density-derived. Grain density and grading need independent measurement.',
                    'numerical_c=300 Pa is model regularization, not the DP c=280 Pa calibration.',
                    'Low-pressure and cyclic experimental calibration required before ground use.']
        records.append(r)

    def spring(name, family, layer, parameters, ev, scale):
        records.append(dict(id=name,group='uniaxial_comparison',layer=layer,family=family,parameters=parameters,
                            evidence={k:copy.deepcopy(ev) for k in parameters},mars_calibrated=False,
                            status='LITERATURE_REFERENCE_OR_SYNTHETIC_MECHANICS',geometry=dict(type='zeroLength',direction=1 if family in ('PySimple2','TzSimple2') else 3),
                            displacement_scale_m=scale,
                            notes=['Force N versus relative displacement m, not continuum stress versus strain.',
                                   'Independent spring fixture; no automatic normal/shear coupling or tributary ring-area calibration.']))
    ss=evidence('Shallow2008','PDF p.17, first example, center/end spring 101', 'earth_example','Original qzType=1. Static damping set to zero instead of original .05; not a Mars force capacity.')
    for v in (2,1):
        ev=ss if v==2 else evidence('Shallow2008','same force and Z50 as spring101, transferred to QzSimple1','comparison_transfer','qzType=1 and suction=.1; damping zero. Same inputs do not guarantee identical backbones.')
        spring(f'qzsimple{v}_reference',f'QzSimple{v}','footing_vertical_equivalent',dict(soil_type=1,ultimate_force_N=38020.,z50_m=.00274347,suction=.1,damping=0.),ev,.00274347)
    spring('pysimple2_reference','PySimple2','embedded_footing_side',dict(soil_type=1,ultimate_force_N=102000.,y50_m=.0145067,drag=.1,damping=0.),
           evidence('Shallow2008','PDF p.17 spring105','earth_example','Static damping zero; clay-type reference, not Mars calibration'),.0145067)
    spring('tzsimple2_reference','TzSimple2','footing_base_tangential',dict(soil_type=1,ultimate_force_N=100000.,z50_m=.0142222,damping=0.),
           evidence('Shallow2008','PDF p.17 spring106','earth_example','Use documented one optional damping argument, set zero for static test'),.0142222)
    # Conditional mechanics benchmarks: literature E,c,phi provide scale, not measured cyclic curves.
    E,A,h,p0,c,phi=11.48e6,1.,1.,500.,5800.,math.radians(30)
    qy=2*p0*math.sin(phi)/(1-math.sin(phi))+2*c*math.cos(phi)/(1-math.sin(phi))
    fy=qy*A; k=E*A/h; uy=fy/k
    ev=evidence('Spohn2022+Brinkman2022+specimen_assumption','c5800 phi30; E11.48MPa; p=.5kPa; A=h=1', 'synthetic_mechanics', 'k=EA/h; fy=A[2p sin(phi)+2c cos(phi)]/[1-sin(phi)]. Not measured joint strength or footing capacity.')
    spring('elasticppgap_compression_benchmark','ElasticPPGap','conditional_normal_contact_only',dict(stiffness_N_m=k,yield_force_N=-fy,gap_m=0.,hardening_ratio=0.,damage='damage',area_m2=A,height_m=h),ev,uy)
    ev=evidence('Spohn2022+Brinkman2022+specimen_assumption','synthetic symmetric DP-scale envelope','synthetic_mechanics','Backbone points (fy,uy), (1.01fy,4uy), (1.02fy,8uy) are a mechanics benchmark; NO measured crust tension/cyclic/pinching calibration. Symmetry is not crust tensile strength.')
    spring('hysteretic_crust_scale_benchmark','Hysteretic','conditional_duricrust_curve_fitting_only',
           dict(positive_backbone=[[fy,uy],[1.01*fy,4*uy],[1.02*fy,8*uy]],
                negative_backbone=[[-fy,-uy],[-1.01*fy,-4*uy],[-1.02*fy,-8*uy]],pinchX=1.,pinchY=1.,damage1=0.,damage2=0.,beta=0.),ev,uy)

    payload=dict(started_on='2026-11-19',completed_on='2026-11-23',units='N, m, kg, s, Pa',scope='independent material specimens, not field validation or ground analysis',
                 source_policy='Each parameter explicitly separates reference, inference, transfer and assumption.',materials=records)
    (HERE/'materials.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    with (HERE/'parameters_and_sources.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['specimen','group','family','layer','parameter','value_SI','source','locator','evidence_kind','note'])
        for r in records:
            for key,val in r['parameters'].items():
                ev=r['evidence'][key];w.writerow([r['id'],r['group'],r['family'],r['layer'],key,json.dumps(val,ensure_ascii=False),ev['source'],ev['locator'],ev['kind'],ev['note']])
    cases=HERE/'cases';cases.mkdir(exist_ok=True)
    for r in records:
        (cases/(r['id']+'.json')).write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
        (cases/(r['id']+'.py')).write_text('"""Run this independent material specimen; SI units."""\nfrom pathlib import Path\nimport sys\nsys.path.insert(0, str(Path(__file__).resolve().parents[1]))\nfrom specimens import main\nif __name__ == "__main__":\n    main(["--case", '+repr(r['id'])+'])\n')
    src=[]
    for name,(file,url) in sources.items():
        p=LIBRARY/file
        src.append(dict(id=name,file=str(p),url=url,exists=p.is_file(),sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None))
    (HERE/'source_manifest.json').write_text(json.dumps(src,ensure_ascii=False,indent=2)+'\n')
    print(f'Created {len(records)} independent specimen records')


if __name__=='__main__':
    build()
