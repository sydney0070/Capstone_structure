"""Rebuild SI material selections from identified evidence; no silent optimization."""
import hashlib,json,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'material_tuning_20261008'/'material_calibration'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def isotropic(E,nu):return dict(E_Pa=E,nu=nu,G_Pa=E/(2*(1+nu)),K_Pa=E/(3*(1-2*nu)))
sources={
 's355_coupon':dict(url='https://doi.org/10.5281/zenodo.6965147',scope='S355 J2+M HEM320 web selected room-temperature prefracture prefixes |strain| <= 0.02',local=str(OLD/'s355_hem320_web'/'reviewed_config.json')),
 'hydro6061':dict(url='https://www.hydro.com/globalassets/01-products--services/extruded-profiles/americas/ena-resources/alloy-data-sheets/hydro_2019_data_sheet_6061.pdf',location='p2 tensile properties, T6/T6511 extruded products',scope='minimum strength; no weld HAZ reduction or cyclic data'),
 'nist6061':dict(url='https://www.nist.gov/mml/acmd/aluminum-6061-t6-uns-aa96061',location='Young modulus polynomial, 2-295 K; fit error 1%'),
 'niar_as4':dict(url='https://www.wichita.edu/industry_and_defense/NIAR/Research/hexcel-8552/AS4-Unitape-2.pdf',location='CAM-RP-2010-002 RevA, Table2-1 PDFp30 RTD normalized means; nominal ply .0074in; density specimens',scope='dry 75F; mean qualification results, NOT design allowables; downloaded PDF blocked, table accessed through indexed primary PDF'),
 'sintered2016':dict(url='https://www.hou.usra.edu/meetings/leag2016/eposter/5082.pdf',local=str(ROOT/'research'/'sintered_jsc1a_leag2016.pdf'),scope='Indyk/Benaroya 2016 JSC-1A sintered small cylinders, same-source groups; no tensile data'),
 'lasp_jsc1a':dict(local=str(OLD/'jsc1a_lasp_review_20261019'/'triaxial_observations.json'),scope='LASP/JSC-1A Tons1-3 Earth-gravity dry drained triaxial summary; measured sigma3 15-60kPa; simulant not in-situ Moon'),
 'jez1':dict(url='https://www.mdpi.com/2076-3417/11/15/7052',scope='dry tamped JEZ-1 Dr86.28%, direct shear .5mm/min, sigma_n25-75kPa soil,25-100kPa interfaces; simulant not in-situ Mars'),
 'pdmy02':dict(url='https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html',scope='suggested Dr50% sand parameters, illustrative not site calibration; G/K listed at80kPa whereas reference pressure listed101kPa'),
 'dp':dict(url='https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/DruckerPrager.html',scope='triaxial compression matching convention; same phi does not reproduce MC in all Lode-angle paths'),
}
for s in sources.values():
 if 'local' in s:s['sha256']=digest(Path(s['local']))
steel=json.loads((OLD/'s355_hem320_web'/'reviewed_config.json').read_text())['params']
steel['Fy']*=1000;steel['E']*=1000
T=293.;Eal=sum(a*T**i for i,a in enumerate([77.71221,.01030646,-.0002924100,8.993600e-7,-1.070900e-9]))*1e9
M={
 'steel_s355_coupon':dict(role='frame uniaxial cyclic; selected coupon scope',law='Steel02',params=steel,density_kg_m3=7850,nu=.3,status='PREVIOUS_COUPON_CALIBRATION_CURRENT_RUNTIME_REPLAY_IN_RESULTS',sources=['s355_coupon'],assumptions=['density7850 and nu0.30 are conventional assumptions; not identified from coupon']),
 'steel_s355_nominal':dict(role='general steel baseline beam/shell/solid',law='J2PerfectPlastic',**isotropic(200e9,.3),fy_Pa=355e6,density_kg_m3=7850,status='NOMINAL_ASSUMPTION_NOT_GRADE_THICKNESS_CERTIFIED',sources=[],assumptions=['S355 nominal reference; actual steel grade/product/thickness not received; EPP postyield is sensitivity model']),
 'al6061_t6_extrusion':dict(role='frame/skin/solid',law='J2PerfectPlastic',**isotropic(Eal,.33),fy_Pa=240e6,fu_min_Pa=260e6,density_kg_m3=.098*0.45359237/0.0254**3,temperature_K=T,status='LITERATURE_ELASTIC_AND_YIELD_REFERENCE',sources=['hydro6061','nist6061'],assumptions=['nu0.33 assumed; EPP afteryield NOT calibrated; no welding/heat exposure degradation; sheet products require new product reference']),
 'concrete_c30_prototype':dict(role='RC unconfined concrete fiber',law='Concrete02',params=dict(fpc=-30e6,epsc0=-.002,fpcu=-6e6,epsU=-.006,lam=.1,ft=3e6,Ets=1.5e9),**isotropic(30e9,.2),density_kg_m3=2400,status='FIXED_PROTOTYPE_NOT_EXPERIMENTALLY_CALIBRATED',sources=[],assumptions=['C30/SD400 inherited prototype, NOT explicitly chosen grade by user','softening not objective under arbitrary element/integration changes; no confinement enhancement']),
 'rebar_sd400_prototype':dict(role='RC reinforcing fibers',law='Steel02',params=dict(Fy=400e6,E=200e9,b=.01,R0=18,cR1=.925,cR2=.15),density_kg_m3=7850,status='NOMINAL_PROTOTYPE_NOT_UCSD_A6_CALIBRATION',sources=[]),
 'cfrp_as4_8552_ud_rtd':dict(role='orthotropic ply; laminate axes and stacking explicit',law='ElasticOrthotropic',E1_Pa=18.46e6*6894.757293168,E2_Pa=1.34e6*6894.757293168,G12_Pa=.70e6*6894.757293168,nu12=.302,nu23=.35,density_kg_m3=1590,ply_thickness_m=.0074*.0254,
 Xt_mean_Pa=289.47e3*6894.757293168,Xc_mean_Pa=202.70e3*6894.757293168,Yt_mean_Pa=9.27e3*6894.757293168,Yc_mean_Pa=38.85e3*6894.757293168,S_offset02_Pa=8e3*6894.757293168,
 status='ELASTIC_PLY_WITH_MEAN_FIRST_ONSET_SCREENING',sources=['niar_as4'],assumptions=['E3=E2,G13=G12,nu23=.35,G23=E2/(2(1+nu23)) transverse isotropy; out-of-plane constants NOT measured here','no progressive damage/delamination/impact; compression and tension moduli differ in source; single linear E uses tension normalized value']),
}
for label,E,nu,rho,fc,por,nE in [('high_porosity',5475e6,.246,2200,84.3e6,.1178,5),('low_porosity',8358e6,.266,2600,202.3e6,.0144,6)]:
 M['sintered_jsc1a_'+label]=dict(role='structural regolith candidate',law='ElasticIsotropic',**isotropic(E,nu),density_kg_m3=rho,fc_mean_Pa=fc,ft_Pa=None,porosity=por,E_sample_n=nE,compression_sample_n=8,status='ELASTIC_COMPRESSION_SCREEN_ONLY_TENSION_CAPACITY_UNKNOWN',sources=['sintered2016'],assumptions=['mean properties not allowable; no extrapolation to bonded regolith; 2016 source only'])
obs=json.loads((OLD/'jsc1a_lasp_review_20261019'/'triaxial_observations.json').read_text())
fits={}
for den in [1.7,1.8,1.9]:
 rows=[r for r in obs if r['nominal_density_g_cm3']==den]
 x=np.array([r['sigma3_peak_kPa'] for r in rows]);y=np.array([r['sigma1_peak_kPa'] for r in rows])
 slope=float(x@y/(x@x));phi=math.degrees(math.asin((slope-1)/(slope+1)))
 # Do not combine stiffness from different confinement states: reference at30kPa tests.
 ref=[r for r in rows if r['nominal_confinement_kPa']==30 and r['Etan_MPa'] is not None]
 E=float(np.median([r['Etan_MPa'] for r in ref]))*1e6
 residual=[r['phi_residual_c0_deg'] for r in rows if r['phi_residual_c0_deg'] is not None]
 fits[str(den)]=dict(count=len(rows),test_ids=[r['test_id'] for r in rows],flagged_ids=[r['test_id'] for r in rows if r['stress_angle_review_flag']],slope_c0=slope,phi_deg=phi,rms_sigma1_kPa=float(np.sqrt(np.mean((slope*x-y)**2))),E_reference_test_ids=[r['test_id'] for r in ref],E_median_Pa=E,phi_residual_median_deg=float(np.median(residual)),method='all observations retained; through-origin sigma1 vs sigma3 least squares c=0; descriptive fit NOT independent validation; static tangent E median of nominal30kPa subgroup; nu .30 assumed')
 for state,angle in [('peak',phi),('residual',float(np.median(residual)))]:
  M[f'moon_jsc1a_{int(den*10)}_{state}']=dict(role='dry simulant ground reference',law='DruckerPrager',**isotropic(E,.3),density_kg_m3=float(np.median([r['density_g_cm3'] for r in rows]))*1000,phi_deg=angle,c_Pa=0.,psi_deg=0.,status='SUMMARY_STRENGTH_FIT_WITH_ASSUMED_ELASTIC_AND_FLOW',sources=['lasp_jsc1a','dp'],assumptions=['constant static E at30kPa subgroup; not Gmax','psi0 non-dilatant baseline, not measured dilation fit','no postpeak transition; peak/residual are SEPARATE bounding laws','reference test confinement15-60kPa; shallow lunar confinement often below evidence','temperature,vacuum,reduced-gravity grain response not calibrated'])
M['mars_jez1_peak']=dict(role='dry simulant ground reference',law='DruckerPrager',**isotropic(30e6,.3),density_kg_m3=1840,phi_deg=40.6,c_Pa=3100,psi_deg=0,status='PUBLISHED_SHEAR_STRENGTH_STIFFNESS_UNIDENTIFIED',sources=['jez1','dp'],assumptions=['E30MPa and nu.30 explicitly assumed; E sensitivity0.5/1/2 mandatory','triaxial-matched DP cannot exactly reproduce measured direct shear at every path','no liquefaction parameters inferred from dry direct shear'])
M['mars_jez1_residual']={**M['mars_jez1_peak'],'phi_deg':39.7,'c_Pa':3200}
# Resolve the documented 80/101 kPa ambiguity by scaling G,K using d=0.5,
# keeping gamma_max at the table's101kPa reference; record this interpretation.
factor=math.sqrt(101/80)
M['earth_sand_dr50_pdmy02']=dict(role='Earth reference sand; drained effective stress skeleton',law='PressureDependMultiYield02',density_kg_m3=1900,G_Pa=100e6*factor,K_Pa=233e6*factor,phi_deg=33.5,p_ref_Pa=101000,peak_shear_strain=.1,pressure_exponent=.5,PTAngle_deg=25.5,contrac1=.045,contrac3=.15,dilat1=.06,dilat3=.15,contrac2=5.,dilat2=3.,liquefac1=1.,liquefac2=0.,e=.7,cs1=.9,cs2=.02,cs3=.7,pa_Pa=101000,numerical_c_Pa=300.,n_surfaces=20,status='OFFICIAL_SUGGESTED_EXAMPLE_NOT_SITE_CALIBRATED',sources=['pdmy02'],assumptions=['rho is saturated total density; bbarBrick has NO pore-pressure DOF','water table/buoyancy/pore pressure NOT implicit; earth saturated model must supply initial effective stress or use u-p formulation for site use','G/K scaled80to101kPa with pressure exponent0.5; source table ambiguity explicitly interpreted','numerical_c is NOT soil Mohr-Coulomb cohesion'])
M['earth_sand_dr50_dry_reference']={**M['earth_sand_dr50_pdmy02'],'density_kg_m3':2650/1.7,'role':'Earth dry reference sand; dry total stress equals effective stress','status':'ASSUMED_DRY_DENSITY_OFFICIAL_SKELETON_PARAMETERS','assumptions':['grain density2650kg/m3 assumed, e.70 gives dry bulk1558.82kg/m3','constitutive skeleton suggested sand parameters; not a site or dry-soil experimental calibration','G/K scaled80to101kPa with exponent.5 to resolve source reference-pressure ambiguity','no groundwater/pore-pressure DOF; numerical_c300Pa NOT physical cohesion']}
for key in ['earth_sand_dr50_pdmy02','earth_sand_dr50_dry_reference']:
 K,G=M[key]['K_Pa'],M[key]['G_Pa'];M[key]['E_Pa']=9*K*G/(3*K+G);M[key]['nu']=(3*K-2*G)/(2*(3*K+G))
M['shield_loose_jsc1a']=dict(role='loose shielding mass plus gravity load on retaining shell',law='MassOnly',density_kg_m3=1700,status='BULK_DENSITY_REFERENCE_WITH_SETTLEMENT_SENSITIVITY',sources=['lasp_jsc1a'],assumptions=['nominal1.7g/cc test condition; thickness explicit; density range1700-1900 for sensitivity','no fictitious shell tensile/flexural stiffness; retaining shell is separate material','soil redistribution/contact/settlement not captured by fixed attached mass'])
interfaces={
 'jez1_smooth_steel_peak':dict(mu=math.tan(math.radians(18.9)),adhesion_Pa=1500,source='jez1',scope='measured smooth steel, normal25-100kPa'),
 'jez1_smooth_steel_residual':dict(mu=math.tan(math.radians(16.7)),adhesion_Pa=1500,source='jez1',scope='measured smooth steel, normal25-100kPa'),
 'generic_unbonded':dict(mu=None,adhesion_Pa=0,scope='requires chosen contact pair; do not set mu=tan(soil phi) automatically')}
save(ROOT/'configs'/'catalog.json',dict(schema=1,units='m N s kg Pa K radians(unless _deg)',gravity_m_s2=dict(Earth=9.81,Moon=1.62,Mars=3.711),materials=M,interfaces=interfaces,sources=sources))
save(ROOT/'research'/'moon_strength_fits.json',fits)
print('Wrote',len(M),'material records. 6061 E at293K GPa=',Eal/1e9)
