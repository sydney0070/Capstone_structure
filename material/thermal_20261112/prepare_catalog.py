"""Versioned thermal overlay; preserves the mechanically verified release."""
import argparse, csv, hashlib, json, math
from pathlib import Path

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]
BASE=HERE.parent/'material_release_20261101/configs/catalog.json'
SOURCES={
 'steel':{'url':'https://eurocodes.jrc.ec.europa.eu/sites/default/files/2021-12/handbook3.pdf','location':'Chapter IV, IV-7; conventional structural steel 12e-6/K, concrete 10e-6/K','scope':'ordinary building thermal-action reference; not measured Martian cryogenic product data'},
 'al_nist':{'url':'https://trc.nist.gov/cryogenics/materials/6061%20Aluminum/6061_T6Aluminum_rev.htm','location':'Linear expansion fit, 4–300 K, relative to L293; coefficients a=-412.77,b=-.30389,c=.0087696,d=-9.9821e-6, scale1e-5','scope':'published 6061-T6 fit; no mechanical strength/temper degradation calibration'},
 'concrete':{'url':'https://www.fhwa.dot.gov/publications/research/infrastructure/pavements/pccp/thermal.cfm','location':'PCC typical8–12 microstrain/K, basalt aggregate6–8, steel11–12','scope':'terrestrial reference; specific mix/moisture/cryogenic expansion not identified'},
 'cfrp':{'url':'https://www.mdpi.com/2504-477x/5/2/63','location':'Table1: AS4-8552 model alphaX1e-20/K, alphaY/Z3.261e-5/K, Vf.5742','scope':'curing-process research model; transferred thermal coefficients, not NIAR RTD qualification CTE or cryogenic measurement'},
 'regolith_gap':{'url':'https://doi.org/10.1016/j.jnoncrysol.2010.04.049','scope':'glass made from melted JSC-1A differs from porous sintered and loose simulants; no numeric coefficient transferred from that glass paper'},
 'nd_wrapper':{'url':'https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/InitStrain.html','scope':'modern3D wrapper sigma=f(eps+eps0); negative thermal strain required; unavailable in installed3.4.0'},
}
def environment():
    files=sorted((PROJECT/'ground/ground_modeling/MARS_DATA').glob('twins_calib_*.csv'))
    samples={}
    for f in files:
        with f.open() as h:
            for r in csv.DictReader(h):
                values=[]
                for k in ('BMY_TIP_ROD_TEMP','BPY_TIP_ROD_TEMP'):
                    try:v=float(r[k])
                    except (ValueError,TypeError):continue
                    if math.isfinite(v) and 160<=v<=300:values.append(v)
                if values:samples[float(r['AOBT'])]=sum(values)/len(values)
    if not samples:raise RuntimeError('No valid original TWINS proxy temperatures')
    ordered=sorted(samples.items());values=[v for _,v in ordered]
    result=dict(reference_temperature_K=sum(values)/len(values),minimum_K=min(values),maximum_K=max(values),
        source_sample_count=len(values),valid_sensor_range_K=[160,300],reference_basis='arithmetic mean of valid deduplicated native-rate TWINS tip-temperature proxies; explicit research stress-free-temperature assumption',
        thermal_transfer_scope='air-proxy envelope only; NOT measured temperature of every material/soil layer; no solar/radiative/conductive heat-transfer solution',
        sources=[dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files])
    # Export a 5-minute-bin proxy history, without pretending it is synchronized
    # with the short seismic window or with temperatures inside every component.
    bins={};t0=ordered[0][0]
    for t,T in ordered:bins.setdefault(int((t-t0)//300),[]).append((t,T))
    with (HERE/'mars_temperature_proxy.csv').open('w',newline='') as h:
        w=csv.writer(h);w.writerow(['time_s','temperature_K'])
        for rows in bins.values():w.writerow([sum(t-t0 for t,T in rows)/len(rows),sum(T for t,T in rows)/len(rows)])
    return result

def main(argv=None):
    parser=argparse.ArgumentParser()
    parser.add_argument('--minimum-temperature-K',dest='minimum_K',type=float)
    parser.add_argument('--maximum-temperature-K',dest='maximum_K',type=float)
    parser.add_argument('--reference-temperature-K',dest='reference_K',type=float)
    parser.add_argument('--structure-environment',type=Path,help='JSON extracted from existing structure environment functions')
    args=parser.parse_args(argv)
    supplied=[args.minimum_K,args.maximum_K,args.reference_K]
    if any(x is not None for x in supplied) and not all(x is not None for x in supplied):parser.error('supply minimum, maximum AND reference temperatures together')
    if all(x is not None for x in supplied) and (not all(math.isfinite(x) and x>0 for x in supplied) or args.minimum_K>=args.maximum_K):parser.error('positive finite temperatures and min<max required')
    if all(x is not None for x in supplied) and (min(supplied)<4 or max(supplied)>300):parser.error('current NIST model covers4–300K; outside this range requires a new documented temperature law')
    base=json.loads(BASE.read_text());env=environment();overlay={}
    env['temperature_assignment_status']='PROVISIONAL_AWAITING_USER_NUMERIC_INPUT'
    env['observed_proxy_minimum_K']=env['minimum_K'];env['observed_proxy_maximum_K']=env['maximum_K']
    env['observed_proxy_mean_K']=env['reference_temperature_K']
    if args.structure_environment is not None:
        if any(x is not None for x in supplied):parser.error('choose existing structure environment OR explicit numeric temperatures')
        env=json.loads(args.structure_environment.read_text())
        if env.get('temperature_assignment_status')!='EXTRACTED_EXISTING_STRUCTURE_ENVIRONMENT':parser.error('expected extracted structure environment')
        env['structure_environment_source']=str(args.structure_environment.resolve())
        if min(env['minimum_K'],env['reference_temperature_K'])<4 or max(env['maximum_K'],env['reference_temperature_K'])>300:parser.error('structure temperatures outside current NIST4–300K model')
    if all(x is not None for x in supplied):
        env.update(minimum_K=args.minimum_K,maximum_K=args.maximum_K,reference_temperature_K=args.reference_K,
            temperature_assignment_status='USER_SPECIFIED_NUMERIC_TEMPERATURES',reference_basis='user-specified stress-free reference temperature')
    for name,m in base['materials'].items():
        rec=dict(reference_temperature_K=env['reference_temperature_K'],research_temperature_range_K=[min(env['minimum_K'],env['reference_temperature_K']),max(env['maximum_K'],env['reference_temperature_K'])],
            experimental_thermal_calibration=False,temperature_dependent_strength_calibrated=False,
            assumptions_acknowledged_for_research=True)
        if name.startswith('steel_') or name.startswith('rebar_'):
            alpha=12e-6;status='CONVENTIONAL_REFERENCE_NOT_THERMAL_FIT';sources=['steel','concrete']
            note='12e-6/K retained from legacy dome steel; rebar transfer is conventional, not SD400 product test'
        elif name.startswith('al6061'):
            alpha=None;status='PUBLISHED_TEMPERATURE_FUNCTION';sources=['al_nist'];note='integrate published length curve between material temperature and explicit reference; do not use room-temperature23.6e-6/K at all Mars temperatures'
            rec.update(model='nist_6061_length_curve',reference_length_temperature_K=293.,data_temperature_range_K=[4.,300.],fit_error_relative_to_source_data_percent=4.,length_fit={'coefficients':[-412.77,-.30389,.0087696,-9.9821e-6,0.],'scale':1e-5,'low_break_K':18.,'low_constant':-415.45})
        elif name.startswith('concrete_'):
            alpha=10e-6;status='CONVENTIONAL_REFERENCE_NOT_MIX_CALIBRATED';sources=['steel','concrete'];note='normal-weight concrete reference, mix/aggregate/moisture and cold-service effects unmeasured'
        elif name.startswith('cfrp_'):
            alpha=None;status='TRANSFERRED_RESEARCH_ORTHOTROPIC_MODEL';sources=['cfrp'];note='ply axes1=fiber,2=transverse,3=thickness; cross-ply/laminate expansion must be calculated, not replaced by scalar0'
            rec.update(model='constant_orthotropic',alpha_principal_per_K=[0.,32.61e-6,32.61e-6],source_alpha1_per_K=1e-20,source_fiber_volume_fraction=.5742)
        elif name.startswith('sintered_'):
            alpha=8e-6;status='PROVISIONAL_ANALOGUE_ASSUMPTION';sources=['concrete','regolith_gap'];note='upper-end terrestrial basalt-aggregate analogue chosen for sensitivity, NOT a measured sintered JSC-1A coefficient; porosity-specific CTE unidentified'
            rec['sensitivity_alpha_per_K']=[4e-6,8e-6,12e-6]
        elif name.startswith('moon_') or name.startswith('mars_') or name.startswith('shield_'):
            alpha=8e-6;status='PROVISIONAL_EFFECTIVE_BULK_ASSUMPTION';sources=['concrete','regolith_gap'];note='effective dry bulk continuum trial coefficient; grain CTE does not identify bulk thermal strain, compaction/ratcheting not calibrated; peak/residual keep same thermal input'
            rec['sensitivity_alpha_per_K']=[4e-6,8e-6,12e-6]
        elif name.startswith('earth_'):
            alpha=11e-6;status='PROVISIONAL_EFFECTIVE_SKELETON_ASSUMPTION';sources=['concrete'];note='terrestrial quartz-rich aggregate analogue; effective sand skeleton CTE unmeasured; saturated record has NO water thermal expansion or thermal pore-pressure model'
            rec['sensitivity_alpha_per_K']=[5.5e-6,11e-6,16.5e-6]
        else:raise ValueError(name)
        if alpha is not None:rec.update(model='constant_isotropic',alpha_linear_per_K=alpha)
        rec.update(status=status,sources=sources,scope=note)
        if m['law']=='MassOnly':rec.update(mechanical_thermal_stress_applicable=False,scope=note+'; coefficient is bookkeeping only; preserve mass and do not create stiffness or thermal load')
        else:rec['mechanical_thermal_stress_applicable']=True
        overlay[name]=rec
    doc=dict(schema=1,created='2026-11-12',units='m N s kg Pa K; alpha is LINEAR per K; engineering strain order xx yy zz xy yz xz',
        mechanical_base_path=str(BASE),mechanical_base_sha256=hashlib.sha256(BASE.read_bytes()).hexdigest(),
        status='RESEARCH_THERMAL_INPUTS_NOT_EXPERIMENTAL_RETUNING',environment=env,materials=overlay,sources=SOURCES)
    (HERE/'thermal_catalog.json').write_text(json.dumps(doc,ensure_ascii=False,indent=2))
    full=json.loads(BASE.read_text());full['schema']='thermal-overlay-v1';full['mechanical_base_sha256']=doc['mechanical_base_sha256']
    for name,m in full['materials'].items():m['thermal']=overlay[name]
    full['thermal_environment']=env;full['thermal_sources']=SOURCES
    (HERE/'catalog_with_thermal.json').write_text(json.dumps(full,ensure_ascii=False,indent=2))
    print(json.dumps(dict(records=len(overlay),environment=env),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
