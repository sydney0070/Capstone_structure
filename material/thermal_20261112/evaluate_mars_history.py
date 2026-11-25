"""Run the original TWINS proxy history through all19 material-point laws.

Quasistatic prescribed total-strain boundary, not a coupled SSI or component
heat-transfer analysis. No air-to-material temperature equality is inferred.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from thermal import catalog,thermal_strain,mass_only_state
from verify_thermal import native,MECH

def run(args):
    if not args.accept_research_assumptions:raise ValueError('Explicit --accept-research-assumptions required')
    if not np.isfinite(args.assumed_cte_scale) or args.assumed_cte_scale<=0:raise ValueError('positive finite assumed CTE scale required')
    d=catalog();out=args.output_dir;out.mkdir(parents=True,exist_ok=False)
    ref=d['environment']['reference_temperature_K'] if args.reference_temperature_K is None else args.reference_temperature_K
    source=Path(d['environment'].get('profile_csv',Path(__file__).parent/'mars_temperature_proxy.csv'))
    with source.open() as h:history=list(csv.DictReader(h))
    time=np.array([float(r['time_s']) for r in history]);temperature=np.array([float(r['temperature_K']) for r in history])
    summary=dict(boundary=args.boundary,reference_temperature_K=ref,assumed_cte_scale=args.assumed_cte_scale,rows=len(history),
        temperature_assignment_status=d['environment'].get('temperature_assignment_status'),temperature_history_path=str(source),
        study='common existing structure-environment proxy temperature applied to material-point benchmarks, NOT actual component/soil temperature; unattenuated envelope for comparison',
        quasistatic=True,rate_or_creep_effects=False,full_SSI_coupling=False,experimental_thermal_calibration=False,mechanical_strength_temperature_dependence=False,records=[])
    for name,r in d['materials'].items():
        scale=args.assumed_cte_scale if r['status'].startswith('PROVISIONAL') else 1.
        eth=np.array([thermal_strain(name,float(T),ref,allow_research_assumptions=True)*scale for T in temperature])
        if MECH[name]['law']=='MassOnly':
            # Mass conservation is retained. No false thermal force/stiffness
            # is introduced merely because the metadata contains a coefficient.
            stress=np.zeros_like(eth);total=eth if args.boundary=='free' else np.zeros_like(eth)
            strains=eth;temps=temperature;times=time
        else:
            pre=np.zeros(6)
            if MECH[name]['law'] in ('DruckerPrager','PressureDependMultiYield02'):pre[:3]=-50000/(3*MECH[name]['K_Pa'])
            # Exactly20 preload steps then history, for native PDMY stage1.
            ramp=np.linspace(0,1,21)[:,None]*pre
            total=np.vstack([ramp,pre+(eth if args.boundary=='free' else np.zeros_like(eth))])
            strains=np.vstack([np.zeros_like(ramp),eth]);mech=total-strains
            stress=native(name,mech,soil_plastic=True)
            temps=np.r_[np.full(21,ref),temperature];times=np.r_[np.full(21,np.nan),time]
        if not np.isfinite(stress).all():raise RuntimeError('Nonfinite native material stress')
        with (out/(name+'.csv')).open('w',newline='') as h:
            w=csv.writer(h);w.writerow(['stage','time_s','temperature_K']+['total_'+k for k in ('xx','yy','zz','xy','yz','xz')]+['thermal_'+k for k in ('xx','yy','zz','xy','yz','xz')]+['stress_'+k+'_Pa' for k in ('xx','yy','zz','xy','yz','xz')])
            for t,T,e,et,s in zip(times,temps,total,strains,stress):w.writerow(['preload' if np.isnan(t) else 'thermal_history','' if np.isnan(t) else t,T,*e,*et,*s])
        summary['records'].append(dict(name=name,status=r['status'],applied_assumed_CTE_scale=scale,finite=True,
            thermal_stress_applicable=r['mechanical_thermal_stress_applicable'],maximum_stress_vector_norm_Pa=float(np.max(np.linalg.norm(stress,axis=1))),
            maximum_thermal_strain_components=np.max(abs(eth),axis=0).tolist(),
            mechanical_stage='50kPa isotropic preload then PDMYstage1' if MECH[name]['law']=='PressureDependMultiYield02' else 'original mechanical law'))
        print(name,flush=True)
    summary['completed_records']=len(summary['records']);summary['complete']=len(summary['records'])==19
    summary['source_hashes']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in Path(__file__).parent.iterdir() if f.suffix in ('.py','.json','.csv')}
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    return summary
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--boundary',choices=['free','restrained'],required=True)
    p.add_argument('--reference-temperature-K',dest='reference_temperature_K',type=float);p.add_argument('--assumed-cte-scale',type=float,default=1.)
    p.add_argument('--accept-research-assumptions',action='store_true');r=run(p.parse_args());print('completed',r['completed_records'])
