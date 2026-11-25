"""RC Beam verification: zero lateral load residual displacement and discretization sensitivity.
N-m-Pa; provisional C30/SD400 properties; not an experimental calibration.
"""
from pathlib import Path
from dataclasses import replace, asdict
import csv, json, traceback
import numpy as np
import openseespy.opensees as ops
from rc_section import SectionGeometry, create_rc_fiber_section

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results_rc_final_validation'; OUT.mkdir(exist_ok=True)
TARGETS=(.003,-.003,.006,-.006,.012,-.012,.018,-.018,0.)


def solver():
    ops.constraints('Plain'); ops.numberer('RCM'); ops.system('BandGeneral')
    ops.test('NormDispIncr',1e-9,80); ops.algorithm('Newton')
    ops.integrator('LoadControl',1.0); ops.analysis('Static')


def save_csv(path, rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['phase','tip_displacement_m','applied_horizontal_force_N']);w.writerows(rows)


def run_case(case_name, ny=24, nz=12, nip=5, subdivisions=60, zero_force_unloading=False):
    ops.wipe(); ops.model('basic','-ndm',2,'-ndf',3)
    geom=replace(SectionGeometry(),ny=ny,nz=nz)
    create_rc_fiber_section(ops,geom=geom)
    ops.node(1,0,0);ops.node(2,0,3.0);ops.fix(1,1,1,1)
    ops.geomTransf('PDelta',1)
    ops.beamIntegration('Lobatto',1,10,nip)
    ops.element('forceBeamColumn',1,1,2,1,1)
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(2,1,0,0)
    solver(); data=[('initial',0.,0.)]; u=0.
    case_dir=OUT/case_name; case_dir.mkdir(exist_ok=True)
    try:
        for target in TARGETS:
            du=(target-u)/subdivisions
            for i in range(subdivisions):
                ops.integrator('DisplacementControl',2,1,du)
                if ops.analyze(1)!=0:raise RuntimeError(f'displacement control failed: target={target}, step={i}')
                u=float(ops.nodeDisp(2,1));F=float(ops.getLoadFactor(1))
                data.append(('cyclic',u,F))
        end_force=float(ops.getLoadFactor(1))
        end_u=float(ops.nodeDisp(2,1))
        # Key distinction: a return to u=0 DOES NOT represent zero-load unloading.
        # Advance load factor (reference force=1 N) from its current end value to zero.
        if zero_force_unloading and abs(end_force)>1e-6:
            increment=-end_force/120.0
            for i in range(120):
                ops.integrator('LoadControl',increment)
                if ops.analyze(1)!=0:raise RuntimeError(f'zero-force unloading failed at step={i}')
                data.append(('zero_force_unload',float(ops.nodeDisp(2,1)),float(ops.getLoadFactor(1))))
        ops.reactions()
        final_f=float(ops.getLoadFactor(1)); final_u=float(ops.nodeDisp(2,1))
        reaction=float(ops.nodeReaction(1,1))
        a=np.asarray([(v,f) for _,v,f in data],dtype=float)
        save_csv(case_dir/'history.csv',data)
        out={'status':'COMPLETED','geometry':asdict(geom),'lobatto_points':nip,'subdivisions':subdivisions,
             'max_abs_force_kN':float(np.max(np.abs(a[:,1]))/1000),
             'end_cyclic_force_kN':end_force/1000,'end_cyclic_displacement_mm':end_u*1000,
             'final_force_N':final_f,'final_tip_displacement_mm':final_u*1000,
             'support_reaction_x_N':reaction,
             'equilibrium_relative_error':abs(reaction+final_f)/max(abs(final_f),1.),
             'zero_force_unloading':zero_force_unloading,
             'points':len(data)}
        return out,a
    except Exception:
        save_csv(case_dir/'partial_history.csv',data)
        raise


def compare(base,other):
    # cycle history only; both use same TARGETS and subdivisions, so aligned index.
    a=base[:541,:];b=other[:541,:]
    if len(a)!=len(b): raise ValueError('Mismatched history lengths')
    d=a[:,1]-b[:,1]
    return {'peak_force_relative_difference':float(abs(np.max(abs(a[:,1]))-np.max(abs(b[:,1])))/max(np.max(abs(a[:,1])),1)),
            'rms_force_difference_kN':float(np.sqrt(np.mean(d*d))/1000),
            'max_path_force_difference_kN':float(np.max(abs(d))/1000)}


def main():
    result={'units':'N,m,Pa','status':'INCOMPLETE','tests':{},
            'warning':'This is numerical verification of a provisional RC beam; no experimental validation.'}
    series={}
    cases=[('baseline_24x12_ip5',24,12,5,True),
           ('refined_fibers_48x24_ip5',48,24,5,False),
           ('more_integration_24x12_ip7',24,12,7,False),
           ('refined_both_48x24_ip7',48,24,7,False)]
    for name,ny,nz,nip,unload in cases:
        try:
            item,a=run_case(name,ny,nz,nip,60,unload)
            result['tests'][name]=item;series[name]=a
        except Exception as e:
            result['tests'][name]={'status':'FAILED','error':str(e),'traceback':traceback.format_exc()}
        (OUT/'final_summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
        print(name,result['tests'][name].get('status'),result['tests'][name].get('error',''))
    key='baseline_24x12_ip5'
    if key in series:
        result['comparison']={name:compare(series[key],a) for name,a in series.items() if name!=key}
    base=result['tests'].get(key,{})
    # Explicit robust evaluation independent of bool identity semantics.
    result['residual_displacement_pass']=bool(base.get('status')=='COMPLETED' and abs(base.get('final_force_N',1e99))<1e-3 and base.get('equilibrium_relative_error',1)<1e-4)
    result['status']='COMPLETED_NUMERICALLY' if all(x.get('status')=='COMPLETED' for x in result['tests'].values()) else 'PARTIAL_OR_FAILED'
    (OUT/'final_summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({'status':result['status'],'residual_displacement_pass':result['residual_displacement_pass'],'comparison':result.get('comparison',{})},ensure_ascii=False,indent=2))
    print('Saved:',OUT)

if __name__=='__main__':main()
