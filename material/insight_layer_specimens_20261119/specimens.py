"""Data -> OpenSees model -> single-element analysis -> plots. SI units.

Run all: python specimens.py
Run one: python specimens.py --case duricrust_dp
Register only: register_material(ops, record, tag)
"""
import argparse
import csv
import json
import math
import os
import time
import traceback
from pathlib import Path

import numpy as np
import openseespy.opensees as ops

HERE = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/capstone_specimen_mpl')
XYZ = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]


# 1. DATA: all specimen parameters and evidence are in materials.json.
def load_materials():
    return json.loads((HERE/'materials.json').read_text())['materials']


# 2. MODEL: explicit material constructors, one 1m cube or a two-node spring.
def register_material(o, record, tag=1):
    p, law = record['parameters'], record['family']
    if law == 'ElasticIsotropic':
        o.nDMaterial(law, tag, p['E_Pa'], p['nu'], p['density_kg_m3'])
    elif law == 'DruckerPrager':
        phi, psi = math.radians(p['phi_deg']), math.radians(p['psi_deg'])
        # Triaxial-compression mapping; this matches the existing verified project DP mapping.
        sy = 6*p['c_Pa']*math.cos(phi)/(3-math.sin(phi))
        r = 2*math.sqrt(2)*math.sin(phi)/(math.sqrt(3)*(3-math.sin(phi)))
        rb = 2*math.sqrt(2)*math.sin(psi)/(math.sqrt(3)*(3-math.sin(psi)))
        o.nDMaterial(law, tag, p['K_Pa'], p['G_Pa'], sy, r, rb,
                     0., 0., 0., 0., 0., 1., p['density_kg_m3'], 101000.)
    elif law == 'PressureDependMultiYield02':
        o.nDMaterial(law, tag, 3, p['density_kg_m3'], p['G_Pa'], p['K_Pa'], p['phi_deg'],
                     p['peak_shear_strain'], p['p_ref_Pa'], p['pressure_exponent'], p['PTAngle_deg'],
                     p['contrac1'], p['contrac3'], p['dilat1'], p['dilat3'], p['n_surfaces'],
                     p['contrac2'], p['dilat2'], p['liquefac1'], p['liquefac2'], p['e'],
                     p['cs1'], p['cs2'], p['cs3'], p['pa_Pa'], p['numerical_c_Pa'])
    elif law in ('QzSimple1','QzSimple2'):
        o.uniaxialMaterial(law, tag, p['soil_type'], p['ultimate_force_N'], p['z50_m'], p['suction'], p['damping'])
    elif law == 'PySimple2':
        o.uniaxialMaterial(law, tag, p['soil_type'], p['ultimate_force_N'], p['y50_m'], p['drag'], p['damping'])
    elif law == 'TzSimple2':
        o.uniaxialMaterial(law, tag, p['soil_type'], p['ultimate_force_N'], p['z50_m'], p['damping'])
    elif law == 'ElasticPPGap':
        o.uniaxialMaterial(law, tag, p['stiffness_N_m'], p['yield_force_N'], p['gap_m'], p['hardening_ratio'], p['damage'])
    elif law == 'Hysteretic':
        branches = [v for point in p['positive_backbone']+p['negative_backbone'] for v in point]
        o.uniaxialMaterial(law, tag, *branches, p['pinchX'], p['pinchY'], p['damage1'], p['damage2'], p['beta'])
    else:
        raise ValueError(f'Unsupported family: {law}')


def cube(record):
    ops.wipe()
    ops.model('basic','-ndm',3,'-ndf',3)
    register_material(ops,record)
    for n,xyz in enumerate(XYZ,1):
        ops.node(n,*xyz)


def configure(increment=.1):
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormUnbalance',1e-5,100,0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl',increment)
    ops.analysis('Static')


def checked(steps=1):
    if ops.analyze(steps):
        raise RuntimeError(f'OpenSees nonconvergence at pseudo-time {ops.getTime()}')


def gp_stress():
    raw=np.asarray(ops.eleResponse(1,'stress'),dtype=float)
    if raw.size in (6,7):
        if not np.isfinite(raw).all():
            raise RuntimeError('Invalid single-point stress response')
        return raw[:6].reshape(1,6)
    s=np.asarray([ops.eleResponse(1,'material',gp,'stress')[:6] for gp in range(1,9)],dtype=float)
    if s.shape!=(8,6) or not np.isfinite(s).all():
        raise RuntimeError('Invalid eight-point stress response')
    return s


def write_csv(path,header,rows):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)


# 3. ANALYSIS: material tests only. No dome, gravity, field profile or global ground run.
def verify_mass(record):
    cube(record)
    for n in range(1,9):
        ops.fix(n,0,1,1)
    for n in range(2,9):
        ops.equalDOF(1,n,1)
    ops.element('bbarBrick',1,*range(1,9),1,0.,0.,0.)
    ops.node(9,0.,0.,0.)
    ops.fix(9,1,1,1)
    ops.uniaxialMaterial('Elastic',2,1e6)
    ops.element('zeroLength',2,9,1,'-mat',2,'-dir',1)
    configure(0.)
    mass=1e6/ops.eigen('-fullGenLapack',1)[0]
    expected=record['parameters']['density_kg_m3']
    if not math.isclose(mass,expected,rel_tol=1e-7):
        raise RuntimeError(f'Incorrect cube mass {mass}; expected {expected}')
    return dict(expected_kg=expected,recovered_kg=mass,relative_error=abs(mass/expected-1),passed=True)


def triaxial(record,pressure,out):
    cube(record)
    for n,(x,y,z) in enumerate(XYZ,1):
        ops.fix(n,int(x==0),int(y==0),int(z==0))
    ops.element('bbarBrick',1,*range(1,9),1,0.,0.,0.)
    ops.timeSeries('Linear',1)
    ops.pattern('Plain',1,1)
    for n,(x,y,z) in enumerate(XYZ,1):
        ops.load(n,-pressure*x/4,-pressure*y/4,-pressure*z/4)
    configure(.1)
    checked(10)
    ops.loadConst('-time',0.)
    u0=ops.nodeDisp(5,3)
    p=record['parameters']
    plastic=record['family']=='DruckerPrager'
    if plastic:
        sine=math.sin(math.radians(p['phi_deg']))
        q_expected=(2*pressure*sine+2*p['c_Pa']*math.sqrt(1-sine*sine))/(1-sine)
        ey=q_expected/p['E_Pa']
        # Resolve the elastic branch even for stiff blocky ejecta at low pressure.
        strains=np.r_[np.linspace(0,.8*ey,41),np.linspace(.8*ey,3*ey,111)[1:],np.linspace(3*ey,2*ey,61)[1:]]
    else:
        strains=np.r_[np.linspace(0,1e-4,81),np.linspace(1e-4,0,81)[1:]]
    ops.timeSeries('Path',2,'-dt',1.,'-values',*(u0-strains).tolist(),'-useLast')
    ops.pattern('Plain',2,2)
    for n in range(5,9):
        ops.sp(n,3,1.)
    ops.wipeAnalysis()
    configure(0.)
    checked()
    ops.integrator('LoadControl',1.)
    rows=[];max_lat=0.;max_uniform=0.;max_balance=0.
    for i in range(len(strains)):
        if i:
            checked()
        s=gp_stress();mean=s.mean(axis=0)
        q=-mean[2]+(mean[0]+mean[1])/2
        max_lat=max(max_lat,float(np.max(abs(s[:,:2]+pressure))))
        max_uniform=max(max_uniform,float(np.max(abs(s-mean))))
        ops.reactions()
        # Unit area: net top SP reaction equals minus deviatoric axial force.
        max_balance=max(max_balance,abs(sum(ops.nodeReaction(n,3) for n in range(5,9))+q))
        rows.append([i,u0-ops.nodeDisp(5,3),q,*mean])
    a=np.asarray(rows)
    E_actual=float(np.polyfit(a[:20,1],a[:20,2],1)[0])
    scale=max(float(np.max(abs(a[:,2]))),pressure,1.)
    checks=dict(elastic_modulus=abs(E_actual/p['E_Pa']-1)<1e-5,
                confinement=max_lat/max(pressure,1)<1e-5,homogeneity=max_uniform/scale<1e-6,
                force_balance=max_balance/scale<1e-5,finite=bool(np.isfinite(a).all()))
    metrics=dict(pressure_Pa=pressure,initial_E_Pa=E_actual,expected_E_Pa=p['E_Pa'],
                 confinement_max_error_Pa=max_lat,force_balance_max_error_N=max_balance,
                 final_strain=float(a[-1,1]),final_q_Pa=float(a[-1,2]))
    if plastic:
        plateau=float(a[150,2])
        metrics.update(expected_yield_q_Pa=q_expected,computed_yield_q_Pa=plateau,
                       yield_relative_error=abs(plateau/q_expected-1),expected_residual_strain=2*ey)
        checks.update(yield_strength=abs(plateau/q_expected-1)<1e-4,elastic_unloading=abs(a[-1,2])/q_expected<1e-4)
    else:
        checks.update(elastic_curve=bool(np.allclose(a[:,2],p['E_Pa']*a[:,1],rtol=1e-6,atol=.01)),return_to_origin=abs(a[-1,2])<.01)
    checks={key:bool(value) for key,value in checks.items()}
    write_csv(out/f'triaxial_p{int(pressure)}.csv',['step','axial_compression_strain','q_Pa','sx_Pa','sy_Pa','sz_Pa','sxy_Pa','syz_Pa','szx_Pa'],rows)
    return dict(test='triaxial_load_unload',**metrics,checks=checks,passed=all(checks.values()))


def cyclic_shear(record,pressure,out):
    cube(record)
    for n in range(1,9):
        ops.fix(n,*((1,1,1) if n<=4 else (0,1,0)))
    for n in (6,7,8):
        ops.equalDOF(5,n,1,3)
    ops.element('SSPbrick',1,*range(1,9),1,0.,0.,0.)
    ops.updateMaterialStage('-material',1,'-stage',0)
    ops.timeSeries('Linear',1)
    ops.pattern('Plain',1,1)
    ops.load(5,0.,0.,-pressure)
    configure(.1)
    # Multi-yield material has a piecewise tangent. Check displacement convergence
    # then independently enforce force/stress balance on the saved response.
    ops.test('NormDispIncr',1e-11,100,0)
    checked(10)
    ops.loadConst('-time',0.)
    ops.updateMaterialStage('-material',1,'-stage',1)
    ops.integrator('LoadControl',0.)
    checked()
    mean_p0=float(-np.mean(gp_stress()[:,:3]))
    ops.timeSeries('Linear',2)
    ops.pattern('Plain',2,2)
    ops.load(5,1.,0.,0.)
    targets=[1e-5,-1e-5,.0001,-.0001,.001,-.001,.005,-.005,0.]
    rows=[];start=0.;refinements=0
    def advance(target,depth=0):
        nonlocal refinements
        previous=ops.nodeDisp(5,1)
        ops.integrator('DisplacementControl',5,1,float(target-previous))
        if ops.analyze(1)==0:
            return
        if depth>=8:
            raise RuntimeError(f'cyclic shear failed after 8 step bisections at target gamma={target:g}')
        refinements+=1
        advance((previous+target)/2,depth+1)
        advance(target,depth+1)
    failure=None
    for target in targets:
        for u in np.linspace(start,target,61)[1:]:
            try:
                advance(u)
            except RuntimeError as e:
                failure=str(e)
                break
            s=gp_stress().mean(axis=0)
            ops.reactions()
            applied=ops.getLoadFactor(2)
            reaction=sum(ops.nodeReaction(n,1) for n in range(1,5))
            rows.append([u,s[5],s[2],applied,reaction,ops.nodeDisp(5,3),-float(np.mean(s[:3]))])
        if failure:
            break
        start=target
    if not rows:
        raise RuntimeError(failure or 'No shear steps completed')
    a=np.asarray(rows);scale=max(float(np.max(abs(a[:,3]))),1.)
    eq=float(np.max(abs(a[:,3]+a[:,4]))/scale)
    vert=float(np.max(abs(a[:,2]+pressure))/pressure)
    stress_force=float(np.max(abs(a[:,1]-a[:,3]))/scale)
    checks=dict(complete=failure is None,finite=bool(np.isfinite(a).all()),force_balance=eq<1e-5,vertical_stress=vert<1e-5,stress_force=stress_force<1e-5)
    write_csv(out/f'cyclic_shear_p{int(pressure)}.csv',['gamma_xz','tau_xz_Pa','sigma_z_Pa','applied_force_N','base_reaction_N','vertical_displacement_m','mean_effective_pressure_Pa'],rows)
    return dict(test='drained_constant_vertical_stress_cyclic_shear',pressure_Pa=pressure,
                element='SSPbrick',initial_mean_pressure_Pa=mean_p0,steps=len(rows),step_bisections=refinements,error=failure,peak_tau_Pa=float(np.max(abs(a[:,1]))),
                final_tau_Pa=float(a[-1,1]),force_balance_relative_error=eq,vertical_stress_relative_error=vert,
                scope='No pore fluid. Lateral strain constrained; mean pressure evolves and is not equal to vertical pressure.',checks=checks,passed=all(checks.values()))


def spring_test(record,out):
    ops.wipe()
    ops.model('basic','-ndm',3,'-ndf',3)
    register_material(ops,record)
    direction=record['geometry']['direction']
    for n in (1,2):
        ops.node(n,0.,0.,0.)
        ops.fix(n,*((1,1,1) if n==1 else tuple(int(d!=direction) for d in (1,2,3))))
    ops.element('zeroLength',1,1,2,'-mat',1,'-dir',direction)
    # The official material tester commits each setStrain. The response is force N
    # versus deformation m for this spring. It is not a brick stress response.
    ops.testUniaxialMaterial(1)
    scale=record['displacement_scale_m']
    targets=np.array([0.,-.1,-1.,-4.,0.,1.,4.,0.,-8.,0.,8.,0.])*scale
    rows=[]
    for start,end in zip(targets[:-1],targets[1:]):
        for u in np.linspace(start,end,101)[1:]:
            ops.setStrain(float(u))
            rows.append([u,ops.getStress(),ops.getTangent()])
    a=np.asarray(rows);p=record['parameters'];law=record['family']
    checks=dict(finite=bool(np.isfinite(a).all()),responds=bool(np.max(abs(a[:,1]))>0))
    if 'ultimate_force_N' in p:
        checks['force_capacity_bound']=bool(np.max(abs(a[:,1]))<=p['ultimate_force_N']*1.001)
    if law=='ElasticPPGap':
        checks.update(compression_only=bool(np.max(a[:,1])<1e-6),yield_capacity=math.isclose(float(np.min(a[:,1])),p['yield_force_N'],rel_tol=1e-6))
    if law=='Hysteretic':
        # First virgin compression point remains in initial elastic range.
        expected=p['negative_backbone'][0][0]/p['negative_backbone'][0][1]
        checks['initial_stiffness']=math.isclose(a[0,1]/a[0,0],expected,rel_tol=1e-5)
    checks={key:bool(value) for key,value in checks.items()}
    write_csv(out/'cyclic_spring.csv',['relative_displacement_m','force_N','tangent_N_m'],rows)
    return dict(test='committed_uniaxial_material_history',steps=len(rows),max_force_N=float(np.max(a[:,1])),min_force_N=float(np.min(a[:,1])),
                scope='Material tester on the same defined zeroLength law; no global spring equilibrium solve or ring area conversion.',checks=checks,passed=all(checks.values()))


def run(record,output):
    out=output/record['id'];out.mkdir(parents=True,exist_ok=True)
    (out/'specimen.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    result=dict(id=record['id'],family=record['family'],status=record['status'],mars_calibrated=False,tests=[])
    start=time.monotonic()
    try:
        if record['geometry']['type']=='cube':
            result['mass_test']=verify_mass(record)
            if record['family']=='PressureDependMultiYield02':
                for p in (5000.,80000.):
                    result['tests'].append(cyclic_shear(record,p,out))
            else:
                for p in (5000.,50000.):
                    result['tests'].append(triaxial(record,p,out))
        else:
            result['tests'].append(spring_test(record,out))
        result['passed']=all(t['passed'] for t in result['tests'])
    except Exception as e:
        result.update(passed=False,error=str(e),traceback=traceback.format_exc())
    finally:
        ops.wipe()
    result['elapsed_seconds']=time.monotonic()-start
    (out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    return result


# 4. VISUALIZATION: independent plots retain the pressure and evidence labels.
def plot_case(record,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=output/record['id']
    fig,ax=plt.subplots(figsize=(7.5,5.2))
    fig.subplots_adjust(left=.13,right=.97,top=.90,bottom=.20)
    if record['geometry']['type']=='zeroLength':
        p=out/'cyclic_spring.csv'
        if p.exists():
            a=np.loadtxt(p,delimiter=',',skiprows=1);ax.plot(a[:,0]*1000,a[:,1]/1000,lw=1)
        ax.set(xlabel='Relative displacement (mm)',ylabel='Force (kN)')
    else:
        for p in sorted(out.glob('*.csv')):
            a=np.loadtxt(p,delimiter=',',skiprows=1)
            pressure=int(p.stem.split('_p')[-1])/1000
            if p.name.startswith('triaxial'):
                ax.plot(a[:,1]*100,a[:,2]/1000,label=f'confinement {pressure:g} kPa')
                ax.set(xlabel='Axial compression strain (%)',ylabel='Deviatoric stress q (kPa)')
            else:
                ax.plot(a[:,0]*100,a[:,1]/1000,label=f'vertical stress {pressure:g} kPa',lw=1)
                ax.set(xlabel='Shear strain (%)',ylabel='Shear stress (kPa)')
        if ax.lines:
            ax.legend(fontsize=8)
    ax.set_title(record['id'],fontsize=11)
    ax.axhline(0,color='#999',lw=.5);ax.grid(alpha=.2)
    fig.text(.5,.035,'Numerical specimen check; Mars field/cyclic calibration is incomplete.',ha='center',fontsize=8,color='#8a3d16')
    fig.savefig(out/'response.png',dpi=140)
    plt.close(fig)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',help='specimen id; omit to run every candidate')
    parser.add_argument('--output',type=Path,default=HERE/'results')
    args=parser.parse_args(argv)
    records=load_materials()
    if args.case:
        records=[r for r in records if r['id']==args.case]
        if not records:
            parser.error('Unknown specimen id')
    args.output.mkdir(parents=True,exist_ok=True)
    results=[]
    for record in records:
        r=run(record,args.output);results.append(r)
        plot_case(record,args.output)
        print(f"{r['id']}: {'PASS' if r['passed'] else 'FAIL'} ({r['elapsed_seconds']:.2f}s)",flush=True)
        if not r['passed']:
            print(r.get('error',r['tests']),flush=True)
    summary=dict(opensees_version=ops.version(),scope='Numerical implementation checks, not experimental validation',
                 case_count=len(results),passed_count=sum(r['passed'] for r in results),all_pass=all(r['passed'] for r in results),cases=results)
    name='summary.json' if not args.case else f'summary_{args.case}.json'
    (args.output/name).write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    if not summary['all_pass']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
