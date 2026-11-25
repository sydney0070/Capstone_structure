"""Native material-point re-evaluation; no full-SSI thermal claim.

Each 3D point is a native SSPbrick driven by affine MECHANICAL strains.
Temperatures and TOTAL strains are retained separately in the raw CSVs.
All source mechanical laws/parameters remain those of the verified release.
"""
import argparse,csv,hashlib,json,sys,traceback
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
from thermal import catalog,alpha,thermal_strain,mechanical_strain,mass_only_state,laminate_thermal_resultants,BASE
from capmat.materials import register_uniaxial,register_nd

XYZ=np.array([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)],float)
MECH=json.loads((BASE/'configs/catalog.json').read_text())['materials']
UNI={'Steel02','Concrete02'}
def stiff(m):
    if m['law']=='ElasticOrthotropic':
        E1,E2,G12,nu12,nu23=m['E1_Pa'],m['E2_Pa'],m['G12_Pa'],m['nu12'],m['nu23']
        S=np.diag([1/E1,1/E2,1/E2,1/G12,2*(1+nu23)/E2,1/G12])
        S[0,1]=S[1,0]=S[0,2]=S[2,0]=-nu12/E1;S[1,2]=S[2,1]=-nu23/E2
        return np.linalg.inv(S)
    E=m.get('E_Pa',m.get('params',{}).get('E'));nu=m.get('nu',.3)
    if 'K_Pa' in m:K,G=m['K_Pa'],m['G_Pa']
    else:K,G=E/(3*(1-2*nu)),E/(2*(1+nu))
    C=np.zeros((6,6));C[:3,:3]=K-2*G/3;C[np.arange(3),np.arange(3)]+=2*G
    C[3:,3:]=np.eye(3)*G;return C

def native(name,strains,*,soil_plastic=False):
    """Commit one strain history to the actual original native material."""
    e=np.asarray(strains,float);m=MECH[name];ops.wipe()
    if e.ndim!=2 or e.shape[1]!=6 or not np.isfinite(e).all() or np.max(abs(e[0]))>1e-14:
        raise ValueError('Native path must start at zero mechanical strain and contain finite6-component rows')
    if m['law'] in UNI:
        register_uniaxial(ops,name,1);ops.testUniaxialMaterial(1);stress=[]
        for v in e:ops.setStrain(float(v[0]));stress.append([ops.getStress(),0,0,0,0,0])
        return np.array(stress)
    ops.model('basic','-ndm',3,'-ndf',3)
    if m['law']=='ElasticOrthotropic':
        E1,E2,G12,nu12,nu23=m['E1_Pa'],m['E2_Pa'],m['G12_Pa'],m['nu12'],m['nu23']
        ops.nDMaterial('ElasticOrthotropic',1,E1,E2,E2,nu12,nu23,nu12*E2/E1,G12,E2/(2*(1+nu23)),G12,m['density_kg_m3'])
    else:register_nd(ops,name,1)
    u=np.empty((len(e),8,3))
    for j,(x,y,z) in enumerate(XYZ):
        u[:,j,0]=e[:,0]*x+e[:,3]*y/2+e[:,5]*z/2
        u[:,j,1]=e[:,3]*x/2+e[:,1]*y+e[:,4]*z/2
        u[:,j,2]=e[:,5]*x/2+e[:,4]*y/2+e[:,2]*z
        ops.node(j+1,x,y,z);ops.fix(j+1,*[int(np.max(abs(u[:,j,k]))<1e-16) for k in range(3)])
    ops.element('SSPbrick',1,*range(1,9),1,0.,0.,0.)
    if m['law']=='PressureDependMultiYield02':ops.updateMaterialStage('-material',1,'-stage',0)
    # One decoupled elastic DOF keeps the equation system nonempty when all
    # physical point kinematics are prescribed. It carries no external load.
    ops.node(9,100,0,0);ops.fix(9,1,1,1);ops.node(10,101,0,0);ops.fix(10,0,1,1)
    ops.uniaxialMaterial('Elastic',99,1.);ops.element('truss',2,9,10,1.,99)
    tag=10
    for j in range(8):
        for k in range(3):
            if np.max(abs(u[:,j,k]))<1e-16:continue
            ops.timeSeries('Path',tag,'-dt',1.,'-values',*u[:,j,k].tolist(),'-useLast')
            ops.pattern('Plain',tag,tag);ops.sp(j+1,k+1,1.);tag+=1
    ops.constraints('Transformation');ops.numberer('RCM');ops.system('BandGeneral')
    ops.test('NormUnbalance',1e-9,30);ops.algorithm('Newton');ops.integrator('LoadControl',1.);ops.analysis('Static')
    rows=[]
    # Initial row is zero strain and zero stress; subsequent rows are actual
    # committed native evaluations, not a scalar elastic surrogate.
    rows.append(np.zeros(6))
    for i in range(1,len(e)):
        if ops.analyze(1):raise RuntimeError(f'{name}: native point step{i} failed')
        if soil_plastic and i==20 and m['law']=='PressureDependMultiYield02':ops.updateMaterialStage('-material',1,'-stage',1)
        s=np.asarray(ops.eleResponse(1,'stress'),float)[:6]
        if s.shape!=(6,) or not np.isfinite(s).all():raise RuntimeError('Invalid native stress')
        rows.append(s)
    return np.array(rows)

def run(out):
    out.mkdir(parents=True,exist_ok=False);d=catalog();env=d['environment'];ref=env['reference_temperature_K']
    report=dict(runtime=ops.version(),temperature_assignment_status=env.get('temperature_assignment_status','PROVISIONAL'),temperature_reference_K=ref,temperature_envelope_K=[env['minimum_K'],env['maximum_K']],
        physical_thermal_calibration_completed=False,full_SSI_thermal_coupling_completed=False,
        native_3d_thermal_wrapper_supported=False,mechanical_constants_refitted=False,records=[],checks=[])
    def check(name,passed,**kw):report['checks'].append(dict(name=name,passed=bool(passed),**kw))
    T=np.r_[ref,np.linspace(ref,env['minimum_K'],21)[1:],np.linspace(env['minimum_K'],env['maximum_K'],41)[1:],np.linspace(env['maximum_K'],ref,21)[1:]]
    for name,r in d['materials'].items():
        row=dict(name=name,status=r['status'],alpha_at_reference_per_K=alpha(name,ref).tolist(),tests=[])
        try:
            eth=np.array([thermal_strain(name,float(t),ref,allow_research_assumptions=True) for t in T])
            check(name+'_Tref_zero',np.max(abs(thermal_strain(name,ref,ref,allow_research_assumptions=True)))==0)
            if MECH[name]['law']=='MassOnly':
                states=[mass_only_state(name,float(t),mass_kg=1700.,volume_ref_m3=1.,T_ref=ref,allow_research_assumptions=True) for t in T]
                ok=all(v['mass_kg']==1700 and v['thermal_force_N']==0 and not v['structural_stiffness_added'] for v in states)
                row.update(thermal_stress_test='NOT_APPLICABLE_MASS_ONLY',assumed_density_range_kg_m3=[min(v['assumed_density_kg_m3'] for v in states),max(v['assumed_density_kg_m3'] for v in states)])
                check(name+'_mass_conservation',ok);row['tests'].append('mass conservation and no fictitious thermal force')
                stress=np.zeros((len(T),6));total=eth
            else:
                # Free expansion: total=thermal, so native mechanical strain=0.
                em=np.array([mechanical_strain(name,e,float(t),ref,allow_research_assumptions=True) for e,t in zip(eth,T)])
                sf=native(name,em)
                check(name+'_free_expansion',np.max(abs(sf))<1e-5,max_stress_Pa=float(np.max(abs(sf))))
                # Small restrained heating: independent tangent response; soil
                # check remains stage0 to isolate the thermal sign/units.
                dt=.01 if ref+.01<=max(env['maximum_K'],ref) else -.01
                es=thermal_strain(name,ref+dt,ref,allow_research_assumptions=True)
                actual=native(name,np.array([np.zeros(6),-es]))[-1]
                if MECH[name]['law'] in UNI:
                    E=MECH[name]['params']['E'] if MECH[name]['law']=='Steel02' else MECH[name]['E_Pa']
                    expected=np.array([-E*es[0],0,0,0,0,0])
                else:expected=-stiff(MECH[name])@es
                error=np.linalg.norm(actual-expected)/max(np.linalg.norm(expected),1e-12)
                check(name+'_restrained_sign_units',error<1e-3,relative_error=float(error),stress_Pa=actual.tolist())
                # Re-evaluate the full observed temperature envelope as a
                # constrained material-point scenario, not a heat-transfer field.
                total=np.zeros_like(eth);mech=-eth
                if MECH[name]['law'] in ('DruckerPrager','PressureDependMultiYield02'):
                    # Establish 50kPa hydrostatic compression first. The
                    # additional constrained thermal cycle then follows the
                    # original nonlinear law at its fixed mechanical parameters.
                    K=MECH[name]['K_Pa'];pre=np.r_[np.full(3,-50000/(3*K)),[0.,0.,0.]]
                    ramp=np.linspace(0,1,21)[:,None]*pre
                    total=np.vstack([ramp,pre+np.zeros_like(eth[1:])]);temp=np.r_[np.full(21,ref),T[1:]]
                    mech=np.array([mechanical_strain(name,e,float(t),ref,allow_research_assumptions=True) for e,t in zip(total,temp)])
                    stress=native(name,mech,soil_plastic=True);tc=temp
                else:stress=native(name,mech);tc=T
                check(name+'_temperature_cycle_finite',np.isfinite(stress).all(),peak_stress_norm_Pa=float(max(np.linalg.norm(stress,axis=1))))
                # Independent thermal strain compensation: the original
                # nonzero mechanical path must replay identically at 3 temps.
                q=np.r_[np.linspace(0,-.001,11),np.linspace(-.001,.0001,12)[1:],np.linspace(.0001,0,5)[1:]]
                path=np.zeros((len(q),6));path[:,0]=q
                if MECH[name]['law'] not in UNI:path[:,1]=path[:,2]=np.minimum(q,0.)
                reference=native(name,path);diff=0.
                for temp in (env['minimum_K'],ref,env['maximum_K']):
                    th=thermal_strain(name,float(temp),ref,allow_research_assumptions=True)
                    restored=np.array([mechanical_strain(name,e+th,float(temp),ref,allow_research_assumptions=True) for e in path])
                    diff=max(diff,float(np.max(abs(native(name,restored)-reference))))
                scale=max(float(np.max(abs(reference))),1.)
                check(name+'_mechanical_path_replay',diff/scale<1e-8,max_difference_Pa=diff)
                row['tests']=['free expansion','restrained heating sign/units','constrained temperature cycle','original native mechanical-history replay at3 temperatures']
                Tsave=tc
            if MECH[name]['law']=='MassOnly':Tsave=T
            with (out/(name+'.csv')).open('w',newline='') as h:
                w=csv.writer(h);w.writerow(['step','temperature_K']+['total_'+k for k in ('xx','yy','zz','xy','yz','xz')]+['thermal_'+k for k in ('xx','yy','zz','xy','yz','xz')]+['stress_'+k+'_Pa' for k in ('xx','yy','zz','xy','yz','xz')])
                for i,(temp,e,s) in enumerate(zip(Tsave,total,stress)):
                    w.writerow([i,float(temp),*e,*thermal_strain(name,float(temp),ref,allow_research_assumptions=True),*s])
            row['passed']=all(v['passed'] for v in report['checks'] if v['name'].startswith(name+'_'))
        except Exception as ex:row.update(passed=False,error=str(ex),traceback=traceback.format_exc())
        report['records'].append(row);(out/'summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(name,row['passed'],flush=True)
    test_temperature=env['maximum_K'] if abs(env['maximum_K']-ref)>1 else env['minimum_K']
    a0=thermal_strain('cfrp_as4_8552_ud_rtd',test_temperature,ref,allow_research_assumptions=True)
    a90=thermal_strain('cfrp_as4_8552_ud_rtd',test_temperature,ref,angle_deg=90,allow_research_assumptions=True)
    check('CFRP90_axis_swap',np.allclose(a90[:2],a0[[1,0]],atol=1e-15))
    result=laminate_thermal_resultants([0,90,90,0],[.00018796]*4,[test_temperature]*4,ref)
    check('symmetric_laminate_no_thermal_bending',np.max(abs(result['M_thermal_N']))<1e-10)
    try:thermal_strain('al6061_t6_extrusion',301,ref);guard=False
    except ValueError:guard=True
    check('no_NIST_out_of_range_extrapolation',guard)
    try:thermal_strain('mars_jez1_peak',test_temperature,ref);guard=False
    except ValueError:guard=True
    check('assumptions_require_explicit_opt_in',guard)
    # Prove that the PDMY stage change acts on the element's material copy.
    for name in ('earth_sand_dr50_pdmy02','earth_sand_dr50_dry_reference'):
        pre=np.r_[np.full(3,-50000/(3*MECH[name]['K_Pa'])),[0.,0.,0.]]
        ramp=np.linspace(0,1,21)[:,None]*pre
        shear=np.tile(pre,(41,1));shear[:,3]=np.linspace(0,.02,41)
        path=np.vstack([ramp,shear[1:]])
        elastic=native(name,path,soil_plastic=False);plastic=native(name,path,soil_plastic=True)
        difference=abs(plastic[-1,3]-elastic[-1,3])/max(abs(elastic[-1,3]),1.)
        check(name+'_native_plastic_stage_active',difference>.05,relative_shear_difference=float(difference))
    check('coverage_exact19',set(MECH)=={v['name'] for v in report['records']} and len(report['records'])==19)
    report['all_numerical_checks_pass']=all(x['passed'] for x in report['checks']) and all(x['passed'] for x in report['records'])
    report['check_count']=len(report['checks']);report['source_hashes']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in Path(__file__).parent.iterdir() if f.suffix in ('.py','.json')}
    (out/'summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));ops.wipe()
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    r=run(a.output_dir);print('checks',r['check_count'],'pass',r['all_numerical_checks_pass']);sys.exit(0 if r['all_numerical_checks_pass'] else 1)
