"""Conditional acoustic-tensor screen for the perfectly-plastic DP records.

The continuum loading tangent is checked against a native single-element
tangent. A negative acoustic determinant for a sampled normal demonstrates a
loss-of-ellipticity possibility. Snapshot screening does NOT prove the cause
of any failed global run or the stability of all unsampled directions/states.
"""
import csv
import json
import math
from pathlib import Path

import numpy as np
import openseespy.opensees as ops
from run_specimens import records,base
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'stability'
PAIRS=[(0,0),(1,1),(2,2),(0,1),(1,2),(0,2)]


def stiffness(E,nu,phi,psi,stress):
    G=E/(2*(1+nu));K=E/(3*(1-2*nu));I=np.eye(3)
    C=(K-2*G/3)*np.einsum('ij,kl->ijkl',I,I)+G*(np.einsum('ik,jl->ijkl',I,I)+np.einsum('il,jk->ijkl',I,I))
    s=stress-np.trace(stress)/3*I;D=s/np.linalg.norm(s)
    def rho(x):
        sine=math.sin(math.radians(x));return 2*math.sqrt(2)*sine/(math.sqrt(3)*(3-sine))
    nf=D+rho(phi)*I;ng=D+rho(psi)*I
    A=np.einsum('ijkl,kl->ij',C,ng);B=np.einsum('ij,ijkl->kl',nf,C)
    Cp=C-np.einsum('ij,kl->ijkl',A,B)/np.einsum('ij,ij',nf,A)
    return C,Cp


def native_check(record):
    #3.4.0 DP's generic factory needs a ThreeDimensional copy, as used by
    # bbarBrick. The identity wrapper requests this copy. InitialStateAnalysis
    # is never switched on/off, so no initial-strain reset occurs.
    base.cube(record);ops.nDMaterial('InitialStateAnalysisWrapper',2,1,3)
    p=record['parameters'];pressure=5000.;phi=math.radians(p['phi_deg'])
    tau=2*math.sqrt(3)*math.cos(phi)/(3-math.sin(phi))*(p['c_Pa']+pressure*math.tan(phi))
    gamma_y=tau/p['G_Pa'];eps=np.array([-pressure/(3*p['K_Pa'])]*3+[0.,0.,0.])
    def response(strain):
        ops.NDTest('SetStrain',2,*strain.tolist())
        return np.asarray(ops.NDTest('GetStress',2))
    initial=response(eps);assert np.allclose(initial[:3],-pressure,rtol=1e-10)
    ops.NDTest('CommitState',2)
    for gamma in np.linspace(0,3*gamma_y,61)[1:]:
        eps[5]=gamma;response(eps);ops.NDTest('CommitState',2)
    # All central perturbations remain on the plastic-loading side of the
    # committed state. Central differences at the yield corner itself would
    # incorrectly average loading and elastic-unloading tangents.
    delta=gamma_y*1e-5;h=delta/1000;trial=eps.copy();trial[5]+=delta
    stress=response(trial);native=np.empty((6,6))
    for j in range(6):
        plus=trial.copy();minus=trial.copy();plus[j]+=h;minus[j]-=h
        native[:,j]=(response(plus)-response(minus))/(2*h)
    S=np.array([[stress[0],stress[3],stress[5]],[stress[3],stress[1],stress[4]],[stress[5],stress[4],stress[2]]])
    Ce,Cp=stiffness(p['E_Pa'],p['nu'],p['phi_deg'],p['psi_deg'],S)
    expected=np.array([[Cp[i,j,k,l] for k,l in PAIRS] for i,j in PAIRS])
    err=float(np.max(abs(native-expected))/p['G_Pa'])
    ops.wipe()
    return dict(specimen=record['id'],method='Finite difference of native NDTest trial stresses; no native tangent recorder is exposed for this3.4.0 material',native_versus_continuum_tangent_max_error_over_G=err,limit_relative=1e-3,passed=err<1e-3,engineering_strain_order='xx yy zz xy yz xz; shear strain is engineering gamma',small_plastic_trial_increment=delta,central_difference_increment=h)


def normals():
    values=[];labels=[]
    for a,b in [(0,1),(0,2),(1,2)]:
        for angle in np.linspace(0,math.pi/2,181):
            n=np.zeros(3);n[a]=math.cos(angle);n[b]=math.sin(angle);values.append(n);labels.append((a,b,math.degrees(angle)))
    return np.array(values),labels


def determinants(principal_deviator,nu,phi,psi=0.):
    # Work in the stress principal basis. Q/G = I + a n.n - An.Bn/den.
    # Matrix determinant lemma gives the determinant without large tensor arrays.
    N,_=normals();k=2*(1+nu)/(3*(1-2*nu));a=k+1/3
    def rho(x):
        sine=math.sin(math.radians(x));return 2*math.sqrt(2)*sine/(math.sqrt(3)*(3-sine))
    rf=rho(phi);rg=rho(psi);den=2+9*k*rf*rg
    V=principal_deviator[:,None,:]*N[None,:,:]
    A=2*V+3*k*rg*N;B=2*V+3*k*rf*N
    dot=np.sum(A*B,axis=2)-a/(1+a)*np.sum(A*N,axis=2)*np.sum(B*N,axis=2)
    return (1+a)*(1-dot/den)


def main():
    OUT.mkdir(exist_ok=True)
    selected=[next(r for r in records() if r['id']==name) for name in ['mgs1_ss','blocky_phi35']]
    native=[native_check(r) for r in selected]
    if not all(r['passed'] for r in native):raise RuntimeError('Continuum tangent derivation must match native tangent before interpreting stability')
    # Independent algebra/control: associated flow has no negative sampled
    # determinant in the two reference stress states; nonassociated shear can.
    D=np.array([[1.,1.,-2.],[1.,0.,-1.]])
    D=D/np.linalg.norm(D,axis=1)[:,None]
    control=[]
    for r in selected:
        phi=r['parameters']['phi_deg'];nu=r['parameters']['nu']
        for psi in [0.,phi]:
            det=determinants(D,nu,phi,psi)
            control.append(dict(specimen=r['id'],nu=nu,phi_deg=phi,psi_deg=psi,triaxial_min_det_Q_over_G3=float(det[0].min()),pure_shear_min_det_Q_over_G3=float(det[1].min())))
    assert all(r['pure_shear_min_det_Q_over_G3']>=-1e-10 for r in control if r['psi_deg']==r['phi_deg'])
    # Check the optimized determinant expression against direct 4th-order
    # contraction in the same stress principal basis.
    N,_=normals();proof=[]
    for r in selected:
        p=r['parameters'];_,Cp=stiffness(p['E_Pa'],p['nu'],p['phi_deg'],0.,np.diag(D[1]))
        direct=np.array([np.linalg.det(np.einsum('j,ijkl,l->ik',n,Cp,n)/p['G_Pa']) for n in N])
        efficient=determinants(D[1:2],p['nu'],p['phi_deg'])[0]
        error=float(np.max(abs(direct-efficient)));assert error<1e-10
        proof.append(dict(specimen=r['id'],direct_tensor_vs_matrix_lemma_max_absolute_error=error,passed=True))
    path_screen=[]
    for r in records():
        p=r['parameters'];det=determinants(D,p['nu'],p['phi_deg'],p['psi_deg'])
        path_screen.append(dict(specimen=r['id'],triaxial_min_det_Q_over_G3=float(det[0].min()),pure_shear_min_det_Q_over_G3=float(det[1].min()),scope='Conditional perfect-plastic loading tangent; two reference stress directions, not all field states'))
    with (OUT/'candidate_path_stability.csv').open('w',newline='',encoding='utf-8-sig') as f:w=csv.DictWriter(f,fieldnames=list(path_screen[0]));w.writeheader();w.writerows(path_screen)
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for ax,r in zip(axes,selected):
        p=r['parameters'];angle=np.linspace(0,90,181)
        for psi,label in [(0.,'Actual model: dilation0'),(p['phi_deg'],'Associated-flow control')]:
            # The0–2 principal plane is the second181-point block.
            values=determinants(D[1:2],p['nu'],p['phi_deg'],psi)[0,181:362]
            ax.plot(angle,values,label=label)
        ax.axhline(0,color='black',lw=.8);ax.set(title=r['id']+' | pure-shear loading tangent',xlabel='Normal angle in principal-stress plane (deg)',ylabel='det(Q) / G^3');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(OUT/'acoustic_determinant.png',dpi=150);plt.close(fig)
    rows=[]
    run=ROOT/'ground/results/nonlinear_layers_20261125/refined_gravity100/nl_ss_block35_r2/run'
    summary=json.loads((run/'summary.json').read_text())
    for stage in ['gravity','depressurized','pressurized']:
        a=np.genfromtxt(run/(stage+'_material_state.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
        for name in np.unique(a['material']):
            d=a[a['material']==name];m=summary['ground_materials'][name];p=d['mean_compression_Pa'];q=d['q_Pa'];phi=math.radians(m['phi_deg']);sine=math.sin(phi);rho=2*math.sqrt(2)*sine/(math.sqrt(3)*(3-sine));sy=6*m['c_Pa']*math.cos(phi)/(3-sine)
            f=np.sqrt(2/3)*(q-sy)-3*rho*p
            mask=(abs(f)<np.maximum(.01,abs(q)*1e-6))&(q>1.)
            near=d[mask];min_det=None;negative=0
            if len(near):
                S=np.zeros((len(near),3,3));S[:,0,0]=near['sx_Pa'];S[:,1,1]=near['sy_Pa'];S[:,2,2]=near['sz_Pa'];S[:,0,1]=S[:,1,0]=near['sxy_Pa'];S[:,1,2]=S[:,2,1]=near['syz_Pa'];S[:,0,2]=S[:,2,0]=near['szx_Pa']
                principal=np.linalg.eigvalsh(S);principal-=principal.mean(axis=1)[:,None];principal/=np.linalg.norm(principal,axis=1)[:,None]
                mins=[]
                for chunk in np.array_split(principal,max(1,math.ceil(len(principal)/200))):
                    det=determinants(chunk,m['nu'],m['phi_deg'],m['psi_deg']);mins.extend(np.min(det,axis=1).tolist())
                min_det=min(mins);negative=sum(x < -1e-6 for x in mins)
            rows.append(dict(stage=stage,material=name,total_GP=len(d),near_yield_GP=len(near),sampled_normals=543,negative_sampled_det_GP=negative,min_det_Q_over_G3=min_det,scope='Prospective continuum plastic-loading tangent at snapshot stress; near-yield status does not prove current plastic loading; negative determinant is a local diagnostic, not global collapse'))
    result=dict(native_checks=native,determinant_formula_checks=proof,analytical_controls=control,candidate_path_screen=path_screen,ground_snapshot_screen=rows,model='3D perfect-plastic DP, nonassociated psi0, continuum loading tangent',source='R31 Sabet/deBorst2019 for general ellipticity/mesh-dependence rationale; DP derivation and native tangent check are project calculations',limitations=['Not a measured material instability experiment','543 normals sampled in three principal planes, including duplicate axis endpoints; absence of negative determinants does not prove stability','Stored phase snapshots, not all intermediate committed stress states','No Cosserat/nonlocal/viscous regularisation implemented','Does not establish the cause of an unreturned native call or the global soil bearing capacity'])
    (OUT/'stability_assessment.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    with (OUT/'ground_snapshot_screen.csv').open('w',newline='',encoding='utf-8-sig') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
