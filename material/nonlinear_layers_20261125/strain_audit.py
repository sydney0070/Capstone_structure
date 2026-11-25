"""Recover committed bbarBrick strains from archived displacements.

Kinematics are verified against a warped native element before applying them
to the ground mesh. This adds strain-range evidence without re-running loads.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import openseespy.opensees as ops
from scipy.optimize import linear_sum_assignment

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'ground/ground_modeling'))
from ground_hybrid.config import default_config
from ground_hybrid.conforming_mesh import build_conforming_ground

SIGNS=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],dtype=float)
POINTS=SIGNS/np.sqrt(3)


def gradients(X):
    natural=np.empty((8,8,3))
    for g,p in enumerate(POINTS):
        for axis in range(3):
            other=[i for i in range(3) if i!=axis]
            natural[g,:,axis]=SIGNS[:,axis]*np.prod(1+SIGNS[:,other]*p[other],axis=1)/8
    J=np.einsum('eni,gnj->egij',X,natural);weight=np.linalg.det(J)
    if np.any(weight<=0):raise ValueError('Nonpositive reference Jacobian')
    return np.einsum('gna,egaj->egnj',natural,np.linalg.inv(J)),weight


def strains(gradient,weight,U):
    du=np.einsum('eni,egnj->egij',U,gradient)
    e=np.stack([du[:,:,0,0],du[:,:,1,1],du[:,:,2,2],du[:,:,0,1]+du[:,:,1,0],du[:,:,1,2]+du[:,:,2,1],du[:,:,0,2]+du[:,:,2,0]],axis=-1)
    volume=e[:,:,:3].sum(axis=2);average=(volume*weight).sum(axis=1)/weight.sum(axis=1)
    e[:,:,:3]+=(average[:,None]-volume)[:,:,None]/3
    return e


def native_patch():
    X=(SIGNS+1)/2;X[6]+=[.12,-.07,.15]
    U=np.array([[2,-3,1],[4,2,-2],[-2,5,3],[1,1,-4],[3,-2,5],[-1,4,2],[5,1,-3],[-3,-4,4]],dtype=float)*1e-5
    ops.wipe();ops.model('basic','-ndm',3,'-ndf',3)
    for n,xyz in enumerate(X,1):ops.node(n,*xyz)
    ops.nDMaterial('ElasticIsotropic',1,1e6,.27,1000.);ops.element('bbarBrick',1,*range(1,9),1,0.,0.,0.)
    ops.node(9,0.,0.,0.);ops.node(10,0.,0.,0.);ops.fix(9,1,1,1);ops.fix(10,0,1,1)
    ops.uniaxialMaterial('Elastic',2,1e6);ops.element('zeroLength',2,9,10,'-mat',2,'-dir',1)
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
    for n,u in enumerate(U,1):
        for axis,value in enumerate(u,1):ops.sp(n,axis,value)
    ops.constraints('Transformation');ops.numberer('RCM');ops.system('BandGeneral');ops.test('NormUnbalance',1e-8,20);ops.algorithm('Newton');ops.integrator('LoadControl',1.);ops.analysis('Static')
    if ops.analyze(1):raise RuntimeError('Warped native kinematic patch failed')
    native=np.array([ops.eleResponse(1,'material',gp,'strain') for gp in range(1,9)])
    if native.shape!=(8,6):raise RuntimeError('Native strain recorder shape')
    grad,w=gradients(X[None,:,:]);offline=strains(grad,w,U[None,:,:])[0]
    cost=np.linalg.norm(native[:,None,:]-offline[None,:,:],axis=2);rows,order=linear_sum_assignment(cost)
    error=float(np.max(abs(native-offline[order])));assert error<1e-12
    # Affine completeness in the same warped reference geometry.
    A=np.array([[1.,.2,-.3],[.4,-.5,.6],[-.2,.8,.9]])*1e-4
    affine=X@A.T;expected=np.array([A[0,0],A[1,1],A[2,2],A[0,1]+A[1,0],A[1,2]+A[2,1],A[0,2]+A[2,0]])
    affine_error=float(np.max(abs(strains(grad,w,affine[None,:,:])-expected)));assert affine_error<1e-12
    ops.wipe()
    return order,dict(native_patch_max_strain_error=error,affine_patch_max_strain_error=affine_error,native_GP_to_offline_order=order.tolist(),reference='Warped8-node cube with prescribed non-affine displacements; exact affine completeness checked independently')


def shear_invariant(e):
    T=np.zeros((*e.shape[:-1],3,3));T[:,:,0,0]=e[:,:,0];T[:,:,1,1]=e[:,:,1];T[:,:,2,2]=e[:,:,2];T[:,:,0,1]=T[:,:,1,0]=e[:,:,3]/2;T[:,:,1,2]=T[:,:,2,1]=e[:,:,4]/2;T[:,:,0,2]=T[:,:,2,0]=e[:,:,5]/2
    principal=np.linalg.eigvalsh(T)
    return principal[:,:,2]-principal[:,:,0]


def main():
    out=HERE/'strain_audit';out.mkdir(exist_ok=True);order,verification=native_patch()
    run=ROOT/'ground/results/nonlinear_layers_20261125/refined_gravity100/nl_ss_block35_r2/run'
    summary=json.loads((run/'summary.json').read_text());a=np.load(run/'stage_displacements.npz')
    config=default_config('Mars');config['layers']=summary['ground_layers'];config.update(width_m=90.,length_m=90.,depth_m=30.,max_vertical_spacing_m=4.,max_aspect_ratio=200.)
    ops.wipe();g=build_conforming_ground(config,radius=20.,width=14.,kind='ring',sectors=24,radial_subdivisions=2,excavation_depth_m=3.)
    tags=np.array(sorted(g['nodes']));assert np.array_equal(tags,a['node_tags']);assert np.allclose(np.array([g['nodes'][t] for t in tags]),a['coordinates'],atol=1e-10,rtol=0)
    index={int(t):i for i,t in enumerate(tags)};elements=sorted(g['elements']);conn=np.array([[index[t] for t in g['elements'][el]['nodes']] for el in elements]);X=a['coordinates'][conn]
    grad,weight=gradients(X);names=np.array([g['elements'][el]['material'] for el in elements]);ops.wipe()
    data={stage:strains(grad,weight,a[stage+'_m'][conn])[:,order,:] for stage in ['gravity','depressurized','pressurized']}
    rows=[]
    for stage,e in data.items():
        gamma=shear_invariant(e);increment=shear_invariant(e-data['gravity'])
        stress=np.genfromtxt(run/(stage+'_material_state.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
        # CSV order is element-major and native-GP-major; verify before pairing.
        assert np.array_equal(stress['element'].reshape(-1,8)[:,0],elements)
        assert np.array_equal(stress['gauss_point'].reshape(-1,8),np.tile(np.arange(1,9),(len(elements),1)))
        pressure=stress['mean_compression_Pa'].reshape(-1,8)
        for name in np.unique(names):
            mask=names==name;K=summary['ground_materials'][name]['K_Pa'];eps=e[mask];vol=eps[:,:,:3].sum(axis=2);expected=-pressure[mask]/K
            absolute=float(np.max(abs(vol-expected)));relative=absolute/max(float(np.max(abs(expected))),1e-12)
            assert relative<1e-8,('Bbar/psi0 volumetric consistency failed',stage,name,relative)
            rows.append(dict(stage=stage,material=name,gauss_points=int(mask.sum()*8),max_engineering_shear_strain_pct=float(gamma[mask].max()*100),max_postgravity_increment_shear_strain_pct=float(increment[mask].max()*100),max_abs_total_volume_strain_pct=float(np.max(abs(vol))*100),bbar_elastic_bulk_pressure_strain_error_relative=relative,scope='Committed phase snapshots; max principal strain difference is engineering maximum shear strain; no independent G/Gmax curve validation'))
    result=dict(source_run=str(run),kinematic_verification=verification,mesh_coordinates_and_tags_match=True,stage_material_strain_ranges=rows,scope='Offline exact small-strain bbar kinematics from saved nodal displacements, not a new load analysis or experimental strain measurement')
    (out/'strain_ranges.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    with (out/'strain_ranges.csv').open('w',newline='',encoding='utf-8-sig') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
