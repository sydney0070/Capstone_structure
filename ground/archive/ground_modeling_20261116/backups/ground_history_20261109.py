"""All-eight-Gauss-point small-strain diagnostics for valid trilinear bbar bricks.

Work bookkeeping is not whole-model energy conservation, damping calibration,
bearing capacity or a mesh-independent measure of damage. No safety verdict.
"""
from pathlib import Path
import csv,math
import numpy as np
import openseespy.opensees as ops
from ground_modeling import brick_jacobian_determinants

class GroundHistory:
    def __init__(self,g,path=None):
        self.path=Path(path) if path else None
        if self.path and self.path.exists():raise ValueError('fresh history output required')
        self.g=g;self.previous=None;self.reference_elastic=None;self.reference_plastic=None;self.max_incremental_plastic=0.;self.internal_work=0.;self.plastic_work=0.;self.rows=[];self.max_plastic=0.;self.max_abs_strain=0.
        self.layout=[];self.Cinv=[];self.weight=[]
        for e,record in g['elements'].items():
            coords=np.array([g['nodes'][n] for n in record['nodes']])
            # OpenSees BbarBrick enumerates xi,eta,zeta with zeta fastest.
            weights=brick_jacobian_determinants(coords)
            if not np.isfinite(weights).all() or np.min(weights)<=0:raise ValueError('positive finite Jacobians required for history')
            if not math.isclose(sum(weights),record['volume_m3'],rel_tol=1e-9,abs_tol=1e-12):raise ValueError('history/mesh quadrature volume mismatch')
            m=g['materials'][record['material']];G,K=m['G_Pa'],m['K_Pa'];lam=K-2*G/3
            C=np.zeros((6,6));C[:3,:3]=lam;C[np.arange(3),np.arange(3)]+=2*G;C[3:,3:]=np.eye(3)*G
            for ip in range(1,9):self.layout.append((e,ip,record['material']));self.Cinv.append(np.linalg.inv(C));self.weight.append(weights[ip-1])
        self.Cinv=np.array(self.Cinv);self.weight=np.array(self.weight)
    def capture(self):
        stress=[];strain=[]
        for e,ip,name in self.layout:
            s=ops.eleResponse(e,'material',ip,'stress');eps=ops.eleResponse(e,'material',ip,'strain')
            if len(s)!=6 or len(eps)!=6:raise RuntimeError(f'Missing all-GP response e={e} ip={ip}')
            stress.append(s);strain.append(eps)
        s=np.asarray(stress);eps=np.asarray(strain)
        if not np.isfinite(s).all() or not np.isfinite(eps).all():raise RuntimeError('Nonfinite all-GP history')
        ee=np.einsum('nij,nj->ni',self.Cinv,s);ep=eps-ee
        # Voigt strain contains engineering shear; elastic work uses stress dot strain.
        elastic=float(np.sum(.5*np.sum(s*ee,axis=1)*self.weight))
        if self.previous is not None:
            s0,e0,ep0=self.previous;avg=(s+s0)/2
            self.internal_work+=float(np.sum(np.sum(avg*(eps-e0),axis=1)*self.weight))
            self.plastic_work+=float(np.sum(np.sum(avg*(ep-ep0),axis=1)*self.weight))
        else:self.reference_elastic=elastic;self.reference_plastic=ep.copy()
        self.previous=(s,eps,ep)
        epnorm=np.sqrt(np.sum(ep[:,:3]**2,axis=1)+.5*np.sum(ep[:,3:]**2,axis=1));i=int(np.argmax(epnorm))
        dep=ep-self.reference_plastic;dep_norm=np.sqrt(np.sum(dep[:,:3]**2,axis=1)+.5*np.sum(dep[:,3:]**2,axis=1));di=int(np.argmax(dep_norm))
        self.max_incremental_plastic=max(self.max_incremental_plastic,float(dep_norm[di]))
        self.max_plastic=max(self.max_plastic,float(epnorm[i]));self.max_abs_strain=max(self.max_abs_strain,float(np.max(abs(eps))))
        row=dict(max_plastic_strain_increment_norm=float(dep_norm[di]),incremental_plastic_GP_count_above_1e_8=int(sum(dep_norm>1e-8)),element_at_max_increment=self.layout[di][0],ip_at_max_increment=self.layout[di][1],time=float(ops.getTime()),gauss_points=len(self.layout),max_plastic_strain_norm=float(epnorm[i]),plastic_GP_count_above_1e_8=int(sum(epnorm>1e-8)),element_at_max=self.layout[i][0],ip_at_max=self.layout[i][1],abs_strain_max=float(np.max(abs(eps))),internal_work_incremental_J=self.internal_work,plastic_work_incremental_J=self.plastic_work,elastic_energy_change_J=elastic-self.reference_elastic,
            work_decomposition_residual_J=self.internal_work-self.plastic_work-(elastic-self.reference_elastic))
        row.update(ground_volume_m3=float(sum(self.weight)),
            plastic_volume_above_1e_8_m3=float(self.weight[epnorm>1e-8].sum()),
            incremental_plastic_volume_above_1e_8_m3=float(self.weight[dep_norm>1e-8].sum()),
            integrated_plastic_strain_norm_m3=float(epnorm@self.weight),
            integrated_incremental_plastic_strain_norm_m3=float(dep_norm@self.weight),
            volume_mean_incremental_plastic_strain_norm=float(dep_norm@self.weight/sum(self.weight)))
        if self.path:
            with self.path.open('a',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(row))
                if not self.rows:writer.writeheader()
                writer.writerow(row)
        self.rows.append(row);return row
    def save(self,path):
        path=Path(path)
        if self.rows:
            with path.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(self.rows[0]));w.writeheader();w.writerows(self.rows)
        return dict(committed_samples=len(self.rows),gauss_points=len(self.layout),max_plastic_strain_norm=self.max_plastic,max_incremental_plastic_strain_norm=self.max_incremental_plastic,max_abs_strain=self.max_abs_strain,
            final=self.rows[-1] if self.rows else None,scope='all 8 GP of positive-Jacobian trilinear bbarBrick; detJ weights; e_plastic=e-C^-1 sigma for constant elastic moduli; incremental trapezoidal work, not full system energy or damage acceptance')
