"""Small-displacement ELASTIC screening model, kN-m-s-tonne-kPa.
No calibrated planetary ground, buckling, fracture, or contact opening model.
"""
from pathlib import Path
import sys, json, argparse
import numpy as np
import openseespy.opensees as ops
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from presets import get_case
from planet_conditions import GRAVITY

BASE=Path(__file__).resolve().parent

def default():
    return dict(status='DEMONSTRATION_NOT_CALIBRATED',units='kN,m,s,tonne,kPa',
        shape='dome',material='shell_aluminum',planet='earth',support='pinned',
        geometry=dict(radius=5.,rise=2.5,height=6.,neck_radius=3.),
        mesh=dict(n=48),thickness=.02,gravity_factor=1.,pressure_kPa=0.,
        foundation=dict(k_horizontal_per_m=10000.,k_vertical_per_m=20000.),
        material_override={},material_config=None)

def mesh(c):
    g=c['geometry'];n=c['mesh']['n'];a=g['radius']
    if not isinstance(n,int) or n<4 or n%2:raise ValueError('n must be even and >=4')
    if a<=0:raise ValueError('positive radius required')
    xyz=[];quads=[]
    if c['shape']=='dome':
        h=g['rise']
        if not 0<h<a:raise ValueError('dome requires 0 < rise < radius')
        R=(a*a+h*h)/(2*h)
        for v in np.linspace(-1,1,n+1):
            for u in np.linspace(-1,1,n+1):
                x=a*u*np.sqrt(1-v*v/2);y=a*v*np.sqrt(1-u*u/2)
                xyz.append((x,y,np.sqrt(R*R-x*x-y*y)-(R-h)))
        for j in range(n):
            for i in range(n):
                k=j*(n+1)+i;quads.append((k,k+1,k+n+2,k+n+1))
        boundary=[j*(n+1)+i for j in range(n+1) for i in range(n+1)
                  if i in (0,n) or j in (0,n)]
    elif c['shape']=='hyperboloid':
        h=g['height'];b=g['neck_radius'];nt=4*n
        if not h>0 or not 0<b<a:raise ValueError('positive height and 0 < neck < radius required')
        for z in np.linspace(0,h,n+1):
            radius=np.sqrt(b*b+(a*a-b*b)*(2*z/h-1)**2)
            for theta in np.arange(nt)*2*np.pi/nt:
                xyz.append((radius*np.cos(theta),radius*np.sin(theta),z))
        for j in range(n):
            for i in range(nt):
                k=j*nt+i;l=j*nt+(i+1)%nt
                quads.append((k,l,l+nt,k+nt))
        boundary=list(range(nt))
    else:raise ValueError('unknown shape')
    xyz=np.asarray(xyz);quads=np.asarray(quads,dtype=int)
    boundary=sorted(boundary,key=lambda i:np.arctan2(xyz[i,1],xyz[i,0]))
    lengths={i:0. for i in boundary}
    for i,j in zip(boundary,boundary[1:]+boundary[:1]):
        half=np.linalg.norm(xyz[i]-xyz[j])/2
        lengths[i]+=half;lengths[j]+=half
    return xyz,quads,boundary,lengths

def surface_loads(xyz,quads):
    """2x2 Gauss integration: consistent area and reference pressure vectors."""
    areas=np.zeros(len(xyz));vectors=np.zeros_like(xyz)
    for nodes in quads:
        coord=xyz[nodes]
        for s in (-1/np.sqrt(3),1/np.sqrt(3)):
            for t in (-1/np.sqrt(3),1/np.sqrt(3)):
                N=np.array([(1-s)*(1-t),(1+s)*(1-t),(1+s)*(1+t),(1-s)*(1+t)])/4
                ds=np.array([-(1-t),1-t,1+t,-(1+t)])/4
                dt=np.array([-(1-s),-(1+s),1+s,1-s])/4
                cross=np.cross(ds@coord,dt@coord);jac=np.linalg.norm(cross)
                if jac<1e-12:raise ValueError('degenerate shell element')
                areas[nodes]+=N*jac;vectors[nodes]+=N[:,None]*cross
    return areas,vectors

def run(c):
    if c['units']!='kN,m,s,tonne,kPa':raise ValueError('unit mismatch')
    if c['support'] not in ('pinned','linear_springs'):raise ValueError('unknown support')
    material=get_case(c['material'])
    if c.get('material_config'):
        path=Path(c['material_config'])
        if not path.is_absolute():path=BASE.parent/path
        material=json.loads(path.read_text())
        if material['name']!=c['material'] or material['units']!=c['units']:
            raise ValueError('material config name/units mismatch')
    p=material['params'].copy();p.update(c.get('material_override',{}))
    if 'E' not in p or 'nu' not in p:raise ValueError('isotropic E/nu required; CFRP orientation not implemented')
    if min(p['E'],p['rho'],c['thickness'])<=0 or not -1<p['nu']<.5:raise ValueError('invalid material')
    if c['gravity_factor']<0:raise ValueError('negative gravity factor')
    if c['shape']=='hyperboloid' and c['pressure_kPa']!=0:
        raise ValueError('open hyperboloid is not a closed pressure vessel; add end closures first')
    if min(c['foundation'].values())<=0:raise ValueError('positive line stiffness required')
    xyz,quads,boundary,lengths=mesh(c);areas,pvec=surface_loads(xyz,quads)
    mass=areas*p['rho']*c['thickness'];g=GRAVITY[c['planet']]*c['gravity_factor']
    loads=pvec*c['pressure_kPa'];loads[:,2]-=mass*g
    ops.wipe()
    try:
        ops.model('basic','-ndm',3,'-ndf',6)
        for i,pos in enumerate(xyz,1):
            ops.node(i,*pos);ops.mass(i,*([float(mass[i-1])]*3),0.,0.,0.)
        # rho=0 here: nodal masses already assigned; avoid duplicate inertial mass.
        ops.section('ElasticMembranePlateSection',1,p['E'],p['nu'],c['thickness'],0.)
        for e,nodes in enumerate(quads,1):ops.element('ShellMITC4',e,*map(int,nodes+1),1)
        supports=[]
        for j,i in enumerate(boundary):
            node=i+1
            if c['support']=='pinned':
                ops.fix(node,1,1,1,0,0,0);supports.append(node)
            else:
                ground=len(xyz)+j+1;ops.node(ground,*xyz[i]);ops.fix(ground,1,1,1,1,1,1)
                tags=[]
                for d,key in enumerate(('k_horizontal_per_m','k_horizontal_per_m','k_vertical_per_m'),1):
                    tag=3*j+d;tags.append(tag)
                    ops.uniaxialMaterial('Elastic',tag,c['foundation'][key]*lengths[i])
                ops.element('zeroLength',len(quads)+j+1,ground,node,'-mat',*tags,'-dir',1,2,3)
                supports.append(ground)
        ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
        for i,load in enumerate(loads,1):ops.load(i,*map(float,load),0.,0.,0.)
        ops.constraints('Transformation');ops.numberer('RCM');ops.system('UmfPack')
        ops.test('NormDispIncr',1e-10,20);ops.algorithm('Linear')
        ops.integrator('LoadControl',1.);ops.analysis('Static')
        if ops.analyze(1)!=0:raise RuntimeError('structure solve failed')
        u=np.array([ops.nodeDisp(i+1)[:3] for i in range(len(xyz))]);ops.reactions()
        reactions=np.array([ops.nodeReaction(i)[:3] for i in supports])
        if not np.isfinite(u).all():raise RuntimeError('nonfinite displacement')
        residual=np.sum(reactions,axis=0)+np.sum(loads,axis=0)
        moment=np.sum(np.cross(xyz[boundary],reactions),axis=0)+np.sum(np.cross(xyz,loads),axis=0)
        tension=float(np.sum(np.maximum(-reactions[:,2],0.)))
        out=dict(status='ELASTIC_SCREENING_NOT_DESIGN_VALIDATION',opensees=ops.version(),
            shape=c['shape'],material=c['material'],planet=c['planet'],support=c['support'],
            source_material_status=material.get('status','UNKNOWN'),
            elastic_params_used=dict(E=p['E'],nu=p['nu'],rho=p['rho'],thickness=c['thickness']),
            nodes=len(xyz),elements=len(quads),surface_area_m2=float(sum(areas)),
            mass_tonne=float(sum(mass)),weight_kN=float(sum(mass)*g),
            max_displacement_m=float(max(np.linalg.norm(u,axis=1))),
            min_uz_m=float(min(u[:,2])),max_uz_m=float(max(u[:,2])),
            mean_base_uz_m=float(np.mean(u[boundary,2])),
            required_hold_down_kN=tension,
            force_residual_kN=residual.tolist(),moment_residual_kNm=moment.tolist(),
            load_resultant_kN=np.sum(loads,axis=0).tolist(),
            relative_force_residual=float(np.linalg.norm(residual)/max(sum(np.linalg.norm(loads,axis=1)),1e-12)),
            warnings=['No buckling, cracking, plasticity, environmental effects or calibrated planetary ground.'])
        if tension>1e-7:out['warnings'].append('Tension support required: linear springs/pins assume hold-down; unanchored contact invalid.')
        if out['max_displacement_m']>c['thickness']*.2:
            out['warnings'].append('Displacement exceeds 0.2 thickness screening flag; geometric nonlinear check required.')
        return out,dict(xyz=xyz,quads=quads,displacements=u,loads=loads,reactions=reactions,
                        boundary=np.array(boundary),nodal_mass=mass)
    finally:ops.wipe()

def save(c,outdir):
    result,data=run(c);outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    for name,obj in (('config',c),('summary',result)):
        (outdir/f'{name}.json').write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    np.savez_compressed(outdir/'model_arrays.npz',**data)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path);parser.add_argument('--out',type=Path,default=BASE/'single_run')
    parser.add_argument('--write-default',type=Path)
    args=parser.parse_args()
    if args.write_default:args.write_default.write_text(json.dumps(default(),indent=2)+'\n')
    else:print(json.dumps(save(json.loads(args.config.read_text()) if args.config else default(),args.out),indent=2))
