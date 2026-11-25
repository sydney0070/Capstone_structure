"""Standalone RC ring foundation; SI, explicit example geometry/reinforcement.

C30/SD400 catalogue prototypes, not experimentally calibrated concrete/rebar.
Rectangular ring beam, not the original inverted-T sulfur concrete section.
Builder does not wipe the domain. Demo springs are removable for 3D-soil coupling.
"""
from pathlib import Path
from dataclasses import dataclass,asdict
import argparse,json,math,sys,time,hashlib
import numpy as np
import openseespy.opensees as ops

ROOT=Path(__file__).resolve().parents[1]
RELEASE=ROOT/'material'/'material_release_20261101'
sys.path.insert(0,str(RELEASE))
from capmat.rc import register_rc_section
from capmat.materials import material

@dataclass
class RingConfig:
    radius_m:float=20.
    segments:int=48
    width_m:float=1.
    depth_m:float=.6
    cover_to_bar_surface_m:float=.05
    bar_diameter_m:float=.02
    bars_per_face:int=4
    hinge_length_m:float=.15
    gravity_m_s2:float=3.71
    subgrade_N_m3:float=2e7
    horizontal_to_vertical_ratio:float=.5
    tag_base:int=600000


def build_ring(config=None,*,elastic=False,demo_supports=False):
    """Return 6DOF ring ports/masses. No wipe, no loads; no fake soil if supports=False.

    Ring X/Y are East/North. Local beam y is horizontal radial; local z is Up.
    Nonlinear HingeRadau end regions with elastic transformed-section interiors.
    Reinforcement consists of top/bottom longitudinal bars; stirrups/confinement
    are NOT modeled. Hinge length is an explicit uncalibrated example parameter.
    """
    c=config or RingConfig()
    vals=[c.radius_m,c.width_m,c.depth_m,c.cover_to_bar_surface_m,c.bar_diameter_m,c.hinge_length_m,c.gravity_m_s2,c.subgrade_N_m3,c.horizontal_to_vertical_ratio]
    if not np.isfinite(vals).all() or min(vals)<=0:raise ValueError('Positive finite geometry/material/support parameters required')
    if type(c.segments)!=int or c.segments<12 or c.segments%12:raise ValueError('segments must be a positive multiple of 12, >=12')
    if type(c.bars_per_face)!=int or c.bars_per_face<2:raise ValueError('At least two bars on each face')
    if type(c.tag_base)!=int or c.tag_base<1:raise ValueError('Positive integer tag_base')
    h=c.depth_m/2-c.cover_to_bar_surface_m-c.bar_diameter_m/2
    w=c.width_m/2-c.cover_to_bar_surface_m-c.bar_diameter_m/2
    if min(h,w)<=0:raise ValueError('Cover/bar diameter do not fit section')
    # register_rc_section: local y dimension is depth_m argument.
    # With vecxz=Up, local y is horizontal/radial, z is vertical.
    bars=[(float(y),float(z),c.bar_diameter_m) for z in [-h,h] for y in np.linspace(-w,w,c.bars_per_face)]
    B=c.tag_base;nodes=[B+i+1 for i in range(c.segments)];eles=[B+1000+i for i in range(c.segments)]
    if set(nodes)&set(ops.getNodeTags()) or set(eles)&set(ops.getEleTags()):raise ValueError('Tag collision; choose a different tag_base')
    chord=2*c.radius_m*math.sin(math.pi/c.segments)
    if 8*c.hinge_length_m>=chord:raise ValueError('HingeRadau requires 8*hinge_length < element chord; revise geometry/hinge assumption')
    ops.model('basic','-ndm',3,'-ndf',6)
    for i,n in enumerate(nodes):
        t=2*math.pi*i/c.segments;ops.node(n,c.radius_m*math.cos(t),c.radius_m*math.sin(t),0.)
    mc,mr=material('concrete_c30_prototype'),material('rebar_sd400_prototype')
    Ec,Es=mc['E_Pa'],mr['params']['E'];Gc=mc['G_Pa']
    a,b=max(c.width_m,c.depth_m),min(c.width_m,c.depth_m)
    J=a*b**3*(1/3-.21*(b/a)*(1-b**4/(12*a**4)))
    GJ=Gc*J
    sec,esec,cm,sm,integ,trans=[B+i for i in [2000,2001,2002,2003,2004,2005]]
    meta=register_rc_section(ops,sec,cm,sm,width_m=c.depth_m,depth_m=c.width_m,bars=bars,ny=24,nz=4,GJ_N_m2=GJ,elastic=elastic)
    EA=Ec*meta['concrete_area_m2']+Es*meta['steel_area_m2']
    EIz=Ec*meta['concrete_Iz_m4']+Es*sum(math.pi*d*d/4*y*y for y,z,d in bars)
    EIy=Ec*meta['concrete_Iy_m4']+Es*sum(math.pi*d*d/4*z*z for y,z,d in bars)
    ops.section('Elastic',esec,Ec,EA/Ec,EIz/Ec,EIy/Ec,Gc,GJ/Gc)
    ops.geomTransf('PDelta',trans,0.,0.,1.)
    ops.beamIntegration('HingeRadau',integ,sec,c.hinge_length_m,sec,c.hinge_length_m,esec)
    for i,e in enumerate(eles):ops.element('forceBeamColumn',e,nodes[i],nodes[(i+1)%len(nodes)],trans,integ,'-iter',40,1e-12)
    # Exact chord lengths of this polygon, not a mixture of arc mass and chord elements.
    masses={n:meta['mass_per_length_kg_m']*chord for n in nodes}
    for n,m in masses.items():ops.mass(n,m,m,m,0.,0.,0.)
    fixed=[];support_eles=[]
    if demo_supports:
        kz=c.subgrade_N_m3*c.width_m*chord;kh=kz*c.horizontal_to_vertical_ratio
        ops.uniaxialMaterial('Elastic',B+3000,kh)
        ops.uniaxialMaterial('ENT',B+3001,kz)
        for i,n in enumerate(nodes):
            f=B+4000+i; e=B+5000+i
            if f in ops.getNodeTags() or e in ops.getEleTags():raise ValueError('Demo support tag collision')
            ops.node(f,*ops.nodeCoord(n));ops.fix(f,1,1,1,1,1,1)
            ops.element('zeroLength',e,f,n,'-mat',B+3000,B+3000,B+3001,'-dir',1,2,3)
            fixed.append(f);support_eles.append(e)
    return dict(nodes=nodes,elements=eles,dome_ports=nodes[::c.segments//12],fixed=fixed,support_elements=support_eles,
                nodal_mass_kg=masses,section=meta,bars_y_z_d_m=bars,chord_m=chord,
                total_mass_kg=sum(masses.values()),config=asdict(c),elastic=elastic,
                integration='HingeRadau',torsion='elastic GJ',
                coupling='6DOF ring ports; caller must supply actual soil/contact; no springs by default')


def run_demo(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,default=20.)
    p.add_argument('--segments',type=int,default=48)
    p.add_argument('--width',type=float,default=1.)
    p.add_argument('--depth',type=float,default=.6)
    p.add_argument('--bar-diameter',type=float,default=.02)
    p.add_argument('--bars-per-face',type=int,default=4)
    p.add_argument('--hinge-length',type=float,default=.15)
    p.add_argument('--vertical-load-per-port',type=float,default=100000.,help='Example downward force per 12 ports [N], not actual dome reactions')
    p.add_argument('--elastic',action='store_true')
    p.add_argument('--no-html',action='store_true')
    p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'results_ring_rc')
    a=p.parse_args(argv)
    if not math.isfinite(a.vertical_load_per_port) or a.vertical_load_per_port<0:p.error('Nonnegative finite port load required')
    out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);started=time.perf_counter()
    status=dict(complete=False,status='running',safe_verdict=False)
    def write_status():(out/'run_status.json').write_text(json.dumps(status,indent=2))
    write_status();ops.wipe()
    try:
        c=RingConfig(radius_m=a.radius,segments=a.segments,width_m=a.width,depth_m=a.depth,
                     bar_diameter_m=a.bar_diameter,bars_per_face=a.bars_per_face,hinge_length_m=a.hinge_length)
        m=build_ring(c,elastic=a.elastic,demo_supports=True)
        loads={n:np.array([0.,0.,-mass*c.gravity_m_s2]) for n,mass in m['nodal_mass_kg'].items()}
        for n in m['dome_ports']:loads[n][2]-=a.vertical_load_per_port
        ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
        for n,f in loads.items():ops.load(n,*f,0.,0.,0.)
        ops.constraints('Transformation');ops.numberer('RCM');ops.system('UmfPack')
        ops.test('NormUnbalance',.001,60);ops.algorithm('Newton');ops.integrator('LoadControl',.05);ops.analysis('Static')
        for i in range(20):
            if ops.analyze(1):raise RuntimeError(f'Static nonconvergence at step {i+1}; no safety verdict')
        ops.reactions()
        applied=np.sum(list(loads.values()),axis=0);reaction=np.sum([ops.nodeReaction(n)[:3] for n in m['fixed']],axis=0)
        moment=sum((np.cross(np.array(ops.nodeCoord(n)),f) for n,f in loads.items()),np.zeros(3))
        r_moment=sum((np.cross(np.array(ops.nodeCoord(n)),np.array(ops.nodeReaction(n)[:3]))+np.array(ops.nodeReaction(n)[3:6]) for n in m['fixed']),np.zeros(3))
        force_error=float(np.linalg.norm(applied+reaction)/max(np.linalg.norm(applied),1.))
        moment_error=float(np.linalg.norm(moment+r_moment)/max(np.linalg.norm(applied)*c.radius_m,1.))
        if max(force_error,moment_error)>1e-6:raise RuntimeError(f'Equilibrium failed: F={force_error}, M={moment_error}')
        xyz=np.array([ops.nodeCoord(n) for n in m['nodes']]);u=np.array([ops.nodeDisp(n) for n in m['nodes']])
        if not np.isfinite(u).all():raise RuntimeError('Nonfinite displacement')
        np.savetxt(out/'ring_displacements.csv',np.column_stack((m['nodes'],xyz,u)),delimiter=',',header='node,x_m,y_m,z_m,ux_m,uy_m,uz_m,rx_rad,ry_rad,rz_rad',comments='')
        forces=np.array([ops.eleResponse(e,'localForce') for e in m['elements']])
        np.savetxt(out/'ring_local_forces.csv',np.column_stack((m['elements'],forces)),delimiter=',',header='element,Ni,Vyi,Vzi,Ti,Myi,Mzi,Nj,Vyj,Vzj,Tj,Myj,Mzj',comments='')
        source=RELEASE/'configs'/'catalog.json'
        result=dict(model=m,materials={n:material(n) for n in ['concrete_c30_prototype','rebar_sd400_prototype']},
                    catalog_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),max_settlement_m=float(-u[:,2].min()),
                    total_vertical_load_N=float(-applied[2]),force_balance_relative_error=force_error,moment_balance_relative_error=moment_error,
                    completed_steps=20,elapsed_s=time.perf_counter()-started,
                    assumptions=['Illustrative rectangular RC ring, not original inverted-T sulfur foundation',
                    '8 longitudinal 20mm bars by default; example layout, NOT designed reinforcement or confinement',
                    'C30/SD400 catalogue nominal prototypes, not experiment-fitted materials',
                    'Elastic torsion; shear failure, anchorage, local bearing and rebar buckling not modeled',
                    'Hinge lengths assumed; softening/mesh sensitivity not resolved',
                    'Demo ENT vertical supports and elastic horizontal springs, not 3D soil',
                    'Port loads are examples, not extracted dome loads; no SAFE or optimum-design claim'])
        (out/'summary.json').write_text(json.dumps(result,indent=2))
        if not a.no_html:
            import plotly.graph_objects as go
            scale=100.;closed=np.r_[np.arange(len(xyz)),0];moved=xyz+scale*u[:,:3]
            fig=go.Figure()
            for label,points in [('Initial ring',xyz),('Deformed x100',moved)]:
                fig.add_trace(go.Scatter3d(x=points[closed,0],y=points[closed,1],z=points[closed,2],mode='lines+markers',name=label))
            fig.update_layout(title='C30/SD400 RC ring — illustrative loads/supports; no safety verdict',scene=dict(aspectmode='data',xaxis_title='East [m]',yaxis_title='North [m]',zaxis_title='Up [m]'))
            fig.write_html(out/'ring_rc.html',include_plotlyjs=True)
        status.update(complete=True,status='complete',elapsed_s=time.perf_counter()-started)
        print(json.dumps({k:result[k] for k in ['max_settlement_m','total_vertical_load_N','force_balance_relative_error','moment_balance_relative_error','elapsed_s']},indent=2))
    except Exception as exc:
        status.update(status='failed',error=str(exc));raise
    finally:write_status();ops.wipe()

if __name__=='__main__':run_demo()
