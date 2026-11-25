"""Offline HTML viewer of real FE geometry and committed nodal translations."""
import json
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
from .mesh import FACE_INDEX

COLORS=['#d8ac64','#639cb3','#bf886a','#82aa80','#ab91c4','#dbb68e','#62bdb2','#c67b93']


def geometry(g,frame=None):
    """Per-material boundary surfaces; internal same-material faces removed.
    Half view uses whole elements on the lower-Y half, not a resampled cut.
    """
    center=(min(p[1] for p in g['nodes'].values())+max(p[1] for p in g['nodes'].values()))/2
    result={}
    def mesh(name,faces,kind='soil'):
        tags=sorted({t for f in faces for t in f});lookup={t:i for i,t in enumerate(tags)}
        triangles=[];edges=set()
        for f in faces:
            q=[lookup[t] for t in f]
            triangles.extend([[q[0],q[1],q[2]],[q[0],q[2],q[3]]])
            edges.update(tuple(sorted((q[i],q[(i+1)%4]))) for i in range(4))
        return dict(name=name,kind=kind,tags=tags,triangles=triangles,edges=sorted(edges))
    for mode in ['full','half','surface']:
        groups={}
        if mode=='surface':
            for f in g['faces']:
                if f['side']=='top':groups.setdefault(f['material'],{})[tuple(sorted(f['nodes']))]=f['nodes']
        else:
            for e in g['elements'].values():
                if mode=='half' and e['center_m'][1]>center:continue
                faces=groups.setdefault(e['material'],{})
                for order in FACE_INDEX.values():
                    f=[e['nodes'][q] for q in order];key=tuple(sorted(f))
                    if key in faces:del faces[key]
                    else:faces[key]=f
        result[mode]=[mesh(name,list(faces.values())) for name,faces in sorted(groups.items()) if faces]
    attachments=[]
    for i,raft in enumerate(g['attachments']):
        faces=[list(ops.eleNodes(t)) for t in raft['shell_elements']]
        attachments.append(mesh(f'기초 판 {i+1}',faces,'raft'))
    if frame:
        attachments.append(dict(name='예제 강재 기둥',kind='frame',
                                segments=[list(p) for p in zip(frame['ports'],frame['top_nodes'])]))
    return result,attachments,center


def visualization_node_tags(g,frame=None):
    modes,attachments,_=geometry(g,frame)
    tags={t for meshes in modes.values() for m in meshes for t in m['tags']}
    for m in attachments:
        tags.update(m.get('tags',[]))
        for pair in m.get('segments',[]):tags.update(pair)
    return sorted(tags)


def viewer_payload(g,*,extra=None,visualization=None):
    extra=extra or {};modes,attachments,center=geometry(g,extra.get('frame'))
    tags=visualization_node_tags(g,extra.get('frame'));index={t:i for i,t in enumerate(tags)}
    for meshes in [*modes.values(),attachments]:
        for m in meshes:
            if 'tags' in m:m['nodes']=[index[t] for t in m.pop('tags')]
            if 'segments' in m:m['segments']=[[index[t] for t in pair] for pair in m['segments']]
    xyz=[list(ops.nodeCoord(t)) for t in tags]
    gravity=[list(g['displacement_reference'].get(t,np.zeros(3))[:3]) for t in tags]
    if visualization is not None:
        lookup={t:i for i,t in enumerate(visualization['node_tags'])}
        if not set(tags)<=set(lookup):raise ValueError('animation is missing visible nodes; record visualization_node_tags before transient')
        indices=[lookup[t] for t in tags]
        u=np.asarray(visualization['displacement_m'],float)[:,indices,:]
        times=np.asarray(visualization['time_s'],float)
        dynamic=np.asarray(visualization['reference_displacement_m'],float)[indices,:]
        if u.shape!=(len(times),len(tags),3) or len(times)<2 or not np.all(np.diff(times)>0):raise ValueError('invalid animation dimensions/times')
    else:
        u=np.array([[ops.nodeDisp(t)[:3] for t in tags]],float);times=np.array([ops.getTime()]);dynamic=np.zeros((len(tags),3))
    if not all(np.all(np.isfinite(a)) for a in [xyz,gravity,u,times,dynamic]):raise ValueError('nonfinite viewer data')
    materials=[]
    for i,name in enumerate(sorted(g['material_volume_fractions'])):
        m=g['materials'][name]
        materials.append(dict(name=name,color=COLORS[i%len(COLORS)],volume_fraction=g['material_volume_fractions'][name],
                              E_MPa=m['E_Pa']/1e6,density=m['density_kg_m3'],law=m['active_law']))
    return dict(schema=1,body=g['config']['body'],stage=g['stage'],node_tags=tags,xyz=xyz,
                geometry=modes,attachments=attachments,cut_y_m=center,materials=materials,
                frames=u.tolist(),times=times.tolist(),gravity=gravity,dynamic_reference=dynamic.tolist(),
                has_dynamic=visualization is not None,integration_steps=(visualization or {}).get('integration_steps',0),
                default_reference='dynamic' if visualization is not None else ('gravity' if g['stage']=='static_loaded' else 'total'),
                metadata=dict(terrain=extra.get('terrain'),input=extra.get('input'),dynamic=extra.get('dynamic'),
                              frequency_guidance=extra.get('frequency_guidance'),config_scope=g['config'].get('mesh_scope',g['config']['assumption'])),
                total_nodes=len(g['nodes']),total_elements=len(g['elements']))


def export_viewer(g,path,*,extra=None,visualization=None):
    from plotly.offline import get_plotlyjs
    from .io import serial
    data=viewer_payload(g,extra=extra,visualization=visualization)
    payload=json.dumps(serial(data),ensure_ascii=False,allow_nan=False,separators=(',',':')).replace('<','\\u003c')
    template=Path(__file__).with_name('viewer.html').read_text()
    # Replace the payload last so material names cannot be interpreted as template markers.
    html=template.replace('__PLOTLY_BUNDLE__',get_plotlyjs()).replace('__GROUND_PAYLOAD__',payload)
    Path(path).write_text(html)
    return dict(file=str(path),frames=len(data['frames']),recorded_nodes=len(data['node_tags']),
                has_dynamic=data['has_dynamic'],offline=True)
