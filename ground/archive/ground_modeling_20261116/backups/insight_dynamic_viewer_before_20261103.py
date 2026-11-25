"""Interactive Plotly viewer for sampled OpenSees FEM displacement frames.

Deformation exaggeration is VISUAL ONLY. This module does not create seismic inputs.
"""
from pathlib import Path
import numpy as np
import plotly.graph_objects as go

FACES = ((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))


def build_exterior_by_material(model_data, tags):
    mapping = {int(tag): i for i, tag in enumerate(tags)}
    registry = {}
    for item in model_data['elements'].values():
        nodes = item['nodes']
        for face in FACES:
            fn = tuple(int(nodes[k]) for k in face)
            key = tuple(sorted(fn))
            if key in registry:
                registry[key] = None
            else:
                registry[key] = (item['material'], fn)
    triangles = {}
    for value in registry.values():
        if value is None:
            continue
        material, quad = value
        ids = [mapping[k] for k in quad]
        triangles.setdefault(material, []).extend(((ids[0],ids[1],ids[2]),(ids[0],ids[2],ids[3])))
    return triangles


def save_dynamic_viewer(model_data, frames, history, sample_node, html_path, visual_scale=100.0):
    """frames=[{'time':seconds, 'disp': ndarray (N,3)}], all nodes ordered by tag."""
    if not frames:
        raise ValueError('No sampled displacement frames')
    node_tags = sorted(model_data['node_coords'])
    coords = np.asarray([model_data['node_coords'][tag] for tag in node_tags], dtype=float)
    expected = (len(node_tags), 3)
    for item in frames:
        if np.asarray(item['disp']).shape != expected:
            raise ValueError('Frame displacement shape does not match node tags')
    tris = build_exterior_by_material(model_data, node_tags)
    if not tris:
        raise ValueError('No exterior surface elements found')
    materials = model_data['materials']
    names = [name for name in model_data['material_order'] if name in tris]
    hist = np.asarray(history, dtype=float)
    if hist.ndim != 2 or hist.shape[1] < 7:
        raise ValueError('History should contain t, ux,uy,uz, ax,ay,az')
    scale = float(visual_scale)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('visual_scale must be positive')
    fig = go.Figure()
    def coords_at(item):
        return coords + scale * np.asarray(item['disp'], dtype=float)
    first = coords_at(frames[0])
    for material in names:
        triangles = np.asarray(tris[material], dtype=int)
        mat = materials[material]
        fig.add_trace(go.Mesh3d(x=first[:,0],y=first[:,1],z=first[:,2],
                 i=triangles[:,0],j=triangles[:,1],k=triangles[:,2],
                 color=mat.get('color','#999999'),opacity=0.92,flatshading=True,
                 name=mat.get('label',material),showscale=False,
                 hovertemplate='X=%{x:.3f} m<br>Y=%{y:.3f} m<br>Z=%{z:.3f} m<extra></extra>'))
    nmesh = len(names)
    fig.add_trace(go.Scatter(x=hist[:,0],y=np.linalg.norm(hist[:,1:4],axis=1)*1000,
        mode='lines',line=dict(color='#147d86',width=2),name='중심점 변위 [mm]',xaxis='x2',yaxis='y2'))
    fig.add_trace(go.Scatter(x=[frames[0]['time']],y=[np.interp(frames[0]['time'],hist[:,0],np.linalg.norm(hist[:,1:4],axis=1)*1000)],
        mode='markers',marker=dict(size=11,color='#d14a30'),name='현재 프레임',xaxis='x2',yaxis='y2'))
    plot_frames=[]
    amp = np.linalg.norm(hist[:,1:4],axis=1)*1000
    for k,record in enumerate(frames):
        p = coords_at(record)
        t=float(record['time'])
        updates=[go.Mesh3d(x=p[:,0],y=p[:,1],z=p[:,2]) for _ in names]
        updates.append(go.Scatter(x=[t],y=[float(np.interp(t,hist[:,0],amp))]))
        plot_frames.append(go.Frame(name=str(k),data=updates,traces=list(range(nmesh))+[nmesh+1],
            layout=go.Layout(annotations=[dict(x=0.5,y=1.04,xref='paper',yref='paper',showarrow=False,
                 text=f'시간 {t:.3f} s  ·  변형 표시 배율 ×{scale:g} (실제 변위가 아님)',font=dict(size=13))])))
    fig.frames=plot_frames
    steps=[dict(method='animate',label=f'{float(rec["time"]):.2f}',args=[[str(i)],
        {'mode':'immediate','frame':{'duration':0,'redraw':True},'transition':{'duration':0}}]) for i,rec in enumerate(frames)]
    fig.update_layout(
      title='InSight 지반 동적 변형 | OpenSees 해석 결과 (실험용)',
      template='plotly_white',height=850,
      scene=dict(domain=dict(x=[0,1],y=[0.30,1]),aspectmode='data',
                 xaxis_title='X [m]',yaxis_title='Y [m]',zaxis_title='Z [m]'),
      xaxis2=dict(domain=[0.11,0.93],anchor='y2',title='시간 [s]'),
      yaxis2=dict(domain=[0.06,0.23],anchor='x2',title='중심점 |u| [mm]'),
      legend=dict(x=0.01,y=0.99),
      margin=dict(l=55,r=15,t=90,b=65),
      annotations=[dict(x=0.5,y=1.04,xref='paper',yref='paper',showarrow=False,
                        text=f'시간 {frames[0]["time"]:.3f} s  ·  변형 표시 배율 ×{scale:g}',font=dict(size=13))],
      sliders=[dict(active=0,steps=steps,x=0.12,y=-0.06,len=0.79,
                    currentvalue=dict(prefix='프레임 시간: ',suffix=' s'))],
      updatemenus=[dict(type='buttons',showactive=False,x=0.01,y=-0.06,direction='right',buttons=[
        dict(label='▶ 재생',method='animate',args=[None,{'fromcurrent':True,
             'frame':{'duration':70,'redraw':True},'transition':{'duration':0}}]),
        dict(label='❚❚ 일시정지',method='animate',args=[[None],{'mode':'immediate',
             'frame':{'duration':0,'redraw':False},'transition':{'duration':0}}])])],
    )
    html_path=Path(html_path)
    html_path.parent.mkdir(parents=True,exist_ok=True)
    fig.write_html(str(html_path),include_plotlyjs=True,full_html=True,auto_open=False)
    return html_path
