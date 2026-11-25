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
            key = (item['material'], key)
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
    """Export one offline HTML; displacements remain physical metres in the payload."""
    import json
    from plotly.offline import get_plotlyjs
    if not frames:
        raise ValueError('No displacement frames')
    tags = sorted(model_data['node_coords'])
    coords = np.asarray([model_data['node_coords'][t] for t in tags], dtype=float)
    disp = np.asarray([f['disp'] for f in frames], dtype=float)
    times = np.asarray([f['time'] for f in frames], dtype=float)
    if (coords.shape != (len(tags), 3) or not np.isfinite(coords).all()
            or disp.shape != (len(frames), len(tags), 3) or not np.isfinite(disp).all()):
        raise ValueError('Invalid displacement frames')
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError('Frame times must increase')
    if not np.isfinite(visual_scale) or visual_scale <= 0:
        raise ValueError('visual_scale must be positive')
    triangles = build_exterior_by_material(model_data, tags)
    mats = []
    for name in model_data['material_order']:
        if name not in triangles:
            continue
        material = model_data['materials'][name]
        mats.append(dict(name=name, color=material.get('color', '#999999'),
                         rho=material['rho'], E=material.get('E'), law=material.get('law'), triangles=triangles[name]))
    if not mats:
        raise ValueError('No exterior surface')
    payload = dict(coords=coords.tolist(), tags=tags, disp=disp.tolist(), times=times.tolist(),
                   materials=mats, scale=float(visual_scale), sample_node=int(sample_node),
                   material_model=model_data.get('dynamic_material_model', 'unspecified'),
                   motion=model_data.get('motion_metadata', {}),
                   history=np.asarray(history).tolist(), body=model_data.get('body_name', 'Mars'),
                   context=model_data.get('viewer_context', {}), structures=model_data.get('viewer_structures', []))
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace('<', r'\u003c')
    html = _HTML.replace('__PLOTLY__', get_plotlyjs()).replace('__PAYLOAD__', encoded)
    path = Path(html_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding='utf-8')
    if model_data.get('save_viewer_payload'):
        path.with_suffix('.json').write_text(encoded, encoding='utf-8')
    return path


_HTML = r'''<!doctype html>
<html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>화성 지반 · 동적 응답</title><script>__PLOTLY__</script>
<style>
*{box-sizing:border-box}[hidden]{display:none!important}body{margin:0;background:#edf3f6;color:#1b3040;font:15px system-ui,sans-serif}
main{display:grid;grid-template-columns:310px minmax(0,1fr);gap:18px;padding:20px;max-width:2000px;margin:auto}
aside,article{background:white;border:1px solid #d5e1e9;border-radius:15px;overflow:hidden}aside{padding:24px;height:fit-content}
h1{font-size:23px;margin:0 0 12px}h2{font-size:17px;margin:22px 0 16px}p{color:#647b8d;line-height:1.8}hr{border:0;border-top:1px solid #d5e1e9;margin:25px 0}
label{display:block;margin:12px 0 6px}select,input[type=number],button{font:inherit;border:1px solid #b9cfdc;border-radius:8px;background:white;padding:11px;color:#1b3040}select,input[type=number]{width:100%}button{cursor:pointer}button.primary{background:#147d86;color:white;border-color:#147d86}.buttons{display:flex;gap:8px;margin-top:10px}.buttons button{flex:1}
.material{padding:14px 0;border-bottom:1px solid #e7eef2}.material label{display:flex;gap:8px;align-items:center;overflow-wrap:anywhere}.swatch{width:12px;height:12px;display:inline-block;border-radius:4px;flex-shrink:0}.meta{font-size:12px;color:#647b8d;margin-left:24px}
#view{height:66vh;min-height:470px}.controls{display:flex;gap:12px;align-items:center;padding:18px;border-top:1px solid #d5e1e9}#slider{flex:1;min-width:40px;accent-color:#147d86}#speed{width:135px}#time{min-width:150px;font-variant-numeric:tabular-nums}.note{padding:0 20px;font-size:13px}#chart{height:250px}.footer{padding:16px 20px;border-top:1px solid #d5e1e9;font-size:13px;color:#647b8d}
@media(max-width:900px){main{grid-template-columns:1fr;padding:10px}aside{height:auto}.controls{flex-wrap:wrap}#view{height:55vh}#slider{min-width:130px}}
</style></head><body><main><aside><h1>화성 지반 동적 응답</h1><p>HiRISE DEM × MiniSEED<br>OpenSees 3D 지반 모델</p><hr>
<h2>01 · 시간별 응답</h2><p id="summary"></p><button id="peak">최대 변위 프레임</button><hr>
<h2>02 · 변형 표시</h2><label for="reference">변위 기준</label><select id="reference"><option value="dynamic">동적 해석 시작 이후 변위</option><option value="total">중력 포함 전체 변위</option></select>
<label for="scale">변형 확대 배율</label><input id="scale" type="number" min="0.000001" step="any"><div class="buttons"><button id="real">실제 크기 ×1</button><button id="auto">자동 확대</button></div><p>확대배율은 형상 표시에만 적용됩니다. 마우스로 회전·확대하고 절점에 올려 실제 변위를 확인하세요.</p><hr>
<h2>03 · 재료별 표시</h2><div id="materials"></div><p>재료 색상은 응력·파괴 여부를 나타내지 않습니다.</p></aside>
<article><div id="error" role="alert" style="display:none;padding:18px;color:#a63c28"></div><div id="view"></div><div class="controls"><button id="play" class="primary">▶ 재생</button><button id="reset">처음</button><input id="slider" aria-label="시간 프레임" type="range" min="0" step="1" value="0"><span id="time"></span><select id="speed" aria-label="재생 속도"><option value="0.02">속도 0.02×</option><option value="0.1">속도 0.1×</option><option value="0.5">속도 0.5×</option><option value="1" selected>속도 1×</option><option value="2">속도 2×</option><option value="10">속도 10×</option></select></div>
<p class="note">초기 상태와 저장 프레임을 표시합니다. 프레임 사이의 최대 응답은 생략될 수 있습니다. 아래 그래프는 전체 절점의 저장 프레임별 최대 변위입니다.</p><div id="chart"></div>
<div class="footer">입력: StationXML 기기보정 및 ENU 변환한 SEIS 지표 관측 가속도. 지표 기록을 기반 입력으로 가정한 실험용 응답입니다. 고정 바닥·자유 측면 경계와 가정한 지층 물성을 사용합니다. DEM 기본 위치는 영상 중앙이며 관측소와의 일치를 뜻하지 않습니다.</div></article></main>
<script>
const D=__PAYLOAD__;
const el=id=>document.getElementById(id), N=D.times.length;
let frame=0, playing=false, timer=null, busy=false, pending=false;
el('scale').value=D.scale; el('slider').max=N-1;
const context=D.context||{};
if(context.title){document.title=context.title;document.querySelector('h1').textContent=context.title;}
if(context.subtitle)document.querySelector('aside > p').textContent=context.subtitle;
if(context.scope)document.querySelector('.footer').textContent=context.scope;
if(context.reference_label)el('reference').querySelector('[value=total]').textContent=context.reference_label;
if(context.static_only){el('reference').value='total';el('reference').querySelector('[value=dynamic]').disabled=true;for(const id of ['play','slider','speed'])el(id).disabled=true;}
if(context.geometry_only){
 document.querySelector('aside h2').textContent='모델 형상';
 document.querySelector('aside h2:nth-of-type(2)').hidden=true;
 for(const id of ['peak','reference','scale','real','auto','chart'])el(id).hidden=true;
 document.querySelectorAll('label[for=reference],label[for=scale],aside .buttons,article .controls,article .note').forEach(node=>node.hidden=true);
 el('auto').parentElement.nextElementSibling.textContent='마우스로 회전·확대하고 재료별 표시를 선택해 층을 확인하세요.';
}
function fail(e){pause();el('error').style.display='block';el('error').textContent='표시 오류: '+(e.message||e);}
window.addEventListener('error',e=>fail(e.error||e.message));window.addEventListener('unhandledrejection',e=>fail(e.reason));
const span=Math.max(...[0,1,2].map(a=>Math.max(...D.coords.map(p=>p[a]))-Math.min(...D.coords.map(p=>p[a]))));
const displacement=(f,i)=>D.disp[f][i].map((v,a)=>v-(el('reference').value==='dynamic'?D.disp[0][i][a]:0));
const amplitudeCache=new Map();
function amplitudes(){const ref=el('reference').value;if(!amplitudeCache.has(ref))amplitudeCache.set(ref,D.disp.map((_,f)=>{let max=0;for(let i=0;i<D.coords.length;i++)max=Math.max(max,Math.hypot(...displacement(f,i)));return max*1000}));return amplitudeCache.get(ref)}
function meshData(){const scale=Number(el('scale').value);
 const delta=D.coords.map((p,i)=>displacement(frame,i));
 const pts=D.coords.map((p,i)=>p.map((v,a)=>v+scale*delta[i][a]));
 const traces=D.materials.map((m,k)=>{
 return {type:'mesh3d',x:pts.map(p=>p[0]),y:pts.map(p=>p[1]),z:pts.map(p=>p[2]),
 i:m.triangles.map(t=>t[0]),j:m.triangles.map(t=>t[1]),k:m.triangles.map(t=>t[2]),color:m.color,flatshading:true,name:m.name,visible:el('mat'+k).checked,
 customdata:delta.map((u,i)=>[D.tags[i],Math.hypot(...u)*1000,u[2]*1000]),
 hovertemplate:context.geometry_only?'%{fullData.name}<br>node %{customdata[0]}<br>좌표 %{x:.3f}, %{y:.3f}, %{z:.3f} m<extra></extra>':'%{fullData.name}<br>node %{customdata[0]}<br>표시 좌표 %{x:.3f}, %{y:.3f}, %{z:.3f} m<br>실제 |Δu| %{customdata[1]:.8f} mm<br>실제 Δuz %{customdata[2]:.8f} mm<extra></extra>'};});
 for(const [k,s] of (D.structures||[]).entries()){
  const c=[[],[],[]],custom=[];
  for(const [a,b] of s.segments){for(const i of [a,b]){for(let d=0;d<3;d++)c[d].push(pts[i][d]);custom.push([D.tags[i],Math.hypot(...delta[i])*1000,delta[i][2]*1000]);}for(const v of c)v.push(null);custom.push([null,null,null]);}
  traces.push({type:'scatter3d',mode:'lines',name:s.name,x:c[0],y:c[1],z:c[2],line:{color:s.color,width:5},visible:el('structure'+k).checked,customdata:custom,hovertemplate:'%{fullData.name}<br>node %{customdata[0]}<br>실제 |Δu| %{customdata[1]:.8f} mm<br>실제 Δuz %{customdata[2]:.8f} mm<extra></extra>'});
 }
 return traces;
}
D.materials.forEach((m,k)=>{const row=document.createElement('div');row.className='material';const label=document.createElement('label');const check=document.createElement('input');check.type='checkbox';check.checked=true;check.id='mat'+k;check.onchange=render;const sw=document.createElement('span');sw.className='swatch';sw.style.background=m.color;label.append(check,sw,document.createTextNode(m.name));const meta=document.createElement('div');meta.className='meta';meta.textContent=(m.E?'E '+(m.E/1e6).toLocaleString()+' MPa · ':'')+'ρ '+m.rho.toLocaleString()+' kg/m³';row.append(label,meta);el('materials').append(row)});
for(const [k,s] of (D.structures||[]).entries()){const row=document.createElement('div');row.className='material';const label=document.createElement('label');const cb=document.createElement('input');cb.type='checkbox';cb.checked=true;cb.id='structure'+k;cb.onchange=render;const sw=document.createElement('span');sw.className='swatch';sw.style.background=s.color;label.append(cb,sw,document.createTextNode(s.name));row.append(label);el('materials').append(row);}
el('summary').textContent=context.geometry_only?context.summary:(context.summary||D.material_model)+' · '+N+' 프레임 · '+D.times[0].toFixed(2)+'–'+D.times[N-1].toFixed(2)+' s';
const layout={margin:{l:0,r:0,t:24,b:0},paper_bgcolor:'white',uirevision:'camera',showlegend:false,scene:{aspectmode:'data',xaxis:{title:{text:'동 X [m]'}},yaxis:{title:{text:'북 Y [m]'}},zaxis:{title:{text:'상 Z [m]'}},camera:{eye:{x:1.6,y:-1.8,z:1.1}}}};
async function render(){if(busy){pending=true;return}busy=true;do{pending=false;const scale=Number(el('scale').value);if(!Number.isFinite(scale)||scale<=0)break;
 const amp=amplitudes(),t=D.times[frame];el('slider').value=frame;el('time').textContent=t.toFixed(3)+' s · '+(frame+1)+' / '+N;
 await Promise.all([Plotly.react('view',meshData(),layout,{responsive:true,displaylogo:false}),Plotly.react('chart',[{x:D.times,y:amp,type:'scatter',mode:'lines',line:{color:'#147d86'},hovertemplate:'%{x:.3f} s<br>%{y:.8f} mm<extra></extra>'}],{margin:{l:80,r:25,t:30,b:50},xaxis:{title:'시간 [s]'},yaxis:{title:'최대 |Δu| [mm]'},shapes:[{type:'line',x0:t,x1:t,y0:0,y1:1,yref:'paper',line:{color:'#df784e',width:2}}]},{responsive:true,displaylogo:false})]);
 }while(pending);busy=false;}
function pause(){playing=false;clearTimeout(timer);el('play').textContent='▶ 재생'}
async function tick(){if(!playing)return;if(frame>=N-1){pause();return}const delay=(D.times[frame+1]-D.times[frame])*1000/Number(el('speed').value);timer=setTimeout(async()=>{if(!playing)return;frame++;await render();tick()},delay)}
el('play').onclick=()=>{if(playing){pause();return}if(frame===N-1)frame=0;playing=true;el('play').textContent='❚❚ 일시정지';render();tick()};
el('slider').oninput=()=>{pause();frame=Number(el('slider').value);render()};el('reset').onclick=()=>{pause();frame=0;render()};
el('scale').onchange=render;el('reference').onchange=render;el('real').onclick=()=>{el('scale').value=1;render()};
el('auto').onclick=()=>{const max=Math.max(...amplitudes())/1000;el('scale').value=max>0?Number((span*0.06/max).toPrecision(5)):1;render()};
el('peak').onclick=()=>{pause();const amp=amplitudes();frame=amp.indexOf(Math.max(...amp));render()};el('speed').onchange=()=>{if(playing){clearTimeout(timer);tick()}};
render();
</script></body></html>'''
