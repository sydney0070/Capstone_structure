"""Check visualization against the solver, including capture/no-capture equivalence."""
import json
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import openseespy.opensees as ops
from .verification import compact
from . import build_ground,add_raft,initialize_gravity
from .dynamics import activate_absorbing_boundaries,add_incident_base_velocity,run_transient
from .viewer import viewer_payload,visualization_node_tags,export_viewer
from .io import write_json


def main():
    rows=[];runs=[]
    for capture in [False,True]:
        ops.wipe();g=build_ground(compact('Mars',elastic=True));raft=add_raft(g);initialize_gravity(g)
        activate_absorbing_boundaries(g)
        t=np.arange(21)*.001;v=1e-4*np.sin(np.pi*t/t[-1])**2
        add_incident_base_velocity(g,v,.001)
        probe=raft['nodes'][0];initial=list(ops.nodeDisp(probe)[:3])
        r=run_transient(g,dt_s=.001,steps=20,record_nodes=[probe],
                        visualization_nodes=visualization_node_tags(g) if capture else None,max_visualization_frames=5)
        runs.append(np.array([row['nodes'][probe]['u'] for row in r['rows']]))
    assert np.array_equal(runs[0],runs[1])
    rows.append(dict(name='visual capture does not change solver response',**{'pass':True}))
    vis=r['visualization'];idx=vis['node_tags'].index(probe)
    assert len(vis['time_s'])==5 and vis['time_s'][0]==0 and abs(vis['time_s'][-1]-.020)<1e-12
    assert vis['displacement_m'][0][idx]==initial
    for j,step in enumerate([5,10,15,20],1):
        assert np.array_equal(vis['displacement_m'][j][idx],runs[1][step-1,:3])
    data=viewer_payload(g,visualization=vis);index=data['node_tags'].index(probe)
    assert np.array_equal(data['frames'][-1][index],runs[1][-1,:3])
    assert data['frames'][0]==data['dynamic_reference']
    for meshes in [*data['geometry'].values(),data['attachments']]:
        for m in meshes:
            assert all(0<=i<len(data['xyz']) for i in m.get('nodes',[]))
            assert all(0<=i<len(m['nodes']) for tri in m.get('triangles',[]) for i in tri)
    rows.append(dict(name='initial/final and sampled geometry agree with committed solver nodes',frames=5,**{'pass':True}))
    bad=dict(vis,node_tags=vis['node_tags'][1:])
    try:viewer_payload(g,visualization=bad)
    except ValueError:pass
    else:raise AssertionError('missing visualization node accepted')
    with tempfile.TemporaryDirectory(prefix='ground_viewer_test_') as td:
        out=Path(td)/'test.html';injection='</script><script>not executable</script>'
        export_viewer(g,out,visualization=vis,extra={'input':{'kind':injection},'frequency_guidance':{'pass':np.bool_(True)}})
        html=out.read_text();raw=html.split('<script id="ground-data" type="application/json">')[1].split('</script>')[0]
        decoded=json.loads(raw)
        assert decoded['metadata']['input']['kind']==injection and '\\u003c/script' in raw
        assert 'cdn.plot.ly' not in html[:html.index('<script>')]
        from .cli import main as cli
        with patch('webbrowser.open',return_value=True) as opened:
            assert cli(['--body','Earth','--gravity-only','--no-raft','--output',td])==0
            opened.assert_called_once_with((Path(td)/'ground_viewer.html').resolve().as_uri())
        plain=(Path(td)/'ground_viewer.html').read_text()
        static=json.loads(plain.split('<script id="ground-data" type="application/json">')[1].split('</script>')[0])
        assert not static['has_dynamic'] and len(static['frames'])==1
        assert plain==(Path(td)/'ground_cutaway.html').read_text()
    rows.append(dict(name='offline embedding, safe JSON, static mode and default browser-open dispatch (mocked)',**{'pass':True}))
    out=Path(__file__).resolve().parents[1]/'hybrid_verification'/'viewer_summary.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    write_json(out,dict(cases=rows,all_pass=True,browser_rendering='not verified: local file URL blocked by in-app browser policy'))
    print(json.dumps(rows,ensure_ascii=False,indent=2));return 0


if __name__=='__main__':raise SystemExit(main())
