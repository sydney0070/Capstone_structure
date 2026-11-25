"""JSON, VTK and self-contained HTML outputs; DEM and FE resolutions kept separate."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import openseespy.opensees as ops
from .config import graded_axis


def serial(value):
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v) for v in value]
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,Path):return str(value)
    return value


def write_json(path,data):
    Path(path).write_text(json.dumps(serial(data),ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def load_dem(config,*,prepare_pad=False):
    from .terrain import native_dem
    return native_dem(config,prepare_pad=prepare_pad)


def export_vtu(g,path):
    tags=sorted(g['nodes']);lookup={t:i for i,t in enumerate(tags)};elements=list(g['elements'].items())
    material_ids={n:i for i,n in enumerate(sorted(g['materials']))}
    xyz=np.array([g['nodes'][t] for t in tags]);u=np.array([ops.nodeDisp(t)[:3] for t in tags])
    ref=np.array([g['displacement_reference'].get(t,np.zeros(3))[:3] for t in tags])
    def array(name,values,components=1,typ='Float64'):
        values=np.asarray(values).ravel();body=' '.join(str(v) for v in values)
        return f'<DataArray type="{typ}" Name="{name}" NumberOfComponents="{components}" format="ascii">{body}</DataArray>'
    xml=['<?xml version="1.0"?>','<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian"><UnstructuredGrid>',
         f'<Piece NumberOfPoints="{len(tags)}" NumberOfCells="{len(elements)}">',
         '<Points>'+array('Points',xyz,3)+'</Points>',
         '<Cells>'+array('connectivity',[[lookup[n] for n in e['nodes']] for _,e in elements],typ='Int64')+
         array('offsets',np.arange(1,len(elements)+1)*8,typ='Int64')+array('types',[12]*len(elements),typ='UInt8')+'</Cells>',
         '<PointData>'+array('displacement_m',u,3)+array('increment_from_gravity_m',u-ref,3)+array('node_tag',tags,typ='Int64')+'</PointData>',
         '<CellData>'+array('material_id',[material_ids[e['material']] for _,e in elements],typ='Int32')+
         array('element_tag',[t for t,e in elements],typ='Int64')+array('volume_m3',[e['volume_m3'] for _,e in elements])+'</CellData>',
         '</Piece></UnstructuredGrid></VTKFile>']
    Path(path).write_text('\n'.join(xml));return material_ids


def export_html(g,path,*,extra=None,visualization=None):
    from .viewer import export_viewer
    return export_viewer(g,path,extra=extra,visualization=visualization)


def export_results(g,out,*,raft=None,extra=None,html=True,visualization=None):
    from .solver import initial_stress_summary
    from .interface import raft_response
    from .materials import catalog
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    write_json(out/'config.json',g['config']);write_json(out/'materials_used.json',g['materials'])
    write_json(out/'initial_stress_or_final_stress.json',initial_stress_summary(g))
    summary=dict(body=g['config']['body'],stage=g['stage'],boundary=g['boundary'],engine=ops.version(),
                 mesh=g['mesh'],nodes=len(g['nodes']),elements=len(g['elements']),soil_mass_kg=g['mass_kg'],
                 soil_volume_m3=g['volume_m3'],material_volume_fractions=g['material_volume_fractions'],
                 gravity=g.get('gravity_result'),extra=extra or {},
                 material_scope='literature references and explicit assumptions; numerical checks are not field calibration')
    if raft:
        summary['raft_response']=raft_response(g,raft)
        write_json(out/'structure_ports.json',dict(units='m N s kg Pa',ndm=3,ndf=6,coordinates='global model reference coordinates',
            ports=[dict(node=t,coord_m=list(ops.nodeCoord(t)),tributary_area_m2=raft['weights_m2'][t],
                        gravity_displacement_m_rad=g['displacement_reference'].get(t,[])) for t in raft['nodes']],
            connection='reuse these nodes for beams/shells; build structure BEFORE initialize_gravity; caller provides section/member mass',
            constraint_handler='Plain; no rigidLink/rigidDiaphragm/chained MP constraints',raft=raft))
    write_json(out/'summary.json',summary)
    ids=export_vtu(g,out/'ground.vtu');write_json(out/'vtk_material_ids.json',ids)
    if html:
        summary['viewer']=export_html(g,out/'ground_viewer.html',extra=extra,visualization=visualization)
        # Preserve existing links as a standalone identical copy.
        import shutil
        shutil.copyfile(out/'ground_viewer.html',out/'ground_cutaway.html')
        write_json(out/'summary.json',summary)
    root=Path(__file__).parent
    paths=[*root.rglob('*.py'),*root.rglob('*.json'),*root.rglob('*.html'),root.parent/'ground_modeling.py']
    write_json(out/'manifest.json',dict(engine=ops.version(),python=sys.version,
               hashes={str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
               materials_source_hash=catalog()['base_sha256']))
    return summary
