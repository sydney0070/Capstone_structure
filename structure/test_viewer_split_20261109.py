"""Geometry, physical-reference and entry-point regression checks."""
from pathlib import Path
import copy, json, sys, tempfile, unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'ground/ground_modeling'),str(ROOT/'structure')]
from insight_dynamic_viewer import build_exterior_by_material
from ground_20261107 import save_response_viewer


class ViewerSplit(unittest.TestCase):
    def test_material_interfaces_are_closed_when_upper_layer_is_hidden(self):
        # Two stacked cubes: unlike materials retain BOTH interface surfaces.
        data={'elements':{1:dict(nodes=list(range(1,9)),material='a'),
                          2:dict(nodes=list(range(5,13)),material='b')}}
        result=build_exterior_by_material(data,list(range(1,13)))
        self.assertEqual([len(result[k]) for k in ('a','b')],[12,12])
        data['elements'][2]['material']='a'
        self.assertEqual(len(build_exterior_by_material(data,list(range(1,13)))['a']),20)

    def test_sorted_tags_preserve_real_total_and_dynamic_displacements(self):
        name='material</script><script>bad()</script>'
        xyz=np.array([(0,0,0),(1,0,0),(1,1,0),(0,1,0),
                      (0,0,1),(1,0,1),(1,1,1),(0,1,1)],float)
        g=dict(nodes={i+1:p for i,p in enumerate(xyz)},
            elements={1:dict(nodes=list(range(1,9)),material=name)},
            materials={name:dict(E_Pa=3e7,density_kg_m3=1840,active_law='DruckerPrager')},
            config=dict(body='Mars'))
        order=np.array([7,0,6,1,5,2,4,3]);tags=(order+1).tolist()
        base=np.arange(24,dtype=float).reshape(8,3)*1e-5
        total=np.array([base,base+np.array([.001,0.,-.002])])
        original=total.copy()
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'ground.html'
            save_response_viewer(g,tags,xyz[order],total[:,order,:],[0.,.5],path,
                context=dict(title='Ground only',static_only=False))
            data=json.loads(path.with_suffix('.json').read_text())
            self.assertEqual(data['tags'],list(range(1,9)))
            np.testing.assert_allclose(data['disp'],total)
            np.testing.assert_allclose(np.array(data['disp'][1])-data['disp'][0],
                                       np.tile([.001,0.,-.002],(8,1)),atol=1e-18)
            self.assertEqual(data['structures'],[])
            self.assertNotIn(name,path.read_text())
            self.assertEqual(data['materials'][0]['name'],name)
        np.testing.assert_array_equal(total,original)

    def test_old_integrated_import_resolves_to_new_implementation(self):
        import mars_dome_ring_ground_20261107 as old
        import mars_ssi_integrated_20261109 as new
        self.assertIs(old,new)
        self.assertEqual(Path(new.__file__).name,'mars_ssi_integrated_20261109.py')

    def test_ground_cli_excludes_integrated_model_options(self):
        from ground_20261107 import parser
        options={o for a in parser()._actions for o in a.option_strings}
        self.assertIn('--ground-history',options)
        for key in ('--mass-ton','--ring-segments','--interface','--anchor-capacity'):
            self.assertNotIn(key,options)


if __name__=='__main__':unittest.main()
