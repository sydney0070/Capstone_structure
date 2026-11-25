"""Relocation contracts: canonical ownership and CLI imports outside the project."""
from pathlib import Path
import json,os,subprocess,sys,unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))


class Layout(unittest.TestCase):
    def test_legacy_imports_resolve_to_canonical_implementations(self):
        sys.path[:0]=[str(ROOT/'ground/ground_modeling'),str(ROOT/'structure'),str(ROOT/'material/material_release_20261101')]
        import capmat.ssi,capmat.rc,capmat.elastic_rc,capmat.ground
        import ground_hybrid.interface
        import mars_dome_ring_ground_20261107
        import planet_foundation_demo_20261109
        expected=[(capmat.ssi,ROOT/'SSI'),(ground_hybrid.interface,ROOT/'SSI'),
                  (capmat.rc,ROOT/'structure/sections'),(capmat.elastic_rc,ROOT/'structure/sections'),
                  (capmat.ground,ROOT/'ground'),(mars_dome_ring_ground_20261107,ROOT/'integrated'),
                  (planet_foundation_demo_20261109,ROOT/'integrated')]
        for module,folder in expected:
            self.assertTrue(Path(module.__file__).is_relative_to(folder),module.__file__)

    def test_cli_imports_do_not_depend_on_working_directory(self):
        env=dict(os.environ);env.pop('PYTHONPATH',None)
        commands=[('ground/ground_modeling/ground_20261107.py',['--help'],0),
                  ('integrated/mars_ssi_integrated_20261109.py',['--help'],0),
                  ('structure/mars_dome_ring_ground_20261107.py',['--dt','.1','--freqmax','1'],2)]
        for path,args,code in commands:
            with self.subTest(path=path):
                job=subprocess.run([sys.executable,str(ROOT/path),*args],cwd=ROOT.parent.parent,
                    env=env,capture_output=True,text=True,timeout=20)
                self.assertEqual(job.returncode,code,job.stderr)


if __name__=='__main__':unittest.main()
