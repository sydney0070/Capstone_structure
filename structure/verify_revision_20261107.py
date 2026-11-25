"""Regression checks with independent kinematic and data-preservation expectations."""
from pathlib import Path
import json,sys,tempfile,subprocess,unittest,copy
import numpy as np
from unittest.mock import patch
HERE=Path(__file__).resolve().parent;GROUND=HERE.parent/'ground'/'ground_modeling'
sys.path.insert(0,str(GROUND));sys.path.insert(0,str(HERE.parent/'material'))
sys.path.insert(0,str(HERE.parent/'material'/'material_release_20261101'))
from ground_20261107 import remove_base_rigid_motion,prepare_motion,station_georeference,require_finite,insight_reference_site
from insight_motion import INSIGHT_ACC_PATH, load_acceleration_csv
import materials_20261107 as audited_materials
from materials_20261107 import audit_catalogs,validate_record

class Revision(unittest.TestCase):
    def test_material_rejection(self):
        from capmat.materials import material
        m=material('steel_s355_coupon');m['params']['Fy']=-1
        with self.assertRaises(ValueError):validate_record('bad',m)
        m=material('mars_jez1_peak');m['G_Pa']*=2
        with self.assertRaises(ValueError):validate_record('bad',m)
        self.assertGreater(audit_catalogs()['checked_records'],20)

    def test_nested_material_conflict(self):
        d=json.loads(audited_materials.RELEASE.read_text())
        with tempfile.TemporaryDirectory() as directory:
            a=Path(directory)/'a.json';b=Path(directory)/'b.json'
            a.write_text(json.dumps(d))
            for change in ('Fy','missing','law'):
                other=copy.deepcopy(d)
                if change=='Fy':other['materials']['steel_s355_coupon']['params']['Fy']*=1.1
                elif change=='missing':del other['materials']['steel_s355_coupon']['params']['R0']
                else:other['materials']['steel_s355_coupon']['law']='UnknownLaw'
                b.write_text(json.dumps(other))
                with patch.object(audited_materials,'RELEASE',a),patch.object(audited_materials,'GROUND',b):
                    with self.assertRaises(ValueError):audit_catalogs()

    def test_short_warmup_requires_opt_in(self):
        with tempfile.TemporaryDirectory() as directory:
            source=INSIGHT_ACC_PATH
            with self.assertRaisesRegex(ValueError,'warmup'):
                prepare_motion(source,directory,start=.5,warmup=5,duration=1)
            _,_,m=prepare_motion(source,directory,start=.5,warmup=5,duration=1,allow_short_warmup=True)
            self.assertTrue(m['warmup_truncated'])
            self.assertAlmostEqual(m['effective_warmup_s'],.5)
            self.assertEqual(m['requested_warmup_s'],5)

    def test_nonfinite_raw_responses_rejected(self):
        for bad in ([1,np.nan],[np.nan,1],[1,np.inf],[[0,1],[2,-np.inf]]):
            with self.assertRaisesRegex(RuntimeError,'Nonfinite'):require_finite(bad,'test')
        np.testing.assert_array_equal(require_finite([1,2],'test'),[1,2])

    def test_contact_nan_rejected_before_extrema(self):
        import ground_hybrid.interface as interface
        contact=dict(element_tag=1,ground_node=2,foundation_node=3,auxiliary_nodes=[4,5])
        with patch.object(interface.ops,'eleForce',return_value=[0.]*12),patch.object(interface.ops,'nodeDisp',return_value=[0.,np.nan,0.]):
            with self.assertRaisesRegex(RuntimeError,'Nonfinite contact'):
                interface.contact_summary({'contacts':[contact]})

    def test_deadline_retains_incomplete_status(self):
        with tempfile.TemporaryDirectory() as directory:
            result=subprocess.run([sys.executable,str(HERE/'mars_dome_ring_ground_20261107.py'),
                '--max-runtime','.01','--output-dir',directory],capture_output=True,text=True)
            self.assertEqual(result.returncode,124)
            status=json.loads((Path(directory)/'run_status.json').read_text())
            self.assertEqual(status['status'],'timeout')
            self.assertFalse(status['complete'])

    def test_insight_site_origin_and_relief_preserved(self):
        site,samples=insight_reference_site(60,nx=5,ny=5)
        self.assertAlmostEqual(site['latitude_deg'],4.502384)
        self.assertAlmostEqual(site['longitude_deg'],135.623447)
        self.assertAlmostEqual(samples[12,0],0.)
        self.assertAlmostEqual(samples[12,1],0.)
        self.assertAlmostEqual(samples[12,3],0.)
        self.assertTrue(np.all(samples[:,4]==0))
        self.assertGreater(np.ptp(samples[:,2]),0)
        self.assertAlmostEqual(site['dem_minus_station_elevation_m'],-3.478,places=2)
        self.assertEqual(len(site['source_files']),2)

    def test_finite_rigid_rotation_and_deformation(self):
        q=np.arange(12)*2*np.pi/12;B=np.column_stack((20*np.cos(q),20*np.sin(q),np.zeros(12)));P=np.r_[B,[[0.,0.,20.]]]
        t=.25;R=np.array([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
        C=B@R+[1,2,3];Q=P@R+[1,2,3]
        self.assertLess(np.max(np.abs(remove_base_rigid_motion(B,C,P,Q))),1e-12)
        Q[-1,2]+=.01
        self.assertAlmostEqual(np.linalg.norm(remove_base_rigid_motion(B,C,P,Q)[-1]),.01,places=12)

    def test_target_not_retapered(self):
        with tempfile.TemporaryDirectory() as d:
            dt,a,m=prepare_motion(INSIGHT_ACC_PATH,Path(d),duration=1,warmup=2,freqmax=1)
            _,full=load_acceleration_csv(Path(d)/'motion_full_enu.csv')
            first=round(m['source_start_s']/dt);target=round(m['target_start_model_s']/dt)
            np.testing.assert_allclose(a[target:],full[first+target:first+len(a)],rtol=0,atol=0)
            self.assertAlmostEqual(m['effective_warmup_s'],2.)

    def test_station_location_inside_dem(self):
        d=station_georeference(GROUND/'DTEEC_037262_1845_036761_1845_U01.IMG',GROUND/'stationxml.xb.elyse.2019.151.xml')
        self.assertAlmostEqual(d['latitude_deg'],4.502384)
        self.assertGreater(d['row_0based'],7000);self.assertLess(d['row_0based'],8000)

    def test_unsafe_and_stale_output_guards(self):
        code=str(HERE/'mars_dome_ring_ground_20261107.py')
        for args in [['--dt','.1','--freqmax','1'],['--net-pressure','100000'],['--cover-thickness','2']]:
            r=subprocess.run([sys.executable,code,*args],capture_output=True,text=True)
            self.assertEqual(r.returncode,2);self.assertNotIn('[assembly]',r.stdout)
        with tempfile.TemporaryDirectory() as d:
            sentinel=Path(d)/'old.html';sentinel.write_text('do not overwrite')
            r=subprocess.run([sys.executable,code,'--output-dir',d],capture_output=True,text=True)
            self.assertEqual(r.returncode,2);self.assertEqual(sentinel.read_text(),'do not overwrite')

if __name__=='__main__':unittest.main(verbosity=2)
