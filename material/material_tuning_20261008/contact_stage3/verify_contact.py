"""Mechanical checks for unilateral contact and elastic ties, not design approval."""
import unittest
from copy import deepcopy
import numpy as np
from contact_model import default,run
from model import default as linear_default,run as linear_run

def small():
    c=default();c['mesh']['n']=12;return c

class Checks(unittest.TestCase):
    def test_closed_contact_matches_linear_foundation(self):
        c=small();c.update(support='contact_only',pressure_path_kPa=[0.,.01])
        s,h,a=run(c)
        d=linear_default();d['mesh']['n']=12;d.update(support='linear_springs',pressure_kPa=.01)
        r,b=linear_run(d)
        self.assertEqual(s['status'],'COMPLETED_REQUESTED_PATH')
        self.assertEqual(h[-1]['contact_length_fraction'],1.)
        np.testing.assert_allclose(a['u'],b['displacements'],atol=1e-10)

    def test_equal_contact_anchor_stiffness_matches_bilateral(self):
        c=small();c['anchor']['stiffness_per_m']=c['foundation']['k_vertical_per_m']
        c['pressure_path_kPa']=[0.,1.]
        s,h,a=run(c)
        d=linear_default();d['mesh']['n']=12;d.update(support='linear_springs',pressure_kPa=1.)
        r,b=linear_run(d)
        np.testing.assert_allclose(a['u'],b['displacements'],atol=1e-10)
        self.assertAlmostEqual(h[-1]['anchor_tension_kN'],r['required_hold_down_kN'],places=7)

    def test_no_tension_contact_and_no_compression_anchor(self):
        c=small();c.update(planet='moon',horizontal_force_kN=10.)
        s,h,a=run(c)
        self.assertEqual(s['status'],'COMPLETED_REQUESTED_PATH')
        self.assertLess(max(r['constitutive_error_kN'] for r in h),1e-8)
        self.assertLess(max(r['max_contact_tension_kN'] for r in h),1e-10)
        self.assertLess(max(r['max_anchor_compression_kN'] for r in h),1e-10)
        self.assertLess(s['max_force_residual'],1e-8)
        self.assertGreater(s['peak_anchor_tension_kN'],0.)

    def test_elastic_cycle_returns_to_initial_equilibrium(self):
        c=small();c.update(planet='moon',horizontal_force_kN=10.)
        s,h,a=run(c);start=next(r for r in h if r['phase']=='lateral' and r['horizontal_load_factor']==1.)
        for k in ('max_u_m','contact_compression_kN','anchor_tension_kN','max_gap_m'):
            self.assertAlmostEqual(start[k],h[-1][k],places=8)

    def test_no_anchor_global_impossibility_is_explicit(self):
        c=small();c.update(support='contact_only',planet='moon')
        s,h,a=run(c)
        self.assertEqual(s['status'],'NO_STATIC_CONTACT_EQUILIBRIUM')
        self.assertEqual(s['attempted_stop']['reason'],'nonpositive_required_compression')
        self.assertLess(s['last_converged_pressure_kPa'],s['attempted_stop']['pressure_kPa'])

    def test_partial_opening_with_net_compression(self):
        c=small();c.update(support='contact_only',planet='moon',horizontal_force_kN=10.,pressure_path_kPa=[0.,.075])
        s,h,a=run(c)
        self.assertEqual(s['status'],'COMPLETED_REQUESTED_PATH')
        self.assertGreater(h[-1]['contact_length_fraction'],0.)
        self.assertLess(h[-1]['contact_length_fraction'],1.)
        self.assertGreater(h[-1]['max_gap_m'],0.)
        self.assertGreater(h[-1]['contact_compression_kN'],0.)

    def test_anchor_allowable_exceedance_stops(self):
        c=small();c.update(planet='moon');c['anchor']['allowable_tension_per_m']=.1
        s,h,a=run(c)
        self.assertEqual(s['status'],'ANCHOR_ALLOWABLE_EXCEEDED')
        self.assertGreater(s['peak_anchor_utilization'],1.)

    def test_pressure_increment_refinement(self):
        c=small();c.update(planet='moon',horizontal_force_kN=10.,pressure_path_kPa=[0.,.075,.4])
        s,h,a=run(c);d=deepcopy(c);d['max_pressure_step_kPa']/=2
        r,hh,b=run(d)
        np.testing.assert_allclose(a['u'],b['u'],atol=1e-10)
        np.testing.assert_allclose(a['anchor'],b['anchor'],atol=1e-8)

    def test_invalid_configuration_rejected(self):
        c=small();c['pressure_path_kPa']=[1.,2.]
        with self.assertRaises(ValueError):run(c)
        c=small();c['shape']='hyperboloid'
        with self.assertRaises(ValueError):run(c)

if __name__=='__main__':unittest.main(verbosity=2)
