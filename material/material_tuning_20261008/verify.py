"""Implementation checks against analytic laws; not material validation."""
import unittest
import numpy as np
from copy import deepcopy
from presets import get_case
from experiments import simulate
from materials import pile_capacities
from planet_conditions import vertical_stress

class PhysicsChecks(unittest.TestCase):
    def test_elastic_shear_and_reaction(self):
        c=get_case('elastic_soil');r=simulate(c);p=c['params']
        np.testing.assert_allclose(r['y'],p['E']/(2*(1+p['nu']))*r['x'],atol=1e-7)
        np.testing.assert_allclose(r['reaction'],r['y'],atol=1e-7)
        np.testing.assert_allclose(r['sigma_vertical'],-c['test']['sigma_v'],atol=1e-7)

    def test_uniaxial_concrete_peak(self):
        c=get_case('concrete_cover');p=c['params'];r=simulate(c,np.linspace(0,p['epsc0'],81))
        self.assertAlmostEqual(r['y'][-1],p['fpc'],places=5)
        np.testing.assert_allclose(r['reaction'],r['y'],atol=1e-7)

    def test_elastic_shell_and_poisson(self):
        c=get_case('shell_geopolymer');r=simulate(c)
        np.testing.assert_allclose(r['y'],c['params']['E']*r['x'],atol=1e-7)
        np.testing.assert_allclose(r['transverse_strain'],-c['params']['nu']*r['x'],atol=1e-12)

    def test_metal_elastic_limit(self):
        c=get_case('shell_aluminum');r=simulate(c,np.linspace(0,1e-5,6))
        np.testing.assert_allclose(r['y'],c['params']['E']*r['x'],atol=1e-7)

    def test_orthotropic_shell(self):
        c=get_case('shell_cfrp_elastic');r=simulate(c)
        np.testing.assert_allclose(r['y'],c['params']['Ex']*r['x'],rtol=1e-8,atol=1e-7)

    def test_contact_force_and_penetration(self):
        c=get_case('interface_friction');r=simulate(c);p=c['params'];N=c['test']['normal_force']
        np.testing.assert_allclose(r['normal_reaction'],N,atol=1e-8)
        np.testing.assert_allclose(r['normal_displacement'],-N/p['Kn'],atol=1e-12)
        self.assertAlmostEqual(max(r['y']),p['mu']*N,places=7)
        self.assertAlmostEqual(min(r['y']),-p['mu']*N,places=7)
        # Check the whole stick/slip/reversal path, not only its extrema.
        expected=[0.]
        for dx in np.diff(r['x']):
            expected.append(np.clip(expected[-1]+p['Kt']*dx,-p['mu']*N,p['mu']*N))
        np.testing.assert_allclose(r['y'],expected,atol=1e-8)

    def test_gap_no_tension(self):
        c=get_case('normal_gap');r=simulate(c)
        self.assertLessEqual(max(r['y']),1e-8)
        self.assertAlmostEqual(min(r['y']),c['params']['Fy'],places=6)

    def test_dashpot_rate(self):
        c=get_case('radiation_dashpot');c['params']['alpha']=.7;r=simulate(c)
        x=r['x'];expected=c['params']['C']*np.sign(x)*np.abs(x)**.7
        np.testing.assert_allclose(r['y'],expected,atol=1e-8)

    def test_pile_unit_conversion(self):
        p=pile_capacities(100.,.5,20.,1.,200.)
        self.assertAlmostEqual(p['pult'],50.)
        self.assertAlmostEqual(p['tult'],10*np.pi)
        self.assertAlmostEqual(p['qult'],50*np.pi)

    def test_planet_stress(self):
        self.assertAlmostEqual(vertical_stress('moon',1.5,1.),2.43)
        self.assertAlmostEqual(vertical_stress('mars',1.65,1.),6.1215)

    def test_soil_equilibrium(self):
        for name in ('sand_pdmy02','clay_pimy','moon_regolith_dp','mars_regolith_dp'):
            with self.subTest(name=name):
                c=get_case(name);r=simulate(c)
                np.testing.assert_allclose(r['actual_x'],r['x'],atol=1e-10)
                np.testing.assert_allclose(r['reaction'],r['y'],atol=1e-6)
                np.testing.assert_allclose(r['sigma_vertical'],-c['test']['sigma_v'],atol=1e-6)

    def test_pm4_proxy_labels(self):
        for name in ('sand_pm4','silt_pm4'):
            c=get_case(name);r=simulate(c)
            np.testing.assert_allclose(r['actual_x'],r['x'],atol=1e-10)
            self.assertIn('ru_equivalent',r);self.assertNotIn('pore_pressure',r)
            self.assertGreater(np.max(r['ru_equivalent']),.1)

    def test_bad_input_rejected(self):
        c=get_case('concrete_cover');c['params']['epsU']=.003
        with self.assertRaises(ValueError):simulate(c)
        with self.assertRaises(ValueError):simulate(get_case('pile_py'),[.01,.02])

if __name__=='__main__':unittest.main(verbosity=2)
