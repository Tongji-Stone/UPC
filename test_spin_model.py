"""Independent analytical and physical checks of the local-onset spin model."""
import unittest
from dataclasses import replace
import numpy as np
from scipy.integrate import solve_ivp
import model as m
from spin_model import FrozenSpin,local_response

class SpinTests(unittest.TestCase):
    def setUp(self):
        self.f=FrozenSpin(.05,520,1.24,1.9,190*(1.9**2+3*.3**2)/12,1.,.05)

    def test_analytic_against_independent_numerical_ode(self):
        f=self.f
        # Expand both torques here rather than calling acceleration().
        def rhs(t,y):
            return [(.5*.05*520**2*1.24*1.9*.05-.25*.05*520*1.24*1.9**2*y[0])/f.inertia_kg_m2]
        sol=solve_ivp(rhs,(0,4),[0.],dense_output=True,rtol=1e-11,atol=1e-12)
        t=np.linspace(0,4,101)
        np.testing.assert_allclose(f.solution(t),sol.sol(t)[0],rtol=1e-9,atol=1e-10)

    def test_zero_environment_bias_and_initial_condition(self):
        self.assertEqual(self.f.solution(0),0)
        for kw in [{'density_kg_m3':0},{'speed_mps':0},{'c_n':0}]:
            np.testing.assert_array_equal(replace(self.f,**kw).solution(np.arange(5)),np.zeros(5))

    def test_short_time_and_undamped_limit(self):
        f=self.f;t=1e-8
        self.assertAlmostEqual(f.solution(t)/t,f.bias_torque_nm/f.inertia_kg_m2,places=5)
        f=replace(f,c_nomega=0)
        np.testing.assert_allclose(f.solution(np.arange(5)),f.bias_torque_nm/f.inertia_kg_m2*np.arange(5))

    def test_steady_limit_and_damping_energy(self):
        f=self.f
        self.assertAlmostEqual(f.solution(100),2*f.speed_mps*f.c_n/(f.length_m*f.c_nomega),places=8)
        for omega in [-10,0,10]:self.assertLessEqual(f.torques(omega)[1]*omega,0)

    def test_reduced_rate_normalization_and_vectorized_response(self):
        f=self.f;omega=3
        expected=-.5*f.density_kg_m3*f.speed_mps**2*f.area_m2*f.length_m*f.c_nomega*(omega*f.length_m/(2*f.speed_mps))
        self.assertAlmostEqual(f.torques(omega)[1],expected)
        p=m.DEFAULT_SPIN
        r=local_response([0,.05],[0,520],1.24,190,p)
        np.testing.assert_allclose(r,[0,f.solution(1)])

    def test_invalid_inputs(self):
        for kw in [{'speed_mps':-1},{'c_nomega':-1},{'inertia_kg_m2':0}]:
            with self.assertRaises(ValueError):replace(self.f,**kw)
        with self.assertRaises(ValueError):self.f.solution(-1)

    def test_entire_trajectory_boundary_and_small_reduced_rate(self):
        a=1.2391716561412864
        def evaluate(h):
            sol=m.freefall(h,190,a,stop_height=2566.8,spin_response=m.DEFAULT_SPIN)
            return sol,m.spin_peak_details(sol)
        _,lo=evaluate(66270);sol,hi=evaluate(66290)
        self.assertLess(lo['radial_g'],5);self.assertGreater(hi['radial_g'],5)
        self.assertLess(hi['reduced_rate'],.02)
        sm=m.summarize(sol)
        self.assertGreater(hi['onset_time_s'],sm['speed_mps_peak_time_s'])
        self.assertAlmostEqual(hi['radial_g'],sm['spin_g_peak'],places=8)
        dense=m.spin_peak_details(sol,.025)
        self.assertAlmostEqual(dense['radial_g'],hi['radial_g'],places=6)

if __name__=='__main__':unittest.main(verbosity=2)
