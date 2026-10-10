"""Independent physical baselines and numerical checks; run with unittest."""
import unittest
from unittest.mock import patch
from dataclasses import replace
import numpy as np
from scipy.integrate import solve_ivp
import model as m

class PhysicalTests(unittest.TestCase):
    def test_atmosphere_reference_and_continuity(self):
        self.assertAlmostEqual(m.atmosphere(0)[2],1.225,places=5)
        self.assertAlmostEqual(m.PB[1],22632.06,delta=.04)
        self.assertAlmostEqual(m.atmosphere(86000)[0],186.946,delta=.001)
        for z in m.HB[1:-1]:
            h=m.R_E*z/(m.R_E-z)
            a=np.array(m.atmosphere(h-.001));b=np.array(m.atmosphere(h+.001))
            self.assertLess(np.max(abs(a/b-1)),1e-6)
        for h in [5000,16000,27000,40000,50000,63000,81000]:
            dp=(m.atmosphere(h+.1)[1]-m.atmosphere(h-.1)[1])/.2
            expected=-m.atmosphere(h)[2]*m.gravity(h)
            self.assertLess(abs(dp/expected-1),1e-6)
        with self.assertRaises(ValueError):m.atmosphere(87000)

    def test_cd_endpoints(self):
        self.assertEqual(m.drag_coefficient(.6),.6)
        self.assertAlmostEqual(m.drag_coefficient(1.1),.7375)
        self.assertAlmostEqual(m.drag_coefficient(1.1+1e-9),.74,places=7)
        self.assertAlmostEqual(m.drag_coefficient(1.25),.692)
        self.assertEqual(m.drag_coefficient(8),.692)

    def test_vacuum_energy(self):
        vac=lambda h:(250.,0.,0.,300.)
        sol=m.freefall(60000,190,1.3,stop_height=20000,atmosphere_fn=vac)
        h,v=sol.y[:,-1]
        expected=np.sqrt(2*m.G0*m.R_E**2*(1/(m.R_E+h)-1/(m.R_E+60000)))
        self.assertLess(abs(v/expected-1),1e-8)

    def test_constant_environment_terminal_solution(self):
        # Verify quadratic-drag integration against v=vt*tanh(gt/vt).
        rho=m.atmosphere(0)[2];cd=.6;a=1.3;mass=190.;vt=np.sqrt(2*mass*m.G0/(rho*cd*a))
        sol=solve_ivp(lambda t,y:[m.G0-rho*cd*a*y[0]**2/(2*mass)],(0,60),[0.],rtol=1e-10,atol=1e-11,dense_output=True)
        t=np.linspace(0,60,100)
        self.assertLess(np.max(abs(sol.sol(t)[0]-vt*np.tanh(m.G0*t/vt))),1e-6)
        self.assertAlmostEqual(m.terminal_speed(0,mass,a),vt,places=8)
        uniform=lambda h:(288.15,101325.,rho,340.3)
        with patch.object(m,'gravity',lambda h:m.G0):
            actual=m.freefall(5000,mass,a,stop_height=1000,law='constant',atmosphere_fn=uniform)
        times=np.linspace(0,actual.t[-1],100)
        self.assertLess(np.max(abs(actual.sol(times)[1]-vt*np.tanh(m.G0*times/vt))),1e-5)

    def test_rotation_zero_torque_and_damping(self):
        att=replace(m.DEFAULT_ATTITUDE,yaw_bias=0,damping=0,pitch_stiffness=0,yaw_rate0=.2,pitch_rate0=.1)
        sol=m.freefall(20000,190,1.3,stop_height=5000,attitude=att)
        self.assertLess(np.max(abs(sol.y[3]-.1)),1e-12)
        self.assertLess(np.max(abs(sol.y[5]-.2)),1e-12)
        att=replace(att,damping=.5)
        dwp,dwz=m.angular_rhs(.1,100,1.3,190,0,.1,.2,att)
        self.assertLess(dwp,0);self.assertLess(dwz,0)

    def test_chute_equilibrium_and_exposure(self):
        out=m.chute(2566.8,70,190,1.357)
        self.assertAlmostEqual(out['height_m'][-1],0,places=6)
        self.assertAlmostEqual(out['summary']['landing_mps'],5.5,delta=.002)
        ex=m.exceedance(np.arange(5.),np.array([0.,2.,0.,2.,0.]),1.)
        self.assertAlmostEqual(ex['total_s'],2.)
        self.assertAlmostEqual(ex['longest_s'],1.)

    def test_numerical_convergence(self):
        a=m.summarize(m.freefall(60000,190,1.357,stop_height=2566.8,attitude=m.DEFAULT_ATTITUDE))
        b=m.summarize(m.freefall(60000,190,1.357,stop_height=2566.8,attitude=m.DEFAULT_ATTITUDE,max_step=.1,rtol=1e-10))
        for k in ['speed_mps_peak','drag_g_peak','recovery_k_peak','spin_g_peak']:
            self.assertLess(abs(a[k]/b[k]-1),.001,k)

if __name__=='__main__':unittest.main(verbosity=2)
