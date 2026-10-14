"""Local flat-spin onset response, not accumulated flight rotation (SI).

C_nomega = -C_nr uses reduced rate omega*L/(2*v). Aircraft coefficients
are analogical assumptions, not fitted human aerodynamic derivatives.
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class SpinResponse:
    length_m: float = 1.9
    width_m: float = .6
    head_radius_m: float = .5
    c_n: float = .05
    c_nomega: float = 1.
    response_s: float = 1.

    def __post_init__(self):
        if not all(np.isfinite(x) for x in vars(self).values()):raise ValueError('Finite parameters required')
        if min(self.length_m,self.width_m,self.head_radius_m,self.response_s)<=0 or self.c_nomega<0:
            raise ValueError('Positive geometry/time and nonnegative damping required')

def transverse_inertia(mass,spin):
    if mass<=0:raise ValueError('Positive mass required')
    return mass*(spin.length_m**2+3*(spin.width_m/2)**2)/12

@dataclass(frozen=True)
class FrozenSpin:
    density_kg_m3: float
    speed_mps: float
    area_m2: float
    length_m: float
    inertia_kg_m2: float
    c_nomega: float = 1.
    c_n: float = .05

    def __post_init__(self):
        if not all(np.isfinite(x) for x in vars(self).values()):raise ValueError('Finite parameters required')
        if min(self.area_m2,self.length_m,self.inertia_kg_m2)<=0:raise ValueError('Positive geometry/inertia required')
        if min(self.density_kg_m3,self.speed_mps,self.c_nomega)<0:raise ValueError('Negative environment/damping')

    @property
    def bias_torque_nm(self):
        return .5*self.density_kg_m3*self.speed_mps**2*self.area_m2*self.length_m*self.c_n

    @property
    def damping_nms(self):
        return .25*self.density_kg_m3*self.speed_mps*self.area_m2*self.length_m**2*self.c_nomega

    def torques(self,omega):return self.bias_torque_nm,-self.damping_nms*np.asarray(omega)

    def acceleration(self,omega):
        bias,damping=self.torques(omega)
        return (bias+damping)/self.inertia_kg_m2

    def solution(self,elapsed_s):
        """Exact response from zero angular velocity; elapsed time starts locally."""
        t=np.asarray(elapsed_s,dtype=float)
        if np.any(~np.isfinite(t)) or np.any(t<0):raise ValueError('Nonnegative finite time required')
        if self.damping_nms==0:return self.bias_torque_nm/self.inertia_kg_m2*t
        return self.bias_torque_nm/self.damping_nms*(-np.expm1(-self.damping_nms*t/self.inertia_kg_m2))

def local_response(rho,v,area,mass,spin):
    """Vectorized local response envelope; no translational feedback."""
    rho,v=np.broadcast_arrays(np.asarray(rho,dtype=float),np.asarray(v,dtype=float))
    if np.any(~np.isfinite(rho)) or np.any(~np.isfinite(v)) or np.any(rho<0) or np.any(v<0):raise ValueError('Invalid environment')
    if area<=0:raise ValueError('Positive area required')
    inertia=transverse_inertia(mass,spin)
    if spin.c_nomega==0:
        return .5*rho*v*v*area*spin.length_m*spin.c_n/inertia*spin.response_s
    return 2*v*spin.c_n/(spin.length_m*spin.c_nomega)*(-np.expm1(
        -.25*rho*v*area*spin.length_m**2*spin.c_nomega/inertia*spin.response_s))
