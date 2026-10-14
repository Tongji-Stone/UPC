"""Standalone Python: print three height limits (km); no third-party imports."""
import math
from bisect import bisect_right

# SI units. US Standard Atmosphere 1976.
RE, G0, R, GAMMA = 6356766., 9.80665, 287.05287, 1.4
CP, RECOVERY = GAMMA * R / (GAMMA - 1), .71 ** (1 / 3)
MASS, REFERENCE_MASS, OPENING = 190., 121.2, 2566.8
INFLATION, DESIGN_SPEED, CANOPY_CD = 8., 5.5, 1.5
LENGTH, RADIUS, HEAD_RADIUS = 1.9, .3, .5
CN, CNW, RESPONSE = .05, 1., 1.
INERTIA = MASS * (LENGTH ** 2 + 3 * RADIUS ** 2) / 12
HB = (0., 11000., 20000., 32000., 47000., 51000., 71000., 84852.)
LAPSE = (-.0065, 0., .001, .0028, 0., -.0028, -.002)
TB, PB = [288.15], [101325.]
for i in range(7):
    dz = HB[i + 1] - HB[i]
    temp = TB[i] + LAPSE[i] * dz
    factor = (math.exp(-G0 * dz / (R * TB[i])) if LAPSE[i] == 0
              else (temp / TB[i]) ** (-G0 / (R * LAPSE[i])))
    TB.append(temp)
    PB.append(PB[i] * factor)

# Stratos Scientific Summit Summary Report, 2013: March, July, October.
# Each pair is release height (m), observed peak speed (m/s).
FLIGHTS = ((21828.3, 163.0333), (29610., 241.4722), (38969.4, 377.1))

def atmosphere(h, load_extension=False):
    top = 110000. if load_extension else 86000.
    if not math.isfinite(h) or h > top + .01:
        raise ValueError("Height outside atmosphere domain")
    h = max(h, 0.)
    z = RE * h / (RE + h)
    i = min(6, bisect_right(HB, z) - 1)
    dz = z - HB[i]
    temp = TB[i] + LAPSE[i] * dz
    factor = (math.exp(-G0 * dz / (R * TB[i])) if LAPSE[i] == 0
              else (temp / TB[i]) ** (-G0 / (R * LAPSE[i])))
    return temp, PB[i] * factor / (R * temp), math.sqrt(GAMMA * R * temp)

def gravity(h):
    return G0 * (RE / (RE + max(h, 0.))) ** 2

def drag_coefficient(mach):
    mach = abs(mach)
    if mach <= .6:
        return .6
    if mach <= 1.1:
        return .6 + .55 * (mach - .6) ** 2
    if mach <= 1.25:
        return .74 - .32 * (mach - 1.1)
    return .692

def minimize(function, left, right, tolerance=1e-7):
    # Golden-section minimization in a bracket.
    ratio = (math.sqrt(5.) - 1.) / 2
    a, b = right - ratio * (right - left), left + ratio * (right - left)
    fa, fb = function(a), function(b)
    while right - left > tolerance:
        if fa < fb:
            right, b, fb = b, a, fa
            a = right - ratio * (right - left)
            fa = function(a)
        else:
            left, a, fa = a, b, fb
            b = left + ratio * (right - left)
            fb = function(b)
    return (a, fa) if fa < fb else (b, fb)

def root(function, left, right, tolerance=.01):
    # Bisection, retaining a sign-changing bracket.
    fl, fr = function(left), function(right)
    if fl == 0:
        return left
    if fr == 0:
        return right
    if fl * fr > 0:
        raise ValueError("Root is not bracketed")
    while right - left > tolerance:
        middle = (left + right) / 2
        fm = function(middle)
        if fm == 0:
            return middle
        if fl * fm <= 0:
            right = middle
        else:
            left, fl = middle, fm
    return (left + right) / 2

def rk4(rhs, t, y, dt):
    # Classical fourth-order Runge-Kutta step for (height, speed).
    k1 = rhs(t, y)
    k2 = rhs(t + dt / 2, tuple(y[j] + dt * k1[j] / 2 for j in (0, 1)))
    k3 = rhs(t + dt / 2, tuple(y[j] + dt * k2[j] / 2 for j in (0, 1)))
    k4 = rhs(t + dt, tuple(y[j] + dt * k3[j] for j in (0, 1)))
    return tuple(y[j] + dt * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) / 6
                 for j in (0, 1))

def interpolate(y0, y1, f0, f1, dt, u):
    # Cubic Hermite dense output from endpoint states and derivatives.
    a, b = 2*u**3 - 3*u**2 + 1, u**3 - 2*u**2 + u
    c, d = -2*u**3 + 3*u**2, u**3 - u**2
    return tuple(a*y0[j] + b*dt*f0[j] + c*y1[j] + d*dt*f1[j]
                 for j in (0, 1))

class Trajectory:
    def __init__(self, rhs, initial):
        self.times = [0.]
        self.states = [initial]
        self.slopes = [rhs(0., initial)]

    def at(self, t):
        if not self.times[0] <= t <= self.times[-1]:
            raise ValueError("Time outside trajectory")
        i = min(bisect_right(self.times, t) - 1, len(self.times) - 2)
        dt = self.times[i + 1] - self.times[i]
        return interpolate(self.states[i], self.states[i + 1],
                           self.slopes[i], self.slopes[i + 1], dt,
                           (t - self.times[i]) / dt)

def integrate(rhs, initial, stop_height, max_step=.1):
    # Adaptive RK4 step doubling; accept the Richardson-corrected state.
    path = Trajectory(rhs, initial)
    dt, t, y = min(.05, max_step), 0., initial
    while t < 3000.:
        dt = min(dt, max_step, 3000. - t)
        coarse = rk4(rhs, t, y, dt)
        half = rk4(rhs, t, y, dt / 2)
        fine = rk4(rhs, t + dt / 2, half, dt / 2)
        correction = tuple((fine[j] - coarse[j]) / 15 for j in (0, 1))
        error = max(abs(correction[j]) /
                    (1e-12 + 1e-11 * max(abs(y[j]), abs(fine[j])))
                    for j in (0, 1))
        if not math.isfinite(error):
            raise RuntimeError("Nonfinite integration state")
        factor = 2. if error == 0 else min(2., max(.2, .9 * error**(-.2)))
        if error > 1:
            dt *= factor
            if dt < 1e-8:
                raise RuntimeError("Integration step too small")
            continue
        new = tuple(fine[j] + correction[j] for j in (0, 1))
        slope = rhs(t + dt, new)
        if new[0] <= stop_height:
            f0 = path.slopes[-1]
            fraction = root(lambda u: interpolate(y, new, f0, slope, dt, u)[0]
                            - stop_height, 0., 1., 1e-11)
            final = interpolate(y, new, f0, slope, dt, fraction)
            path.times.append(t + dt * fraction)
            path.states.append((stop_height, final[1]))
            path.slopes.append(rhs(path.times[-1], path.states[-1]))
            return path
        t, y = t + dt, new
        path.times.append(t)
        path.states.append(y)
        path.slopes.append(slope)
        dt *= factor
    raise RuntimeError("Height event not reached")

def freefall(h0, mass, area, load_extension=False):
    atmosphere(h0, load_extension)
    if h0 <= OPENING:
        raise ValueError("Release must precede parachute opening")
    def rhs(t, y):
        h, v = y
        _, rho, sound = atmosphere(h, load_extension)
        drag = .5 * rho * drag_coefficient(v / sound) * area * v * abs(v)
        return (-v, gravity(h) - drag / mass)
    return integrate(rhs, (h0, 0.), OPENING)

def calibrate():
    def peak(area, height):
        path = freefall(height, REFERENCE_MASS, area)
        k = max(range(len(path.states)), key=lambda i: path.states[i][1])
        lo, hi = max(0, k - 1), min(len(path.times) - 1, k + 1)
        return -minimize(lambda t: -path.at(t)[1],
                         path.times[lo], path.times[hi])[1]
    def objective(area):
        return sum((peak(area, h) / v - 1)**2 for h, v in FLIGHTS)
    return minimize(objective, .3, 3., 1e-9)[0]

def metric(h, v, area, kind):
    temp, rho, sound = atmosphere(h, load_extension=(kind == "load"))
    if kind == "load":
        return .5 * rho * drag_coefficient(v / sound) * area * v**2 / (MASS * G0)
    if kind == "heat":
        return temp + RECOVERY * v**2 / (2 * CP)
    if kind == "spin":
        rate = rho * v * area * LENGTH**2 * CNW / (4 * INERTIA)
        omega = 2 * v * CN / (LENGTH * CNW) * (-math.expm1(-rate * RESPONSE))
        return HEAD_RADIUS * omega**2 / G0
    raise ValueError("Unknown constraint")

def peak_value(path, function):
    count = max(500, int(path.times[-1] / .25))
    times = [path.times[-1] * i / count for i in range(count + 1)]
    times[-1] = path.times[-1]
    values = [function(t) for t in times]
    best = max(values)
    for k in range(1, count):
        if values[k] >= values[k-1] and values[k] >= values[k+1]:
            best = max(best, -minimize(lambda t: -function(t),
                                       times[k-1], times[k+1])[1])
    return best

def maximum(h0, area, kind):
    path = freefall(h0, MASS, area, load_extension=(kind == "load"))
    return peak_value(path, lambda t: metric(*path.at(t), area, kind))

def height_limit(area, kind, limit):
    top = 110000 if kind == "load" else 86000
    high = float(top)
    if maximum(high, area, kind) <= limit:
        raise RuntimeError("No finite boundary inside model domain")
    for low in range(top - 1000, int(OPENING), -1000):
        if maximum(float(low), area, kind) <= limit:
            return root(lambda h: maximum(h, area, kind) - limit,
                        float(low), high)
        high = float(low)
    raise RuntimeError("No feasible bracket found")

def canopy_descent(h_open, v_open, area):
    # Fully inflated drag area, calibrated from sea-level equilibrium.
    rho0 = atmosphere(0.)[1]
    bc = 2 * MASS * G0 / (rho0 * DESIGN_SPEED**2) - .6 * area
    def forces(t, h, v):
        _, rho, sound = atmosphere(h)
        u = min(1., max(0., t / INFLATION))
        fraction = u*u*(3 - 2*u)
        body = .5 * rho * drag_coefficient(v / sound) * area * v * abs(v)
        canopy = .5 * rho * bc * fraction * v * abs(v)
        return body + canopy
    def rhs(t, y):
        h, v = y
        return (-v, gravity(h) - forces(t, h, v) / MASS)
    # One continuous ODE covers inflation, full canopy and landing at h=0.
    path = integrate(rhs, (h_open, v_open), 0., max_step=.25)
    load = peak_value(path, lambda t: forces(t, *path.at(t)) / (MASS * G0))
    return {"trajectory": path, "drag_area_m2": bc,
            "canopy_area_m2": bc / CANOPY_CD,
            "opening_speed_mps": v_open, "opening_load_g": load,
            "landing_speed_mps": path.states[-1][1],
            "inflation_completed": path.times[-1] >= INFLATION,
            "full_inflation_height_m": path.at(INFLATION)[0]
            if path.times[-1] >= INFLATION else None}

def complete_descent(h0, area):
    free = freefall(h0, MASS, area)
    parachute = canopy_descent(*free.states[-1], area)
    return {"freefall": free, "parachute": parachute,
            "total_time_s": free.times[-1] + parachute["trajectory"].times[-1]}

def main():
    area = calibrate()
    heights = [height_limit(area, kind, limit) for kind, limit in
               (("load", 5.), ("heat", 390.), ("spin", 5.))]
    descent = complete_descent(min(heights), area)
    for h in heights:
        print(f"{h / 1000:.5f}")
    return {"area_m2": area, "heights_m": heights, "descent": descent}

if __name__ == "__main__":
    results = main()
