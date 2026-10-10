"""Reproducible reduced-order calculations for the UPC space-diving draft.

The atmosphere follows the 1976 U.S. Standard Atmosphere through 86 km.
Flight observations are sparse published summary values, not raw telemetry.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar


ROOT = Path(__file__).resolve().parent
R_E = 6_356_766.0
G0 = 9.80665
R_AIR = 287.05287
GAMMA = 1.4
CP = 1004.5

# Geopotential layer boundaries [m] and temperature slopes [K/m].
HB = np.array([0, 11_000, 20_000, 32_000, 47_000, 51_000, 71_000, 84_852.0])
TB = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946])
L = np.array([-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002, 0])
PB = np.empty(len(HB))
PB[0] = 101_325.0
for i in range(1, len(HB)):
    dz = HB[i] - HB[i - 1]
    if L[i - 1] == 0:
        PB[i] = PB[i - 1] * np.exp(-G0 * dz / (R_AIR * TB[i - 1]))
    else:
        PB[i] = PB[i - 1] * (TB[i] / TB[i - 1]) ** (-G0 / (R_AIR * L[i - 1]))


def atmosphere(h: float) -> tuple[float, float, float, float]:
    """Return T [K], pressure [Pa], density [kg/m3], sound speed [m/s]."""
    h = float(h)
    if h > 86_000:
        raise ValueError("1976 standard-atmosphere implementation ends at 86 km")
    h = max(h, 0.0)
    geop = R_E * h / (R_E + h)
    i = int(np.clip(np.searchsorted(HB, geop, side="right") - 1, 0, len(HB) - 1))
    dz = geop - HB[i]
    temp = TB[i] + L[i] * dz
    if L[i] == 0:
        pressure = PB[i] * np.exp(-G0 * dz / (R_AIR * TB[i]))
    else:
        pressure = PB[i] * (temp / TB[i]) ** (-G0 / (R_AIR * L[i]))
    rho = pressure / (R_AIR * temp)
    sound = np.sqrt(GAMMA * R_AIR * temp)
    return temp, pressure, rho, sound


def gravity(h: float) -> float:
    return G0 * (R_E / (R_E + max(h, 0.0))) ** 2


def freefall(h0: float, beta: float, max_time: float = 500.0):
    """Integrate vertical fall; beta = mass / effective drag area [kg/m2]."""
    if not (0 < h0 <= 86_000 and beta > 0):
        raise ValueError("Requires 0 < exit altitude <= 86 km and beta > 0")
    def rhs(_t, y):
        h, v = y
        rho = atmosphere(h)[2]
        return [-v, gravity(h) - rho * v * abs(v) / (2 * beta)]

    def ground(_t, y):
        return y[0]

    ground.terminal = True
    ground.direction = -1
    sol = solve_ivp(
        rhs, (0, max_time), (h0, 0.0), rtol=1e-8, atol=1e-9,
        max_step=0.25, dense_output=True, events=ground,
    )
    if not sol.success or len(sol.t_events[0]) != 1:
        raise RuntimeError("Freefall integration did not reach the ground")
    return sol


def summarize(sol, beta: float):
    t = np.linspace(0, sol.t[-1], 5001)
    h, v = sol.sol(t)
    rho = np.array([atmosphere(x)[2] for x in h])
    temp = np.array([atmosphere(x)[0] for x in h])
    proper_g = rho * v**2 / (2 * beta * G0)
    imax = np.argmax(v)
    # Continuous initial interval of proper acceleration below 0.1 g.
    idx = np.flatnonzero(proper_g >= 0.1)
    low_g_time = float(t[idx[0]]) if len(idx) else float(t[-1])
    q = 0.5 * rho * v**2
    taw = temp + 0.89 * v**2 / (2 * CP)
    t0 = temp + v**2 / (2 * CP)
    mach = v / np.sqrt(GAMMA * R_AIR * temp)
    sonic_idx = np.flatnonzero(mach >= 1)
    sonic_start = float(t[sonic_idx[0]]) if len(sonic_idx) else None
    sonic_end = float(t[sonic_idx[-1]]) if len(sonic_idx) else None
    return {
        "t": t, "h": h, "v": v, "q": q, "taw": taw, "t0": t0,
        "proper_g": proper_g,
        "v_peak": float(v[imax]), "t_peak": float(t[imax]),
        "h_at_peak": float(h[imax]), "low_g_time": low_g_time,
        "q_peak": float(np.max(q)), "taw_peak": float(np.max(taw)),
        "t0_peak": float(np.max(t0)),
        "mach_peak": float(np.max(mach)),
        "sonic_start_s": sonic_start, "sonic_end_s": sonic_end,
    }


def simulate_chute(h0: float, beta: float, m: float = 190.0,
                   h_deploy: float = 2566.8, cda_chute: float = 100.0,
                   inflation_s: float = 4.0):
    """Smooth inflation model for an illustrative sea-level landing."""
    cda_body = m / beta

    def rhs(_t, y):
        h, v, u = y  # u is time since chute triggering, bounded in dynamics
        rho = atmosphere(h)[2]
        s = np.clip(u / inflation_s, 0, 1)
        cda = cda_body + (cda_chute - cda_body) * (3 * s**2 - 2 * s**3)
        return [-v, gravity(h) - rho * cda * v * abs(v) / (2 * m), 1.0]

    pre = freefall(h0, beta)
    crossing = brentq(lambda t: pre.sol(t)[0] - h_deploy, 0, pre.t[-1])
    v_open = float(pre.sol(crossing)[1])

    def ground(_t, y):
        return y[0]

    ground.terminal = True
    ground.direction = -1
    post = solve_ivp(rhs, (0, 500), (h_deploy, v_open, 0.0),
                     rtol=1e-8, atol=1e-9, max_step=0.025,
                     dense_output=True, events=ground)
    if not post.success or len(post.t_events[0]) != 1:
        raise RuntimeError("Canopy integration did not reach the ground")
    t = np.linspace(0, post.t[-1], 6001)
    h, v, u = post.sol(t)
    s = np.clip(u / inflation_s, 0, 1)
    cda = cda_body + (cda_chute - cda_body) * (3 * s**2 - 2 * s**3)
    rho = np.array([atmosphere(x)[2] for x in h])
    proper_g = rho * cda * v**2 / (2 * m * G0)
    return {
        "v_open": v_open,
        "v_landing": float(v[-1]),
        "peak_proper_g": float(np.max(proper_g)),
        "peak_net_deceleration_g": float(np.max(proper_g - 1)),
        "t_after_opening": float(post.t[-1]),
    }


def main():
    (ROOT / "figures").mkdir(exist_ok=True)
    with (ROOT / "data" / "stratos_summary.csv").open(newline="", encoding="utf-8") as stream:
        observations = list(csv.DictReader(stream))

    # Calibrate to March and July peak speeds, hold October completely aside.
    # Each flight had a suit and varying body posture, so beta is an effective
    # rather than a precisely transferable material constant.
    def loss(beta):
        errs = []
        for row in observations[:2]:
            pred = summarize(freefall(float(row["exit_m"]), beta), beta)["v_peak"]
            errs.append((pred / float(row["peak_mps"]) - 1) ** 2)
        return sum(errs)
    beta = minimize_scalar(loss, bounds=(70, 1000), method="bounded").x

    validation = []
    flight_drag_area = {}
    for row in observations:
        h0 = float(row["exit_m"])
        sim = summarize(freefall(h0, beta), beta)
        validation.append({
            "flight": row["flight"], "exit_m": h0,
            "observed_low_g_s": float(row["low_g_s"]),
            "predicted_low_g_s": sim["low_g_time"],
            "observed_open_g": float(row["opening_g"]),
            "observed_vpeak_mps": float(row["peak_mps"]),
            "predicted_vpeak_mps": sim["v_peak"],
            "predicted_tpeak_s": sim["t_peak"],
            "predicted_hpeak_m": sim["h_at_peak"],
            "predicted_sonic_start_s": sim["sonic_start_s"],
            "predicted_sonic_end_s": sim["sonic_end_s"],
        })
        # Diagnose the spread in effective drag areas after evaluating the
        # held-out flight. These one-point fits are sensitivity scenarios,
        # never calibration inputs to the central 190 kg prediction.
        one_flight_beta = brentq(
            lambda b: summarize(freefall(h0, b), b)["v_peak"]
            - float(row["peak_mps"]),
            70.0, 300.0,
        )
        flight_drag_area[row["flight"]] = 121.2 / one_flight_beta

    # Transfer the calibrated drag area, not the calibrated ballistic
    # coefficient: PLOS ONE reports 121.2 kg for suited Baumgartner.
    cda_felix = 121.2 / beta
    beta_target = 190.0 / cda_felix

    # Eustace used a stabilization drogue. Its *effective* ballistic
    # coefficient is fitted to one FAI datum, with timing left for checking.
    eustace_beta = brentq(
        lambda b: summarize(freefall(41_422.0, b), b)["v_peak"] - 1320 / 3.6,
        70.0, 1000.0,
    )
    eustace = summarize(freefall(41_422.0, eustace_beta), eustace_beta)

    heights = np.arange(40_000, 80_001, 5000)
    sweeps = []
    for h0 in heights:
        sm = summarize(freefall(float(h0), beta_target), beta_target)
        sweeps.append({
            "exit_km": h0 / 1000, "v_peak_mps": sm["v_peak"],
            "q_peak_kpa": sm["q_peak"] / 1000,
            "taw_peak_k": sm["taw_peak"],
            "t0_peak_k": sm["t0_peak"],
            "mach_peak": sm["mach_peak"],
            "low_g_s": sm["low_g_time"],
        })

    scenario = summarize(freefall(60_000.0, beta_target), beta_target)
    chute_cases = {str(tau): simulate_chute(60_000.0, beta_target, inflation_s=tau)
                   for tau in (4.0, 8.0)}

    fig, ax = plt.subplots(figsize=(6.2, 3.5))
    for h0 in (40_000, 50_000, 60_000, 70_000, 80_000):
        sm = summarize(freefall(float(h0), beta_target), beta_target)
        ax.plot(sm["v"], sm["h"] / 1000, label=f"{h0/1000:.0f} km")
    ax.set(xlabel="Downward speed (m s$^{-1}$)", ylabel="Altitude (km)")
    ax.grid(alpha=.25)
    ax.legend(title="Exit altitude", ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "speed_altitude.png", dpi=220)
    plt.close(fig)

    # An empirical *scenario range*, not a statistical confidence interval:
    # the observed flights yield different effective drag areas.
    drag_cases = [
        ("October-only", flight_drag_area["2012-10-14"]),
        ("March+July fit", cda_felix),
        ("July-only", flight_drag_area["2012-07-25"]),
    ]
    sensitivity = []
    fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.2))
    xs = np.arange(40_000, 80_001, 5000, dtype=float)
    for label, cda in drag_cases:
        target_beta_case = 190.0 / cda
        points = [summarize(freefall(h, target_beta_case), target_beta_case)
                  for h in xs]
        qvals = np.array([p["q_peak"] / 1000 for p in points])
        tvals = np.array([p["taw_peak"] for p in points])
        t0vals = np.array([p["t0_peak"] for p in points])
        if not (np.all(np.diff(qvals) > 0) and
                np.all(np.diff(tvals) > 0) and np.all(np.diff(t0vals) > 0)):
            raise RuntimeError("Monotone altitude screen required for interpolation")
        q_cross = float(np.interp(6.0, qvals, xs / 1000))
        t_cross = float(np.interp(400.0, tvals, xs / 1000))
        t0_cross = float(np.interp(400.0, t0vals, xs / 1000))
        sensitivity.append({
            "case": label, "cda_m2": cda,
            "q_6kpa_crossing_km": q_cross,
            "taw_400k_crossing_km": t_cross,
            "t0_400k_crossing_km": t0_cross,
            "joint_screen_km": min(q_cross, t_cross),
            "q60_kpa": float(qvals[4]),
            "taw60_k": float(tvals[4]),
            "t060_k": float(t0vals[4]),
        })
        axs[0].plot(xs / 1000, qvals, "o-", ms=3, label=label)
        axs[1].plot(xs / 1000, tvals, "o-", ms=3, label=label)
        if label == "March+July fit":
            axs[1].plot(xs / 1000, t0vals, "k:", lw=1.4,
                        label="Stagnation proxy, fit")
    axs[0].axhline(6.0, color="0.25", ls="--", lw=1)
    axs[1].axhline(400.0, color="0.25", ls="--", lw=1)
    axs[0].set(xlabel="Exit altitude (km)", ylabel="Peak dynamic pressure (kPa)")
    axs[1].set(xlabel="Exit altitude (km)", ylabel="Peak recovery temperature (K)")
    axs[0].legend(fontsize=7)
    for ax in axs:
        ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "loads_sweep.png", dpi=220)
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.2))
    alt = [a["exit_m"] / 1000 for a in validation]
    comparisons = [
        ("observed_vpeak_mps", "predicted_vpeak_mps",
         "Peak downward speed (m s$^{-1}$)"),
        ("observed_low_g_s", "predicted_low_g_s",
         "Initial interval below 0.1 g (s)"),
    ]
    for ax, (observed, predicted, ylabel) in zip(axs, comparisons):
        ax.plot(alt, [a[observed] for a in validation],
                "ko", label="Published summaries")
        ax.plot(alt, [a[predicted] for a in validation],
                "s--", color="#186399", label="Constant-$\\beta$ model")
        ax.axvspan(35, 41, color="#eeeeee", zorder=-10,
                   label="October speed held out")
        ax.set(xlabel="Exit altitude (km)", ylabel=ylabel)
        ax.grid(alpha=.25)
    axs[0].legend(fontsize=6, loc="upper left")
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "flight_validation.png", dpi=220)
    plt.close(fig)

    output = {"fitted_felix_beta_kg_m2": beta,
              "fitted_felix_cda_m2": cda_felix,
              "target_beta_kg_m2": beta_target,
              "per_flight_effective_cda_m2": flight_drag_area,
              "drag_area_sensitivity": sensitivity,
              "eustace_drogue_beta_kg_m2": eustace_beta,
              "eustace_check": {k: eustace[k] for k in
                                ("v_peak", "t_peak", "h_at_peak",
                                 "sonic_start_s", "sonic_end_s")},
              "validation": validation,
              "sweep": sweeps, "chute_60km": chute_cases,
              "scenario_60km": {k: scenario[k] for k in
                                 ("v_peak", "t_peak", "h_at_peak", "low_g_time",
                                  "q_peak", "taw_peak", "t0_peak",
                                  "mach_peak", "sonic_start_s")}}
    (ROOT / "results.json").write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
