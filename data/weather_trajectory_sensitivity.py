"""October Stratos trajectory sensitivity to nearby 12Z radiosonde profiles.

This is a *weather-input sensitivity* calculation, not a validation against
the in-situ Red Bull Stratos radiosonde. The available soundings are roughly
260-273 km from Roswell airport, around seven hours before or five hours after
the jump. The 00Z nominal reports were released around 23Z the prior day.
Above each sounding top and below the station surface, the model's 1976 USSA
is retained. A 2-km smoothstep blends the two density profiles at each end.
"""

from __future__ import annotations

import csv
import importlib.util
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
spec = importlib.util.spec_from_file_location("space_diving_model", PARENT / "model.py")
model = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(model)

EXIT_M = 38_969.4
BLEND_M = 2_000.0


def smoothstep(x: float) -> float:
    q = float(np.clip(x, 0.0, 1.0))
    return q*q*(3.0 - 2.0*q)


def profile(station: str, nominal_time: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    records = []
    with (ROOT / "wyoming_soundings_selected.csv").open(newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            if row["station"] == station and row["nominal_time_utc"] == nominal_time:
                records.append(row)
    grouped = defaultdict(list)
    for row in records:
        grouped[float(row["geopotential_height_m"])].append(row)
    z = np.array(sorted(grouped))
    p = np.array([np.mean([float(r["pressure_pa"]) for r in grouped[v]]) for v in z])
    t = np.array([np.mean([float(r["temperature_k"]) for r in grouped[v]]) for v in z])
    if len(z) < 20 or np.any(p <= 0):
        raise ValueError(f"Unusable sounding {station}")
    return z, p, t


def hybrid_density(h: float, obs: tuple[np.ndarray, np.ndarray, np.ndarray]) -> float:
    z, p, t = obs
    standard = model.atmosphere(h)[2]
    geopotential = model.R_E * h/(model.R_E+h)
    if geopotential <= z[0] or geopotential >= z[-1]:
        return standard
    p_obs = np.exp(np.interp(geopotential, z, np.log(p)))
    t_obs = np.interp(geopotential, z, t)
    observed = p_obs/(model.R_AIR*t_obs)
    weight = smoothstep((geopotential-z[0])/BLEND_M)
    weight *= smoothstep((z[-1]-geopotential)/BLEND_M)
    return (1-weight)*standard + weight*observed


def run(beta: float, obs: tuple[np.ndarray, np.ndarray, np.ndarray] | None):
    def rhs(_time, y):
        h, v = y
        rho = model.atmosphere(h)[2] if obs is None else hybrid_density(h, obs)
        return [-v, model.gravity(h)-rho*v*abs(v)/(2*beta)]

    def ground(_time, y):
        return y[0]

    ground.terminal = True
    ground.direction = -1
    solution = solve_ivp(rhs, (0, 500), (EXIT_M, 0.0), rtol=1e-8, atol=1e-9,
                         max_step=0.25, dense_output=True, events=ground)
    times = np.linspace(0, solution.t[-1], 5001)
    heights, speeds = solution.sol(times)
    i = int(np.argmax(speeds))
    lower = float(times[max(i-2, 0)])
    upper = float(times[min(i+2, len(times)-1)])
    peak = minimize_scalar(lambda x: -float(solution.sol(x)[1]),
                           bounds=(lower, upper), method="bounded",
                           options={"xatol": 1e-8})
    peak_time = float(peak.x)
    peak_height, peak_speed = map(float, solution.sol(peak_time))
    return {
        "peak_speed_mps": peak_speed,
        "peak_time_s": peak_time,
        "peak_geometric_height_m": peak_height,
        "time_s": times, "height_m": heights, "speed_mps": speeds,
    }


def main() -> None:
    beta = float(json.loads((PARENT / "results.json").read_text(encoding="utf-8"))
                 ["fitted_felix_beta_kg_m2"])
    cases = {"USSA 1976": run(beta, None)}
    tops = {}
    for station, name, nominal in (
        ("72364", "Santa Teresa 14 Oct 12Z", "2012-10-14T12:00:00Z"),
        ("72365", "Albuquerque 14 Oct 12Z", "2012-10-14T12:00:00Z"),
        ("72364", "Santa Teresa 15 Oct 00Z", "2012-10-15T00:00:00Z"),
        ("72365", "Albuquerque 15 Oct 00Z", "2012-10-15T00:00:00Z"),
    ):
        obs = profile(station, nominal)
        cases[name] = run(beta, obs)
        tops[name] = float(obs[0][-1])
    baseline = cases["USSA 1976"]
    rows = []
    for name, case in cases.items():
        rows.append({
            "case": name, "exit_geometric_m": EXIT_M,
            "ballistic_coefficient_kg_m2": round(beta, 5),
            "sounding_top_geopotential_m": "" if name == "USSA 1976" else tops[name],
            "blend_width_m": "" if name == "USSA 1976" else BLEND_M,
            "peak_speed_mps": round(case["peak_speed_mps"], 3),
            "peak_time_s": round(case["peak_time_s"], 3),
            "peak_geometric_height_m": round(case["peak_geometric_height_m"], 1),
            "delta_speed_vs_standard_mps": round(case["peak_speed_mps"]-baseline["peak_speed_mps"], 3),
            "delta_time_vs_standard_s": round(case["peak_time_s"]-baseline["peak_time_s"], 3),
        })
    with (ROOT / "weather_trajectory_sensitivity.csv").open("w", newline="", encoding="utf-8") as dest:
        writer = csv.DictWriter(dest, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    sample_t = np.linspace(0, 100, 801)
    base_v = np.interp(sample_t, baseline["time_s"], baseline["speed_mps"])
    ax.axhline(0, color="black", lw=1, label="USSA baseline")
    for name, case in cases.items():
        if name == "USSA 1976":
            continue
        varied_v = np.interp(sample_t, case["time_s"], case["speed_mps"])
        ax.plot(sample_t, varied_v-base_v, label=name)
    ax.axvline(baseline["peak_time_s"], color="gray", ls=":", lw=1)
    ax.set(xlabel="Time after exit (s)", ylabel="Speed change vs USSA (m/s)",
           title="Weather-input sensitivity for 14 Oct 2012 jump")
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(PARENT / "figures" / "weather_trajectory_sensitivity.png", dpi=220)
    plt.close(fig)
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
