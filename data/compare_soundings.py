"""Compare archival New Mexico radiosondes with the draft's standard atmosphere.

The comparison is at fixed *geopotential* heights. Observed pressure is
interpolated logarithmically, temperature linearly, and density is p/(R*T).
The model's atmosphere function takes geometric altitude, so we convert first.
"""

from __future__ import annotations

import csv
import importlib.util
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
spec = importlib.util.spec_from_file_location("space_diving_model", PARENT / "model.py")
model = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(model)

TARGET_KM = (10.0, 20.0, 25.0)
ROSWELL_AIRPORT = (33.2998775, -104.5294031)


def haversine_km(a: tuple[float, float], b: tuple[float, float]) -> float:
    a1, o1, a2, o2 = map(math.radians, (a[0], a[1], b[0], b[1]))
    q = math.sin((a2-a1)/2)**2 + math.cos(a1)*math.cos(a2)*math.sin((o2-o1)/2)**2
    return 2 * 6371.0 * math.asin(math.sqrt(q))


def standard_at_geopotential(z: float) -> tuple[float, float]:
    h = model.R_E * z / (model.R_E - z)
    t, p, rho, _ = model.atmosphere(h)
    return t, rho


def main() -> None:
    data = defaultdict(list)
    with (ROOT / "wyoming_soundings_selected.csv").open(newline="", encoding="utf-8") as src:
        for row in csv.DictReader(src):
            key = (row["station"], row["nominal_time_utc"])
            data[key].append(row)

    outrows = []
    plotlines = []
    for (station, time), levels in sorted(data.items()):
        # A few duplicate geometric levels may occur. Aggregate those only for
        # interpolation; retain all source observations in the clean input CSV.
        bins = defaultdict(list)
        for r in levels:
            bins[float(r["geopotential_height_m"])].append(r)
        z = np.array(sorted(bins))
        p = np.array([np.mean([float(r["pressure_pa"]) for r in bins[x]]) for x in z])
        t = np.array([np.mean([float(r["temperature_k"]) for r in bins[x]]) for x in z])
        if np.any(p <= 0) or np.any(np.diff(z) <= 0):
            raise ValueError(f"Bad profile {station} {time}")
        loc = (float(levels[0]["latitude"]), float(levels[0]["longitude"]))
        distance = haversine_km(ROSWELL_AIRPORT, loc)
        for km in TARGET_KM:
            height = 1000 * km
            if not z[0] <= height <= z[-1]:
                continue
            ti = float(np.interp(height, z, t))
            pi = float(np.exp(np.interp(height, z, np.log(p))))
            rhoi = pi / (model.R_AIR * ti)
            ts, rhos = standard_at_geopotential(height)
            outrows.append({
                "station": station, "station_name": levels[0]["station_name"],
                "nominal_time_utc": time,
                "distance_from_roswell_airport_km": round(distance, 1),
                "geopotential_height_km": km,
                "observed_temperature_k": round(ti, 2),
                "standard_temperature_k": round(ts, 2),
                "temperature_delta_k": round(ti-ts, 2),
                "observed_density_kg_m3": round(rhoi, 7),
                "standard_density_kg_m3": round(rhos, 7),
                "density_relative_error_pct": round(100*(rhoi/rhos-1), 2),
            })
        if time == "2012-10-14T12:00:00Z":
            plotlines.append((levels[0]["station_name"], z, t, p))

    fields = list(outrows[0])
    with (ROOT / "soundings_vs_standard.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(outrows)

    zgrid = np.linspace(1500, 30000, 300)
    tstd = np.array([standard_at_geopotential(x)[0] for x in zgrid])
    rhostd = np.array([standard_at_geopotential(x)[1] for x in zgrid])
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.0), sharey=True)
    axes[0].plot(tstd, zgrid/1000, color="black", label="USSA 1976", lw=1.6)
    axes[1].axvline(0, color="black", lw=1.6, label="USSA 1976")
    for name, z, t, p in plotlines:
        rho = p/(model.R_AIR*t)
        standard = np.array([standard_at_geopotential(x)[1] for x in z])
        axes[0].plot(t, z/1000, label=name, lw=1.2)
        axes[1].plot(100*(rho/standard-1), z/1000, label=name, lw=1.2)
    axes[0].set(xlabel="Temperature (K)", ylabel="Geopotential height (km)", xlim=(195, 310), ylim=(0, 32))
    axes[1].set(xlabel="Density deviation from USSA (%)", xlim=(-60, 60))
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("14 Oct 2012, 12Z soundings near Roswell")
    fig.tight_layout()
    fig.savefig(PARENT / "figures" / "sounding_validation.png", dpi=220)
    plt.close(fig)

    for row in outrows:
        if row["nominal_time_utc"] == "2012-10-14T12:00:00Z":
            print(row)


if __name__ == "__main__":
    main()
