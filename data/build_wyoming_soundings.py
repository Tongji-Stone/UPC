"""Fetch dated US radiosonde tables from the U. Wyoming sounding archive.

The web service exposes the archived station soundings as small per-launch CSVs.
This is a reproducible date subset, not an in-situ Roswell Stratos sounding.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import requests


ROOT = Path(__file__).resolve().parent
URL = "https://weather.arcc.uwyo.edu/wsgi/sounding"
STATIONS = {"72364": "Santa Teresa", "72365": "Albuquerque"}
TIMES = (
    "2012-03-15 12:00:00", "2012-07-25 12:00:00",
    "2012-10-14 12:00:00", "2012-10-15 00:00:00",
)
FIELDS = ["station", "station_name", "nominal_time_utc", "release_time_utc",
          "geopotential_height_m", "pressure_pa", "temperature_k",
          "latitude", "longitude", "wind_speed_mps"]


def main() -> None:
    profiles = []
    metadata = []
    for station, name in STATIONS.items():
        for when in TIMES:
            params = {"datetime": when, "type": "TEXT:CSV", "id": station, "src": "UNKNOWN"}
            response = requests.get(URL, params=params, timeout=90)
            response.raise_for_status()
            raw = response.content
            source_rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
            clean = []
            for item in source_rows:
                try:
                    p = float(item["pressure_hPa"]) * 100.0
                    z = float(item["geopotential height_m"])
                    t = float(item["temperature_C"]) + 273.15
                except (ValueError, KeyError):
                    continue
                if p <= 0 or z < 0 or t < 100:
                    continue
                clean.append({
                    "station": station, "station_name": name,
                    "nominal_time_utc": when.replace(" ", "T") + "Z",
                    "release_time_utc": item["time"].replace(" ", "T") + "Z",
                    "geopotential_height_m": round(z, 1),
                    "pressure_pa": round(p, 1), "temperature_k": round(t, 2),
                    "latitude": item["latitude"].strip(),
                    "longitude": item["longitude"].strip(),
                    "wind_speed_mps": item["wind speed_m/s"].strip(),
                })
            if not clean:
                raise RuntimeError(f"No valid profile from {response.url}")
            profiles.extend(clean)
            metadata.append({
                "station": station, "name": name, "nominal_time_utc": when,
                "source_url": response.url,
                "source_response_sha256": hashlib.sha256(raw).hexdigest(),
                "n_valid_levels": len(clean),
                "max_geopotential_height_m": max(x["geopotential_height_m"] for x in clean),
            })
            print(station, when, len(clean), "levels; top", metadata[-1]["max_geopotential_height_m"])
    with (ROOT / "wyoming_soundings_selected.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(profiles)
    (ROOT / "wyoming_soundings_provenance.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
