"""Recompute October weather sensitivity with fixed calibrated area and Cd(M)."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run_analysis import weather_check,ROOT
if __name__=='__main__':weather_check(json.loads((ROOT/'results.json').read_text(encoding='utf-8'))['calibration']['area_m2'])
