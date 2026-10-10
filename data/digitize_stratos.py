"""Extract real black-line pixels from Guerster & Walter (2017), Fig. 1.

Guide knots select a corridor only; output values always come from image pixels.
No trajectory model, smoothing fit or synthetic observation is used.
"""
from pathlib import Path
import hashlib
import json
import csv
import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'stratos_fig1_source.png'
# Observed printed tick centres in the 1782 x 1102 embedded source image.
XTICKS = [338., 658., 976., 1296., 1614.]
HTICKS_KM = [35., 30., 25., 20., 15.]
Y_ZERO, Y_400 = 977., 114.
# Pixel guides, not data: corridor isolates the solid speed trace from other traces.
GUIDES = np.array([[98,840],[120,752],[178,610],[250,495],[338,397],
                   [420,316],[500,260],[580,220],[658,196],[730,179],
                   [796,161],[845,174],[900,190],[938,190],[1005,229],
                   [1050,282],[1088,332],[1105,333],[1128,365],
                   [1180,400],[1200,407],[1240,454],[1296,490],
                   [1330,514],[1370,557],[1450,604],[1530,650],
                   [1614,698],[1680,726],[1744,752]], dtype=float)


def main():
    im = Image.open(SOURCE).convert('RGB')
    rgb = np.asarray(im).astype(float)
    assert im.size == (1782, 1102)
    hslope, hintercept = np.polyfit(XTICKS, HTICKS_KM, 1)
    mask = (np.max(rgb, axis=2) - np.min(rgb, axis=2) < 12) & (np.max(rgb, axis=2) < 120)
    rows = []
    occluded_columns = []
    overlay = im.copy()
    draw = ImageDraw.Draw(overlay)
    for x in range(100, 1745, 4):
        center = np.interp(x, GUIDES[:,0], GUIDES[:,1])
        lo, hi = int(center-27), int(center+28)
        ys = np.flatnonzero(mask[max(140,lo):min(965,hi),x]) + max(140,lo)
        if not len(ys):
            occluded_columns.append(x)
            continue  # Coloured traces can cover the speed line: no invented fill.
        groups = np.split(ys, np.flatnonzero(np.diff(ys)>1)+1)
        group = min(groups, key=lambda g: abs(np.median(g)-center))
        y = float(np.median(group))
        h = (hslope*x+hintercept)*1000
        v = (Y_ZERO-y)*400/(Y_ZERO-Y_400)
        rows.append(dict(pixel_x=x,pixel_y=y,height_m=h,speed_mps=v,
                         height_reading_halfwidth_m=abs(hslope)*3000,
                         speed_reading_halfwidth_mps=1200/(Y_ZERO-Y_400)))
        draw.ellipse((x-3,y-3,x+3,y+3),outline=(255,0,255),width=2)
    with (ROOT/'stratos_october_digitized.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader();writer.writerows(rows)
    (ROOT.parent/'build'/'revision').mkdir(parents=True,exist_ok=True)
    overlay.save(ROOT.parent/'build'/'revision'/'digitization_overlay.png')
    meta={'source':'Guerster & Walter 2017 Fig.1, PDF page 7',
          'doi':'10.1371/journal.pone.0187798','source_image_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          'extraction':'embedded image xref142 from local file.pdf; grayscale black-line centre within manually inspected corridor',
          'image_size':list(im.size),'x_ticks_px':XTICKS,'h_ticks_km':HTICKS_KM,
          'h_km_per_pixel':float(hslope),'h_km_intercept':float(hintercept),
          'y_zero_px':Y_ZERO,'y_400_px':Y_400,'reading_halfwidth_px':3,
          'uncertainty_note':'Pixel halfwidth is reading resolution, not statistical confidence or instrument accuracy. Paper section IV estimates velocity digitization error about 5%.',
          'guide_pixel_knots':GUIDES.tolist(),'n_points':len(rows),
          'height_range_m':[rows[-1]['height_m'],rows[0]['height_m']],
          'source_conflict':'Peak speed height 28833m in 2017 paper; 27833m in summit report. Keep 2017 value for same-source comparison.',
          'occluded_columns_skipped':occluded_columns,'calibration_use':False}
    (ROOT/'stratos_october_digitization.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({k:meta[k] for k in ['n_points','height_range_m']},ensure_ascii=False))

if __name__=='__main__': main()
