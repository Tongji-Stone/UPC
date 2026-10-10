"""Render every final PDF page with Poppler and make inspection contact sheets."""
from pathlib import Path
import subprocess
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'build'/'pages';OUT.mkdir(exist_ok=True)
subprocess.run(['pdftoppm','-r','110','-png',str(ROOT/'main.pdf'),str(OUT/'page')],check=True)
from pypdf import PdfReader
count=len(PdfReader(ROOT/'main.pdf').pages)
pages=[p for p in sorted(OUT.glob('page-*.png')) if int(p.stem.split('-')[-1])<=count]
for start in range(0,len(pages),6):
    group=pages[start:start+6];sheet=Image.new('RGB',(1200,1160),'#dddddd');draw=ImageDraw.Draw(sheet)
    for j,p in enumerate(group):
        im=Image.open(p);im.thumbnail((382,548));x=j%3*400+9;y=j//3*580+25
        sheet.paste(im,(x,y));draw.text((x,y-18),p.stem,fill='black')
    sheet.save(ROOT/'build'/'revision'/f'pages_{start+1:02d}_{start+len(group):02d}.png')
print(f'Rendered {len(pages)} pages')
