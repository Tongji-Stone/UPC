"""Create the current portable project archive without repository metadata or caches."""
from pathlib import Path
import zipfile,hashlib,json
ROOT=Path(__file__).resolve().parent
EXCLUDE={'.git','.venv','__pycache__','delivery','build','.pytest_cache'}
def main():
 out=ROOT/'delivery';out.mkdir(exist_ok=True)
 files=[p for p in ROOT.rglob('*') if p.is_file() and not any(part in EXCLUDE for part in p.relative_to(ROOT).parts) and p.suffix not in {'.pyc','.aux','.log','.xdv','.toc','.out','.fls','.fdb_latexmk'}]
 with zipfile.ZipFile(out/'UPC_English.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
  for p in sorted(files):archive.write(p,Path('UPC')/p.relative_to(ROOT))
 with zipfile.ZipFile(out/'UPC_English.zip') as archive:assert archive.testzip() is None
 (out/'archive_manifest.json').write_text(json.dumps({'files':len(files),'archive':'UPC_English.zip','sha256':hashlib.sha256((out/'UPC_English.zip').read_bytes()).hexdigest()},indent=2),encoding='utf-8')
 print('Created delivery/UPC_English.zip')
if __name__=='__main__':main()
