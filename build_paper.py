"""Stage the existing multi-file paper and use the skill's audited XeLaTeX builder."""
from pathlib import Path
import subprocess,sys,shutil,hashlib,json
ROOT=Path(__file__).resolve().parent
TOOL=Path.home()/'.codex'/'skills'/'math-modeling'/'tools'/'latex'/'scripts'/'latex_paper.py'
STAGE=ROOT/'latex_project'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    STAGE.mkdir(parents=True,exist_ok=True)
    mappings=[(ROOT/'main.tex',STAGE/'main.tex')]
    for folder,pattern in [('generated','*.tex'),('figures','*.pdf')]:
        (STAGE/folder).mkdir(exist_ok=True)
        mappings.extend((p,STAGE/folder/p.name) for p in (ROOT/folder).glob(pattern))
    for src,dst in mappings:shutil.copy2(src,dst)
    # Source and assets remain authoritative in UPCchinese; staging avoids
    # copying the Python virtual environment into the TeX sandbox.
    assert all(sha(src)==sha(dst) for src,dst in mappings)
    result=subprocess.run([sys.executable,'-X','utf8',str(TOOL),'build',str(STAGE/'main.tex'),'--engine','xelatex','--timeout','180','--publish',str(ROOT/'main.pdf'),'--overwrite'],cwd=ROOT)
    if result.returncode:raise SystemExit(result.returncode)
    assert all(sha(src)==sha(dst) for src,dst in mappings)
    binding={'main_pdf_sha256':sha(ROOT/'main.pdf'),'resources':[{'source':str(s.relative_to(ROOT)),'staged':str(d.relative_to(ROOT)),'sha256':sha(s)} for s,d in mappings]}
    (ROOT/'build'/'revision'/'paper_source_binding.json').write_text(json.dumps(binding,ensure_ascii=False,indent=2),encoding='utf-8')
    result=subprocess.run([sys.executable,'-X','utf8',str(TOOL),'validate',str(STAGE/'main.tex'),'--pdf',str(ROOT/'main.pdf'),'--contest','generic','--quality-checks','--questions','q1','--min-image-dpi','300'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (ROOT/'build'/'revision'/'latex_validation.json').write_text(result.stdout+result.stderr,encoding='utf-8')
    print(result.stdout)
    if result.returncode:raise SystemExit(result.returncode)
if __name__=='__main__':main()

