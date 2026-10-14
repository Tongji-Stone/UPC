"""Compile the authoritative source using local XeLaTeX; no Perl dependency."""
from pathlib import Path
import subprocess,shutil,hashlib,json,re,os,sys
sys.stdout.reconfigure(encoding='utf-8',errors='replace')
ROOT=Path(__file__).resolve().parent
BUILD=ROOT/'build'/'latex';STAGE=ROOT/'build'/'source_bundle'
PDF=ROOT/'main.pdf'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 from generate_appendix import generate
 generate()
 BUILD.mkdir(parents=True,exist_ok=True);STAGE.mkdir(exist_ok=True)
 cmd=['xelatex','-interaction=nonstopmode','-halt-on-error','-file-line-error','-no-shell-escape','-output-directory=build/latex','main.tex']
 env=dict(os.environ,LC_ALL='C',LANG='C')
 output=[]
 for pass_number in range(1,4):
  result=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace')
  output.append(f'--- XeLaTeX pass {pass_number} ---\n'+result.stdout+result.stderr)
  if result.returncode:
   (BUILD/'console.txt').write_text('\n'.join(output),encoding='utf-8')
   print((result.stdout+result.stderr)[-7000:]);raise SystemExit(result.returncode)
 (BUILD/'console.txt').write_text('\n'.join(output),encoding='utf-8')
 log=(BUILD/'main.log').read_text(encoding='utf-8',errors='replace')
 warnings=[x for x in log.splitlines() if any(k in x for k in ['Warning:','Overfull','Missing character:'])]
 if warnings:
  print('\n'.join(warnings));raise RuntimeError('Resolve LaTeX warnings before publishing; see build/latex/main.log')
 source=(ROOT/'main.tex').read_text(encoding='utf-8')
 figures=sorted(set(re.findall(r'\\paperfig\{([^}]+)\}',source)))
 files=[ROOT/'main.tex',ROOT/'config.json',*sorted((ROOT/'appendices').glob('*.py'))]+sorted((ROOT/'generated').glob('*.tex'))+[ROOT/'figures'/(n+'.pdf') for n in figures]
 mappings=[]
 for p in files:
  q=STAGE/p.relative_to(ROOT);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
  mappings.append({'source':str(p.relative_to(ROOT)),'staged':str(q.relative_to(ROOT)),'sha256':sha(p)})
 PDF.parent.mkdir(parents=True,exist_ok=True)
 shutil.copy2(BUILD/'main.pdf',PDF)
 binding={'builder':'local XeLaTeX, three passes','command':cmd,'passes':3,'pdf_file':str(PDF.relative_to(ROOT)),'main_pdf_sha256':sha(PDF),'resources':mappings,'warnings':warnings}
 (ROOT/'build'/'revision').mkdir(parents=True,exist_ok=True)
 for p in [ROOT/'main.build.json',ROOT/'build'/'revision'/'paper_source_binding.json']:p.write_text(json.dumps(binding,ensure_ascii=False,indent=2),encoding='utf-8')
 print('Compiled and published '+str(PDF.relative_to(ROOT))+'; zero LaTeX warnings')
if __name__=='__main__':main()
