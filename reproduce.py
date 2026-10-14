"""Offline reproduction entry point; all dependencies are ordinary local packages."""
from pathlib import Path
import subprocess,sys,json,hashlib,platform,importlib.metadata
ROOT=Path(__file__).resolve().parent
def run(name,*args):subprocess.run([sys.executable,'-X','utf8',str(ROOT/name),*args],cwd=ROOT,check=True)
def manifest():
 files=[ROOT/'config.json',ROOT/'requirements.txt',ROOT/'build.ps1',*sorted((ROOT/'data').glob('*.py')),*sorted((ROOT/'data').glob('*.csv')),*sorted((ROOT/'data').glob('*.json')),ROOT/'data'/'stratos_fig1_source.png',*sorted(ROOT.glob('*.py')),ROOT/'main.tex',*sorted((ROOT/'appendices').glob('*.py'))]
 obj={'seed':2023,'command':'.\\.venv\\Scripts\\python.exe reproduce.py --compile','python':platform.python_version(),
 'dependencies':{n:importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib','pandas','Pillow','PyMuPDF','pypdf']},
 'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
 'config':json.loads((ROOT/'config.json').read_text(encoding='utf-8')),
 'calculation_runtime_s':json.loads((ROOT/'results.json').read_text(encoding='utf-8'))['runtime_s']}
 (ROOT/'results'/'reproduction_manifest.json').write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
if __name__=='__main__':
 run('data/digitize_stratos.py');run('test_model.py');run('test_spin_model.py');run('run_analysis.py');run('make_figures.py');run('generate_paper_data.py');run('generate_appendix.py');run('verify_appendix.py')
 manifest()
 if '--compile' in sys.argv:run('build_paper.py');run('verify_project.py');run('render_paper.py')
