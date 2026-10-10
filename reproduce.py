"""Single offline entry point: python reproduce.py [--compile]."""
from pathlib import Path
import subprocess,sys,json,importlib.util,time,hashlib
ROOT=Path(__file__).resolve().parent
SKILL=Path.home()/'.codex'/'skills'/'math-modeling'
def run(*args):
    subprocess.run([sys.executable,'-X','utf8',*map(str,args)],cwd=ROOT,check=True)
def manifest():
    source=SKILL/'references'/'roles'/'编程手'/'scripts'/'repro_manifest.py'
    spec=importlib.util.spec_from_file_location('skill_manifest',source);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    inputs=[ROOT/'config.json',ROOT/'data'/'stratos_summary.csv',ROOT/'data'/'stratos_fig1_source.png',ROOT/'data'/'wyoming_soundings_selected.csv',ROOT/'data'/'source_registry.json']
    inputs += [ROOT.parent/n for n in ['19770009539.pdf','file.pdf','kap3a_QR20_ekstra_Report_Final.pdf']]
    inputs += list(ROOT.glob('*.py'))+[ROOT/'data'/'digitize_stratos.py']
    config=json.loads((ROOT/'config.json').read_text(encoding='utf-8'))
    obj=mod.build_manifest(inputs,2023,config,'.\\.venv\\Scripts\\python.exe reproduce.py',['numpy','scipy','matplotlib','pandas','Pillow','PyMuPDF','pypdf','requests'])
    obj['calculation_runtime_s']=json.loads((ROOT/'results.json').read_text(encoding='utf-8'))['runtime_s']
    (ROOT/'results'/'复现清单.json').write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
if __name__=='__main__':
    run(ROOT/'data'/'digitize_stratos.py')
    run(ROOT/'test_model.py')
    run(ROOT/'run_analysis.py')
    run(ROOT/'make_chinese_figures.py')
    run(ROOT/'generate_paper_data.py')
    run(SKILL/'tools'/'figure'/'scripts'/'check_figure.py',ROOT/'figures'/'*.png',ROOT/'figures'/'*.svg',ROOT/'figures'/'*.pdf','--strict')
    run(SKILL/'references'/'roles'/'编程手'/'scripts'/'figure_audit.py',ROOT/'figures','--questions','q1','--strict')
    manifest()
    if '--compile' in sys.argv:run(ROOT/'build_paper.py')

