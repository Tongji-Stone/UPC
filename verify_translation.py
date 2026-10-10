"""Revision audit replacing the obsolete English/Chinese translation equality test."""
from pathlib import Path
import re,json,hashlib
import pymupdf as fitz
ROOT=Path(__file__).resolve().parent
text=(ROOT/'main.tex').read_text(encoding='utf-8-sig')
assert not any(x in text for x in ['而不是','而非'])
expected=['引言','模型假设','符号说明','模型建立','模型求解与检验','危险性分析','结论','模型评价']
assert re.findall(r'\\section\{([^}]+)\}',text)==expected
r=json.loads((ROOT/'results.json').read_text(encoding='utf-8'))
assert r['calibration']['training_flights']==['2012-03-15','2012-07-25']
assert abs(r['cases']['limit']['recovery_k']-400)<.001
assert r['maximum_height_m']<86000
assert max(r['convergence'][k] for k in r['convergence'] if k!='height_difference_m')<1e-5
assert r['convergence']['height_difference_m']<1
bindings=json.loads((ROOT/'build'/'revision'/'paper_source_binding.json').read_text(encoding='utf-8'))
for b in bindings['resources']:
 for k in ['source','staged']:
  assert hashlib.sha256((ROOT/b[k]).read_bytes()).hexdigest()==b['sha256'],b[k]
assert hashlib.sha256((ROOT/'main.pdf').read_bytes()).hexdigest()==bindings['main_pdf_sha256']
doc=fitz.open(ROOT/'main.pdf');pages=[p.get_text() for p in doc]
assert all(len(t.strip())>15 for t in pages)
assert not any(x in ''.join(pages) for x in ['而不是','而非'])
report={'status':'PASS','pages':len(doc),'sections':expected,'figures':len(re.findall(r'\\begin\{figure\}',text)),'tables':len(re.findall(r'\\begin\{table\}',text)),'equation_environments':len(re.findall(r'\\begin\{(?:equation|align|gather)\}',text)),'pdf_sha256':bindings['main_pdf_sha256']}
(ROOT/'build'/'revision'/'revision_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
