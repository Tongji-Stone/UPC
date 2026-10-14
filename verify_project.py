"""Validate the English manuscript, figures, numerical results, and standalone appendix."""
from pathlib import Path
import ast,hashlib,json,re,subprocess,sys
import numpy as np
import pymupdf as fitz
from check_printed_atmosphere import check as check_atmosphere
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 text=(ROOT/'main.tex').read_text(encoding='utf-8')
 assert all(name in text for name in ('骆欣畅','魏贤东','刘奕韬'))
 assert 'pdfauthor={Luo Xinchang, Wei Xiandong, Liu Yitao}' in text
 assert 'Team number:' not in text and 'Problem A' in text
 assert 'Keywords:' not in text
 sections=re.findall(r'\\section\{([^}]+)\}',text)
 assert sections==['Introduction','Assumptions','Notation','Physical Models','Model Calibration and Results','Hazard Assessment','Conclusions','Model Evaluation','Solution Code']
 assert text.count(r'\VerbatimInput')==1 and 'fontsize=\\normalsize' in text
 assert r'\setmainfont{Times New Roman}' in text and r'\setmonofont{Latin Modern Mono}' in text
 assert r'\setstretch{1.2}' in text and r'\arraystretch}{1.15}' in text and r'\parindent}{2em}' in text
 assert r'\setstretch{1.2}\tableofcontents' in text
 assert 'skydiv' not in text.lower() and 'high-altitude diving' not in text.lower()
 labels=re.findall(r'\\label\{([^}]+)\}',text)
 labels+=re.findall(r'\\paperfig\{[^}]+\}\{.*?\}\{([^}]+)\}',text)
 refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text)
 assert set(refs)<=set(labels),(set(refs)-set(labels))
 assert len(labels)==len(set(labels))
 eq_labels=re.findall(r'\\label\{(eq:[^}]+)\}',text)
 assert all(x in refs for x in eq_labels)
 for x in labels:
  if x.startswith(('fig:','tab:')):assert x in refs,x
 bib=re.findall(r'\\bibitem\{([^}]+)\}',text)
 citations=set(k for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}',text) for k in group.split(','))
 assert set(bib)==citations
 ordered=list(dict.fromkeys(k for g in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}',text) for k in g.split(',')))
 assert bib==ordered
 assert 'http' not in text and r'\url{' not in text
 names=re.findall(r'\\paperfig\{([^}]+)\}',text)
 assert len(names)==14 and len(re.findall(r'\\begin\{tikzpicture\}',text))==3
 assert names[8]=='terminal_convergence',names
 fonts={}
 for name in names:
  for ext in ['pdf','svg','png']:assert (ROOT/'figures'/f'{name}.{ext}').is_file()
  doc=fitz.open(ROOT/'figures'/f'{name}.pdf')
  used=sorted({f[3].split('+')[-1] for p in doc for f in p.get_fonts()})
  assert used and all('TimesNewRoman' in f for f in used),(name,used)
  assert all(not re.search('[\u4e00-\u9fff]',p.get_text()) for p in doc)
  sizes={round(sp['size'],2) for page in doc for b in page.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for sp in line['spans']}
  assert 11.0 in sizes and all(any(abs(x-v)<.03 for v in [11,7.7,5.39]) for x in sizes),(name,sizes)
  fonts[name]=used
 for p in (ROOT/'generated').glob('*.tex'):assert not re.search('[\u4e00-\u9fff]',p.read_text(encoding='utf-8'))
 fg=(ROOT/'make_figures.py').read_text(encoding='utf-8')
 assert 'tr.height_m<=20000' in fg and 'np.linspace(0,20000,601)' in fg and 'xlim=(0,20)' in fg
 import pandas as pd
 profiles=[pd.read_csv(ROOT/'results'/f'trajectory_{n}km.csv') for n in [40,50,60,70,80]]
 for p in profiles:assert abs(p.height_m.min()-2566.8)<1e-5
 speeds20=[float(np.interp(20000,p.height_m.iloc[::-1],p.speed_mps.iloc[::-1])) for p in profiles]
 speeds10=[float(np.interp(10000,p.height_m.iloc[::-1],p.speed_mps.iloc[::-1])) for p in profiles]
 assert np.ptp(speeds20)>90 and np.ptp(speeds10)<.1
 assert check_atmosphere(text)<.002
 r=json.loads((ROOT/'results.json').read_text(encoding='utf-8'))
 assert r['calibration']['training_flights']==['2012-03-15','2012-07-25','2012-10-14']
 assert set(r['october_curve_errors'])=={'n','rmse_mps'}
 expected=[105.58,56.36,66.28]
 actual=[r['free_g_extrapolation']['maximum_height_m']/1000,r['maximum_height_m']/1000,r['spin_response']['maximum_height_m']/1000]
 assert [round(v,2) for v in actual]==expected
 assert r['limits']['recovery_k']==390 and r['limits']['spin_g']==5
 frozen=json.loads((ROOT/'docs/numerical_baseline.json').read_text(encoding='utf-8'))
 for path,digest in frozen['source_hashes'].items():assert sha(ROOT/path)==digest,path
 for path,digest in frozen['input_hashes'].items():assert sha(ROOT/path)==digest,path
 assert r['calibration']['area_m2']==frozen['area_m2']
 assert max(abs(a-b) for a,b in zip(actual,frozen['heights_km']))<1e-8
 core=ROOT/'appendices/solve.py'
 imports=sorted(n.names[0].name if isinstance(n,ast.Import) else n.module for n in ast.walk(ast.parse(core.read_text(encoding='utf-8'))) if isinstance(n,(ast.Import,ast.ImportFrom)))
 assert imports==['bisect','math']
 # -I ignores PYTHONPATH and -S disables site-packages: no external numerical solver can be imported.
 execution=subprocess.run([sys.executable,'-I','-S',str(core)],cwd=ROOT,capture_output=True,text=True,check=True)
 values=[float(x) for x in execution.stdout.splitlines()]
 assert len(values)==3 and max(abs(a-b) for a,b in zip(values,actual))<.003
 assert not execution.stderr
 binding=json.loads((ROOT/'main.build.json').read_text(encoding='utf-8'))
 assert sha(ROOT/'main.pdf')==binding['main_pdf_sha256']
 for item in binding['resources']:assert sha(ROOT/item['source'])==item['sha256']
 pdf=fitz.open(ROOT/'main.pdf');joined='\n'.join(p.get_text() for p in pdf)
 assert all(name in joined for name in ('骆欣畅','魏贤东','刘奕韬'))
 assert '\ufffd' not in joined
 assert '??' not in joined
 first_page=pdf[0].get_text()
 assert 'Team number:' not in first_page and 'Problem A' in first_page
 assert all(name in first_page for name in ('骆欣畅','魏贤东','刘奕韬'))
 assert 'Summary' in first_page and 'Keywords:' not in first_page
 assert 'Contents' not in first_page
 assert pdf.metadata.get('author')=='Luo Xinchang, Wei Xiandong, Liu Yitao'
 abstract=first_page.split('Summary',1)[1]
 abstract_words=re.findall(r"[A-Za-z]+(?:[-'][A-Za-z]+)*",abstract)
 assert len(abstract_words)<=300
 abstract_space_tokens=len(abstract.split())
 assert abstract_space_tokens<=300
 assert abstract.count('The total mass of a space diver')==1
 assert text.count(r'\begin{abstract}')==1 and text.count(r'\end{abstract}')==1
 assert 'UPC 2023' not in joined
 assert 'Problem A:' in joined and 'Page 2 of' in pdf[1].get_text()
 assert 'Team number:' not in pdf[1].get_text()
 toc=[sp for b in pdf[1].get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for sp in l['spans'] if 'Introduction' in sp['text']]
 assert toc and all(abs(sp['size']-11.96)<.02 for sp in toc)
 spans=[sp for p in pdf for b in p.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for sp in line['spans']]
 appendix_pages=[p for p in pdf if 'def atmosphere(h, load_extension=False)' in p.get_text()]
 assert len(appendix_pages)==1
 code_sizes={round(sp['size'],2) for b in appendix_pages[0].get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for sp in line['spans'] if sp['font'].startswith('LMMono')}
 body_sizes={round(sp['size'],2) for sp in spans if sp['font']=='TimesNewRomanPSMT' and len(sp['text'])>30}
 assert code_sizes=={11.96},code_sizes
 assert 11.96 in body_sizes,body_sizes
 # All text stays inside the physical page; the build separately rejects overfull boxes.
 assert all(sp['bbox'][0]>=0 and sp['bbox'][2]<=pdf[0].rect.width+.1 for sp in spans)
 log=(ROOT/'build/latex/main.log').read_text(encoding='utf-8',errors='replace')
 assert not any(x in log for x in ['Warning:','Missing character:','Overfull'])
 report={'status':'PASS','abstract_lexical_words':len(abstract_words),'abstract_whitespace_tokens':abstract_space_tokens,'abstract_policy':'UPC first-page summary of at most 300 words, including numeric tokens','pages':len(pdf),'section_count':8,'data_figures':14,'tikz_diagrams':3,'figure_fonts':fonts,'nominal_body_and_code_size_pt':12,'pdf_body_and_code_size_pt':11.96,'figure12_speed_spread_20km_mps':float(np.ptp(speeds20)),'figure12_speed_spread_10km_mps':float(np.ptp(speeds10)),'core_stdout':execution.stdout,'three_core_results_km':values,'paper_heights_km':actual,'core_minus_research_m':[(a-b)*1000 for a,b in zip(values,actual)],'independent_core_command':'python -I -S appendices/solve.py','references_all_used':True,'equation_and_figure_references_valid':True,'pdf_sha256':sha(ROOT/'main.pdf')}
 (ROOT/'build/revision/acceptance.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 print(json.dumps({k:v for k,v in report.items() if k!='figure_fonts'},indent=2))
if __name__=='__main__':main()
