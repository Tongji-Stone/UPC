"""Compare rounded density expressions printed in the paper with SI outputs."""
import re
import numpy as np
import model as m

def check(text):
    errors=[]
    for quantity,index in [(r'\rho',2)]:
        block=text.split(quantity+r'[\text{kg/m}^3](z[\text{km}])=\begin{cases}')[1].split(r'\end{cases}')[0]
        rows=[row.split('&')[0].strip().rstrip(',') for row in block.split(r'\\') if '&' in row]
        assert len(rows)==7
        for i,expression in enumerate(rows):
            expression=re.sub(r'^\[3pt\]\s*','',expression)
            expression=expression.replace(r'\left','').replace(r'\right','').replace(r'\,','*')
            expression=expression.replace(r'\times','*')
            expression=re.sub(r'\\dfrac\{([^{}]+)\}\{([^{}]+)\}',r'((\1)/(\2))',expression)
            expression=re.sub(r'e\^\{([^{}]+)\}',r'exp(\1)',expression)
            expression=re.sub(r'\^\{([^{}]+)\}',r'**(\1)',expression)
            expression=re.sub(r'(?<=[0-9])(?=[z(])','*',expression)
            expression=expression.replace('}(', '}*(')
            expression=re.sub(r'(?<=\))(?=\()', '*',expression)
            for z in np.linspace(m.HB[i]/1000+.001,m.HB[i+1]/1000-.001,11):
                displayed=eval(expression,{'__builtins__':{},'exp':np.exp},{'z':z})
                h=m.R_E*(z*1000)/(m.R_E-z*1000)
                errors.append(abs(displayed/m.atmosphere(h)[index]-1))
    assert max(errors)<2e-3, max(errors)
    return float(max(errors))

if __name__=='__main__':print(check((m.ROOT/'main.tex').read_text(encoding='utf-8')))
