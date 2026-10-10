"""Generate every quantitative paper table and scalar from the current results."""
from pathlib import Path
import json,csv,math
import model as m
ROOT=m.ROOT; OUT=ROOT/'generated';OUT.mkdir(exist_ok=True)
r=json.loads((ROOT/'results.json').read_text(encoding='utf-8'))
f=r['case_summaries']['limit']['free'];c=r['case_summaries']['limit']['canopy'];d=r['cases']['limit']
def table(name,rows):
    (OUT/(name+'.tex')).write_text('\n'.join(' & '.join(map(str,row))+r' \\' for row in rows)+'\n',encoding='utf-8')
def number(x):
    s=f'{x:.6g}'
    if 'e' in s:
        a,b=s.split('e');return '$'+a+r'\times10^{'+str(int(b))+'}$'
    return s
values={'FitArea':(r['calibration']['area_m2'],4),'MaxHeight':(r['maximum_height_m']/1000,2),
 'PeakV':(d['speed_mps'],1),'PeakMach':(d['mach'],3),'FreeG':(d['free_g'],2),'OpenG':(d['opening_g'],2),
 'SpinG':(d['spin_g'],3),'LandingV':(d['landing_mps'],2),'TotalTime':(d['total_time_s'],1),
 'FreeTime':(f['duration_s'],1),'ChuteTime':(c['duration_s'],1),'CanopyArea':(c['canopy_reference_area_m2'],2),
 'CanopyB':(c['canopy_drag_area_m2'],2),'OpenV':(c['opening_speed_mps'],2),'Tzero':(d['stagnation_k'],1),
 'OctV':(r['october']['speed_mps_peak'],1),'OctError':(abs(r['validation'][2]['error_pct']),2),
 'CurveRMSE':(r['october_curve_errors']['rmse_mps'],2),'CurveMAE':(r['october_curve_errors']['mae_mps'],2),
 'CurveBias':(r['october_curve_errors']['bias_mps'],2),'CurveMax':(r['october_curve_errors']['max_abs_error_mps'],2),
 'Energy':(.5*190*d['landing_mps']**2/1000,2),'Stroke':(d['landing_mps']**2/(2*4*m.G0),3),
 'FirstSonic':(f['sonic_times_s'][0],2),'LastSonic':(f['sonic_times_s'][1],2),
 'SuperTime':(f['sonic_times_s'][1]-f['sonic_times_s'][0],2),'HeightConvergence':(r['convergence']['height_difference_m'],4),
 'FastOpenG':(r['canopy_4s']['opening_g_peak'],2),'FastOverTime':(r['canopy_4s']['opening_over_5g']['total_s'],3),
 'MinimumInflation':(r['minimum_inflation_s'],2),'SixtyHeatTime':(r['case_summaries']['60km']['exposures']['recovery_k']['total_s'],2),
 'SeventyHeatTime':(r['case_summaries']['70km']['exposures']['recovery_k']['total_s'],2),
 'EightyHeatTime':(r['case_summaries']['80km']['exposures']['recovery_k']['total_s'],2)}
(OUT/'numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+f'{v:.{n}f}'+'}' for k,(v,n) in values.items())+'\n',encoding='utf-8')
table('layers',[[i+1,f'{m.HB[i]/1000:g}--{m.HB[i+1]/1000:g}',f'{m.L[i]*1000:g}',f'{m.TB[i]:.2f}',number(m.PB[i]),number(m.RHOB[i])] for i in range(7)])
table('flights',[[x['flight'],f"{x['height_m']/1000:.4f}",f"{x['observed_speed_mps']:.2f}",f"{x['predicted_speed_mps']:.2f}",f"{x['error_pct']:+.2f}",f"{x['observed_low_g_s']:.1f}/{x['predicted_low_g_s']:.2f}"] for x in r['validation']])
table('cases',[[f"{x['height_m']/1000:.2f}",f"{x['speed_mps']:.1f}",f"{x['mach']:.3f}",f"{x['free_g']:.3f}",f"{x['recovery_k']:.1f}",f"{x['spin_g']:.3f}"] for x in r['cases'].values()])
table('peaks',[[label,f"{f[key+'_peak']:.3f}",f"{f[key+'_peak_time_s']:.2f}",f"{f[key+'_peak_height_m']/1000:.3f}"] for key,label in [('speed_mps',r'速度 / (m\,s$^{-1}$)'),('mach','马赫数'),('drag_g',r'阻力过载 / $g_0$'),('recovery_k','恢复温度 / K'),('stagnation_k','滞止温度 / K'),('spin_g',r'头部径向载荷 / $g_0$')]])
labels={'free_g':'自由下降过载上限','opening_g':'开伞过载上限','recovery_k':'恢复温度上限','spin_g':'径向载荷上限','landing_mps':'着陆速度上限'}
table('sensitivity',[[labels[k],'/'.join(str(x['value']) for x in r['threshold_sensitivity'] if x['parameter']==k),'/'.join('无可行点' if x['max_height_m'] is None else f"{x['max_height_m']/1000:.2f}" for x in r['threshold_sensitivity'] if x['parameter']==k)] for k in labels])
table('area',[[x['case'],f"{x['area_m2']:.4f}",f"{x['thermal_height_m']/1000:.3f}"] for x in r['area_sensitivity']])
table('atmosphere_check',[[int(h/1000),f'{min(float(x["density_error_pct"]) for x in rows if float(x["z_m"])==h):.2f}',f'{max(float(x["density_error_pct"]) for x in rows if float(x["z_m"])==h):.2f}'] for rows in [list(csv.DictReader((ROOT/'results'/'atmosphere_validation.csv').open(encoding='utf-8-sig')))] for h in [10000,20000,25000]])
print('Generated paper numbers and 7 data tables')

