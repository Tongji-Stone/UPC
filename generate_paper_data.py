"""All numerical paper macros and tables are emitted from the current run."""
from pathlib import Path
import csv,json
import model as m
ROOT=m.ROOT;OUT=ROOT/'generated';OUT.mkdir(exist_ok=True)
r=json.loads((ROOT/'results.json').read_text(encoding='utf-8'));d=r['cases']['limit'];f=r['case_summaries']['limit']['free'];c=r['case_summaries']['limit']['canopy']
def table(name,rows):
 (OUT/(name+'.tex')).write_text('\n'.join(' & '.join(map(str,row))+r' \\' for row in rows)+'\n',encoding='utf-8')
values={'FitArea':(r['calibration']['area_m2'],4),'MaxHeight':(r['maximum_height_m']/1000,2),
'PeakV':(d['speed_mps'],1),'PeakMach':(d['mach'],3),'FreeG':(d['free_g'],3),'SpinG':(d['spin_g'],3),'OpenG':(d['opening_g'],3),
'LandingV':(d['landing_mps'],4),'TotalTime':(d['total_time_s'],1),'FreeTime':(f['duration_s'],1),'ChuteTime':(c['duration_s'],1),
'CanopyArea':(c['canopy_reference_area_m2'],2),'CanopyB':(c['canopy_drag_area_m2'],2),'OpenV':(c['opening_speed_mps'],2),
'CurveRMSE':(r['october_curve_errors']['rmse_mps'],2),'FitLoss':(r['calibration']['relative_sse'],6),
'TerminalGround':(r['terminal_speeds_mps']['0.0'],2),'TerminalOpen':(r['terminal_speeds_mps']['2566.8'],2),'TerminalFive':(r['terminal_speeds_mps']['5000.0'],2),
'FastOpenG':(r['canopy_4s']['opening_g_peak'],3),'FastOverTime':(r['canopy_4s']['opening_exposure']['total_s'],2),
'MinimumInflation':(r['minimum_inflation_s'],2),'FullInflationHeight':(c['height_at_full_inflation_m'],1),
'HeightConvergence':(r['convergence']['height_difference_m'],4),'CylinderInertia':(m.transverse_inertia(190,m.DEFAULT_SPIN),3),
'FreeLimit':(r['limits']['free_g'],0),'HeatLimit':(r['limits']['recovery_k'],0),'SpinLimit':(r['limits']['spin_g'],0),'OpenLimit':(r['limits']['opening_g'],0),'LandingLimit':(r['limits']['landing_mps'],0),
'HeatPeakHeight':(f['recovery_k_peak_height_m']/1000,2),'HeatPeakTime':(f['recovery_k_peak_time_s'],2),'FreePeakHeight':(f['drag_g_peak_height_m']/1000,2),'SpinRPM':((d['spin_g']*m.G0/.5)**.5*60/(2*m.math.pi),2)}
sp=r['spin_response'];spk=sp['boundary']
values.update({'SpinHeight':(sp['maximum_height_m']/1000,2),'SpinOnsetHeight':(spk['onset_height_m']/1000,2),
'SpinOnsetSpeed':(spk['onset_speed_mps'],2),'SpinBoundaryOmega':(spk['omega_response_rad_s'],3),
'SpinReducedRate':(spk['reduced_rate'],4),'SpinResponseTime':(spk['response_s'],0),
'SpinDensityChange':(100*abs(spk['frozen_density_change_fraction']),1),'SpinSpeedChange':(100*abs(spk['frozen_speed_change_fraction']),1)})
low=r['low_altitude_convergence']
ext=r['free_g_extrapolation'];ep=ext['boundary']
values.update({'LoadHeight':(ext['maximum_height_m']/1000,2),'LoadPeakHeight':(ep['peak_height_m']/1000,2),
'LoadPeakTime':(ep['peak_time_s'],2),'LoadAboveFour':(ep['duration_above_4g_s'],2),
'HeatLowLimit':(r['limits']['recovery_k']-20,0),'HeatHighLimit':(r['limits']['recovery_k']+20,0)})
values.update({'OpenSpread':(max(.0001,low['opening_speed_spread_mps']*2),4),'OpenLoadSpread':(max(.0001,low['opening_load_spread_g']*2),4),
'LandSpreadBound':(max(.0001,low['landing_speed_spread_mps']*2),4),
'FamilyOpenMax':(low['opening_load_max_g'],3),
'InflationDistance':(m.CONFIG['opening_height_m']-c['height_at_full_inflation_m'],1)})
for label,name in [('60km','Sixty'),('70km','Seventy'),('80km','Eighty')]:values[name+'HeatTime']=(r['case_summaries'][label]['exposures']['recovery_k']['total_s'],1)
for val,name in [(r['limits']['recovery_k']-20,'HeatLowHeight'),(r['limits']['recovery_k']+20,'HeatHighHeight')]:values[name]=(next(x['max_height_m'] for x in r['threshold_sensitivity'] if x['parameter']=='recovery_k' and x['value']==val)/1000,2)
scan=list(csv.DictReader((ROOT/'results'/'height_sweep.csv').open(encoding='utf-8-sig')))
for key,name in [('free_g','ScanFreeMax'),('spin_g','ScanSpinMax'),('opening_g','ScanOpenMax'),('landing_mps','ScanLandMax')]:values[name]=(max(float(x[key]) for x in scan),3)
four_significant={'FitArea','PeakV','OpenV','TerminalGround','TerminalOpen','TerminalFive',
 'LandingV','ScanLandMax','CylinderInertia','SpinOnsetSpeed','SpinBoundaryOmega'}
scientific={'FitLoss','OpenSpread','OpenLoadSpread','LandSpreadBound'}
def display(k,v,n):
 if k in scientific:
  mantissa,exponent=f'{v:.2e}'.split('e')
  return r'\ensuremath{'+mantissa+r'\times10^{'+str(int(exponent))+'}}'
 if k in four_significant:
  decimals=max(0,3-int(m.math.floor(m.math.log10(abs(v))))) if v else 3
  return f'{v:.{decimals}f}'
 return f'{v:.{n}f}'
values['LandingLag']=(d['landing_mps']-m.CONFIG['landing_design_speed_mps'],2)
scientific.add('LandingLag')
(OUT/'numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+display(k,v,n)+'}' for k,(v,n) in values.items())+'\n',encoding='utf-8')
table('flights',[[x['flight'],f"{x['height_m']/1000:.2f}",f"{x['observed_speed_mps']:.1f}",f"{x['predicted_speed_mps']:.1f}",f"{x['predicted_speed_mps']-x['observed_speed_mps']:+.2f}",f"{x['error_pct']:+.2f}"] for x in r['validation']])
table('cases',[[f"{r['cases'][k]['height_m']/1000:.0f}",f"{r['cases'][k]['speed_mps']:.1f}",f"{r['cases'][k]['mach']:.3f}",f"{r['cases'][k]['free_g']*190*m.G0:.0f}",f"{r['cases'][k]['recovery_k']:.1f}"] for k in ['40km','50km','60km','70km','80km']])
table('peaks',[[label,f"{f[key+'_peak']:.{1 if key in ['speed_mps','recovery_k'] else 3}f}",f"{f[key+'_peak_time_s']:.1f}",f"{f[key+'_peak_height_m']/1000:.2f}"] for key,label in [('speed_mps',r'Speed / (m/s)'),('mach','Mach number'),('drag_g',r'Free-fall load / $g_0$'),('recovery_k','Recovery temperature / K'),('spin_g',r'1 s radial load / $g_0$')]])
labels={'free_g':'Free-fall load','recovery_k':'Recovery temperature','spin_g':'1 s radial load','opening_g':'Opening load','landing_mps':'Landing speed'}
criteria={'free_g':r'$\max n_D\le\FreeLimit$','recovery_k':r'$\max T_{\rm aw}\le\HeatLimit\,\mathrm K$','spin_g':r'$\max n_r\le\SpinLimit$','opening_g':r'$\max n_o\le\OpenLimit$','landing_mps':r'$v_{\rm land}\le\LandingLimit\,\mathrm{m/s}$'}
sources={'free_g':'Spaceflight reference; assumption','recovery_k':'Thermal environment; assumption','spin_g':'1 s response; assumption','opening_g':'Load assumption','landing_mps':'Canopy reference; assumption'}
def display_boundary(x):
 if x['metric']=='free_g':return f"{ext['maximum_height_m']/1000:.2f} (extrapolated)"
 return 'No boundary within 0--86 km' if x['upper_censored'] else ('No feasible interval' if x['maximum_height_m'] is None else f"{x['maximum_height_m']/1000:.2f}")
table('boundaries',[[labels[x['metric']],criteria[x['metric']],sources[x['metric']],display_boundary(x),('Yes' if x['active'] else 'No')] for key in labels for x in r['boundary_table'] if x['metric']==key])
table('sensitivity',[[labels[x['parameter']],f"{x['value']:g}",('No feasible point' if x['max_height_m'] is None else ('No boundary within 0--86 km' if x['upper_censored'] else f"{x['max_height_m']/1000:.2f}"))] for x in r['threshold_sensitivity']])
rows=list(csv.DictReader((ROOT/'results'/'atmosphere_validation.csv').open(encoding='utf-8-sig')))
table('atmosphere_check',[[int(h/1000),f'{min(float(x["density_error_pct"]) for x in rows if float(x["z_m"])==h):.2f}',f'{max(float(x["density_error_pct"]) for x in rows if float(x["z_m"])==h):.2f}'] for h in [10000,20000,25000]])
print('Generated current paper macros and tables')
