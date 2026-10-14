"""Generate paper figures directly with Matplotlib from current calculated CSVs."""
from pathlib import Path
import os,json
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'build'/'mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Arc,Rectangle
import numpy as np
import pandas as pd
import model as m
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],
 'mathtext.fontset':'custom','mathtext.rm':'Times New Roman',
 'mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold',
 'mathtext.sf':'Times New Roman','mathtext.tt':'Times New Roman','mathtext.cal':'Times New Roman:italic','mathtext.fallback':None,
 'font.size':11,'axes.labelsize':11,'axes.titlesize':11,'legend.fontsize':11,
 'xtick.labelsize':11,'ytick.labelsize':11,'lines.linewidth':1.6,
 'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
OUT=ROOT/'results';FIG=ROOT/'figures';FIG.mkdir(exist_ok=True,parents=True)
r=json.loads((ROOT/'results.json').read_text(encoding='utf-8'));area=r['calibration']['area_m2']
colors=['#0072B2','#D55E00','#009E73','#CC79A7','#AA7700'];styles=['-','--','-.',':','-']
names=[]
def load(n):return pd.read_csv(OUT/n)
def save(fig,name):
 fig.tight_layout(pad=1.1,rect=(0,0,1,.85) if fig.legends else (0,0,1,1))
 for ext in ['pdf','svg','png']:fig.savefig(FIG/(name+'.'+ext),dpi=320)
 names.append(name);plt.close(fig)
def grid(ax):ax.grid(alpha=.2);ax.set_axisbelow(True)
h=np.linspace(0,86000,1000);atm=np.array([m.atmosphere(x) for x in h])
for j,name,label in [(0,'atmosphere_temperature','Atmospheric temperature $T$ (K)'),(2,'atmosphere_density',r'Air density $\rho$ (kg/m$^3$)')]:
 fig,ax=plt.subplots(figsize=(150.4/25.4,2.9))
 (ax.semilogy if j==2 else ax.plot)(h/1000,atm[:,j],color=colors[0])
 for z in m.HB[1:-1]:ax.axvline(m.R_E*z/(m.R_E-z)/1000,color='.8',lw=.6)
 ax.set(xlabel='Geometric altitude $h$ (km)',ylabel=label,xlim=(0,86));grid(ax);save(fig,name)
# All force and geometry diagrams are drawn in main.tex with TikZ.
fig,ax=plt.subplots(figsize=(150.4/25.4,2.7));ma=np.linspace(0,3,1201)
ax.plot(ma,[m.drag_coefficient(x) for x in ma],color=colors[0]);ax.axvspan(1.25,3,color='.93',label='$M > 1.25$: constant extrapolation')
ax.axvline(1,color=colors[1],ls=':',label='$M = 1$');ax.set(xlabel='Mach number $M$',ylabel=r'Drag coefficient $C_{\mathrm{D}}$',xlim=(0,3),ylim=(.55,.8));grid(ax);ax.legend();save(fig,'process_q1_drag')
# Observational checks.
fig,axs=plt.subplots(1,2,figsize=(150.4/25.4,3.2));s=load('sounding_profiles.csv')
axs[0].plot(h/1000,atm[:,0],'k-',label='Standard atmosphere');axs[1].semilogy(h/1000,atm[:,2],'k-')
for i,station in enumerate([72364,72365]):
 sub=s[(s.station==station)&(s.date=='2012-10-14T12:00:00Z')]
 axs[0].plot(sub.height_m/1000,sub.temperature_k,styles[i+1],color=colors[i],label=str(station))
 axs[1].plot(sub.height_m/1000,sub.density_kg_m3,styles[i+1],color=colors[i])
for ax in axs:ax.set(xlabel='Geometric altitude $h$ (km)',xlim=(0,40));grid(ax)
axs[0].set_ylabel('Temperature $T$ (K)');axs[1].set_ylabel(r'Density $\rho$ (kg/m$^3$)');fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=3,bbox_to_anchor=(.5,1),columnspacing=.8,handlelength=1.5);save(fig,'raw_q1_atmosphere')
f=load('flight_validation.csv');fig,axs=plt.subplots(1,2,figsize=(150.4/25.4,3));x=np.arange(3)
axs[0].scatter(x,f.observed_speed_mps,color='black',label='Observed peak');axs[0].scatter(x,f.predicted_speed_mps,marker='s',facecolors='none',edgecolors=colors[0],label='Joint fit')
axs[1].bar(x,f.error_pct,color=colors[:3]);axs[1].axhline(0,color='.4',lw=.8)
for ax in axs:ax.set_xticks(x,['March','July','October']);grid(ax)
axs[0].set_ylabel(r'Peak speed $v_{\mathrm{max}}$ (m/s)');axs[1].set_ylabel('Relative error (%)');axs[0].legend();save(fig,'raw_q1_flights')
e=load('october_curve_comparison.csv');fig,ax=plt.subplots(figsize=(150.4/25.4,3.1))
ax.plot(e.height_m/1000,e.observed_mps,'k-',label='Digitized measured profile')
ax.plot(e.height_m/1000,e.modeled_mps,'--',color=colors[0],label='Simulated profile after joint fit')
ax.text(.02,.94,f"RMSE = {r['october_curve_errors']['rmse_mps']:.2f} m/s",transform=ax.transAxes,va='top')
ax.set(xlabel='Geometric altitude $h$ (km)',ylabel='Speed $v$ (m/s)');grid(ax);ax.legend(loc='lower center');save(fig,'result_q1_validation')
# Chapter 5: fixed initial-height family, independent of constraints.
trs=[load('trajectory_'+s+'.csv') for s in ['40km','50km','60km','70km','80km']]
for key,ylabel,name in [('speed_mps','Descent speed $v$ (m/s)','family_velocity'),('drag_force_n','Aerodynamic drag $D$ (N)','family_force'),('recovery_k',r'Recovery temperature $T_{\mathrm{aw}}$ (K)','family_temperature')]:
 fig,ax=plt.subplots(figsize=(150.4/25.4,3.15))
 for i,tr in enumerate(trs):
  ax.plot(tr.height_m/1000,tr[key],color=colors[i],ls=styles[i],label=f'$h_0 = {40+10*i}$ km')
  if key!='recovery_k':
   idx=np.flatnonzero(np.diff(np.sign(tr.mach.to_numpy()-1)))
   for j in idx:ax.plot(tr.height_m.iloc[j]/1000,tr[key].iloc[j],'o',ms=4,color=colors[i])
 if key!='recovery_k':ax.plot([],[],'o',color='.3',ms=4,label='Near $M = 1$')
 ax.set(xlabel='Geometric altitude $h$ (km)',ylabel=ylabel,xlim=(0,80));grid(ax);ax.legend(ncol=3,loc='lower center',bbox_to_anchor=(.5,1.01),columnspacing=1,handlelength=1.8);save(fig,name)
fig,ax=plt.subplots(figsize=(150.4/25.4,3.1))
for i,tr in enumerate(trs):
 sub=tr[tr.height_m<=20000];ax.plot(sub.height_m/1000,sub.speed_mps,ls=styles[i],color=colors[i],label=f'$h_0 = {40+10*i}$ km')
hl=np.linspace(0,20000,601);ax.plot(hl/1000,[m.terminal_speed(x,190,area) for x in hl],'k--',lw=1.3,label=r'Local terminal speed $v_{\mathrm{t}}$')
ax.set(xlabel='Geometric altitude $h$ (km)',ylabel='Speed $v$ (m/s)',xlim=(0,20));grid(ax);ax.legend(ncol=2,loc='lower center',bbox_to_anchor=(.5,1.01),columnspacing=1,handlelength=1.8);save(fig,'terminal_convergence')
# Separate evidence plots for each constraint.
sw=load('height_sweep.csv')
fig,ax=plt.subplots(figsize=(150.4/25.4,3.1));ex=load('free_g_extrapolation_sweep.csv');eb=r['free_g_extrapolation']
ax.plot(sw.height_m/1000,sw.free_g,color=colors[0],label='Standard model: 0–86 km')
ax.axvspan(86,110,color='#EEE3D5',alpha=.7,label='Extrapolated atmosphere')
ax.plot(ex.height_m/1000,ex.free_g,'--',color=colors[0],label='Extrapolation (Sec. 6.1)')
ax.axvline(86,color='.5',ls=':',lw=1)
ax.axhline(r['limits']['free_g'],color=colors[1],ls='--',label='Screening limit: 5')
hx=eb['maximum_height_m']/1000
ax.plot(hx,r['limits']['free_g'],'o',color=colors[1],ms=5)
ax.annotate(f'Boundary {hx:.2f} km',xy=(hx,r['limits']['free_g']),xytext=(42,4.25),arrowprops={'arrowstyle':'->','color':'.3'})
ax.set(xlabel='Release altitude $h_0$ (km)',ylabel=r'Peak free-fall load factor $n_{\mathrm{D}}$',xlim=(0,110),ylim=(0,5.7));grid(ax);ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,columnspacing=.8,handlelength=1.5);save(fig,'risk_free')
for key,ylabel,name in [('recovery_k',r'Peak recovery temperature $T_{\mathrm{aw}}$ (K)','risk_heat'),('spin_g',r'Maximum 1 s radial load factor $n_{\mathrm{r}}$','risk_spin')]:
 fig,ax=plt.subplots(figsize=(150.4/25.4,2.9));ax.plot(sw.height_m/1000,sw[key],color=colors[0],label='Model prediction')
 ax.axhline(r['limits'][key],color=colors[1],ls='--',label='Screening limit')
 bound=r['individual_boundaries'][key]
 for root in bound['crossings_m']:ax.axvline(root/1000,color='.4',ls=':',label=f'Boundary {root/1000:.2f} km')
 if bound['upper_censored']:ax.text(.03,.12,'Upper boundary not reached',transform=ax.transAxes)
 ax.set(xlabel='Release altitude $h_0$ (km)',ylabel=ylabel,xlim=(0,86));grid(ax);ax.legend();save(fig,name)
fig,axs=plt.subplots(1,2,figsize=(150.4/25.4,3.1))
for tau,name,ls in [(8,'canopy_limit.csv','-')]:
 p=load(name);p=p[p.time_s<=15];axs[0].plot(p.time_s,p.drag_g,ls,label=f'{tau} s inflation');axs[1].plot(p.time_s,p.speed_mps,ls,label=f'{tau} s inflation')
axs[0].axhline(r['limits']['opening_g'],color='.4',ls=':');axs[0].set(xlabel=r'Time since deployment $t-t_{\mathrm{o}}$ (s)',ylabel=r'Total aerodynamic load factor $n_{\mathrm{o}}$');axs[1].set(xlabel=r'Time since deployment $t-t_{\mathrm{o}}$ (s)',ylabel='Descent speed $v$ (m/s)')
for ax in axs:ax.legend();grid(ax)
save(fig,'process_q1_opening')
(ROOT/'build'/'revision').mkdir(parents=True,exist_ok=True)
(ROOT/'build'/'revision'/'figure_manifest.json').write_text(json.dumps(names,ensure_ascii=False,indent=2),encoding='utf-8')
print('Generated',len(names),'figures')
