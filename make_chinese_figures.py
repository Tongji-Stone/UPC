"""Chinese figures generated solely from this run's numerical and observed data."""
from pathlib import Path
import os,json
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'build'/'mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils.setup_style import setup_style
from utils.export_figure import export_figure
from utils.visual_qa import audit_layout
import model as m
setup_style(journal='general',lang='zh',use_sciplots=False)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'figure.constrained_layout.use':False,'figure.autolayout':False,'font.size':10,'axes.labelsize':10,'legend.fontsize':8,'axes.titlesize':10,'lines.linewidth':1.6})
OUT=ROOT/'results';FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
r=json.loads((ROOT/'results.json').read_text(encoding='utf-8'))
a=r['calibration']['area_m2'];reviews={}
colors=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00']
def load(name):return pd.read_csv(OUT/name)
def save(fig,name):
    fig.tight_layout(pad=1.2)
    issues=audit_layout(fig);reviews[name]=issues
    if any(sev=='FAIL' for sev,msg in issues):raise RuntimeError(str(issues))
    export_figure(fig,str(FIG/name),formats=['pdf','svg','png'],dpi=320,grayscale_preview=True)
    # Preserve DPI in monochrome preview, which the upstream helper omits.
    from PIL import Image
    p=FIG/(name+'_grayscale.png');im=Image.open(p);im.save(p,dpi=(320,320));im.close()
    gray=ROOT/'build'/'grayscale';gray.mkdir(exist_ok=True)
    p.replace(gray/p.name)
    plt.close(fig)
def axes(n=1,height=3.5):return plt.subplots(1,n,figsize=(6.5,height),squeeze=False)[0:2]
# raw 1: measured temperature profiles and standard-layer reference.
fig,axs=plt.subplots(1,2,figsize=(6.5,3.8),sharey=True)
h=np.linspace(0,86000,700);atm=np.array([m.atmosphere(x) for x in h]);z=m.R_E*h/(m.R_E+h)
axs[0].plot(atm[:,0],z/1000,'k-',label='标准大气');axs[1].semilogx(atm[:,2],z/1000,'k-')
s=load('sounding_profiles.csv')
for i,station in enumerate([72364,72365]):
    sub=s[(s.station==station)&(s.date=='2012-10-14T12:00:00Z')]
    axs[0].plot(sub.temperature_k,sub.z_m/1000,color=colors[i],ls=['--',':'][i],label=str(station))
    axs[1].semilogx(sub.density_kg_m3,sub.z_m/1000,color=colors[i],ls=['--',':'][i])
for ax in axs:
    for zb in m.HB[1:-1]/1000:ax.axhline(zb,color='.85',lw=.6,zorder=0)
    ax.set_ylim(0,86);ax.grid(alpha=.2)
axs[0].set(xlabel='温度 (K)',ylabel='位势高度 (km)');axs[1].set(xlabel='密度 (kg/m³)')
axs[0].legend();save(fig,'raw_q1_atmosphere')
# raw 2: observations and fitted-only evidence.
f=load('flight_validation.csv');fig,axs=plt.subplots(1,2,figsize=(6.5,3.1))
x=np.arange(3);labels=['3 月','7 月','10 月']
axs[0].scatter(x,f.observed_speed_mps,color='black',label='公开峰速',zorder=3)
axs[0].scatter(x,f.predicted_speed_mps,marker='s',facecolors='none',edgecolors=colors[0],label='模型')
axs[1].bar(x,f.observed_low_g_s,color='.65',width=.45,label='公开低 g 时长')
axs[1].scatter(x,f.predicted_low_g_s,color=colors[1],marker='D',label='模型',zorder=3)
for ax in axs:ax.set_xticks(x,labels);ax.legend();ax.set_ylim(bottom=0)
axs[0].set_ylabel('峰值速度 (m/s)');axs[1].set_ylabel('初始低于 0.1 g 时长 (s)')
save(fig,'raw_q1_flights')
# raw 3: the digitized v-h curve and reading resolution.
obs=pd.read_csv(ROOT/'data'/'stratos_october_digitized.csv');fig,ax=plt.subplots(figsize=(6.5,3.3))
ax.plot(obs.height_m/1000,obs.speed_mps,color='black',lw=1.3,label='2017 年图 1 数字化观测')
sel=obs.iloc[::24];ax.errorbar(sel.height_m/1000,sel.speed_mps,xerr=sel.height_reading_halfwidth_m/1000,yerr=sel.speed_reading_halfwidth_mps,fmt='none',ecolor=colors[1],capsize=2,label='±3 像素读数范围')
ax.set(xlabel='几何高度 (km)',ylabel='速度 (m/s)',xlim=(12,40));ax.legend();ax.grid(alpha=.2)
save(fig,'raw_q1_october')
# process 1: the reference and extension of drag.
fig,ax=plt.subplots(figsize=(6.5,2.8));ma=np.linspace(0,3,1000)
ax.plot(ma,[m.drag_coefficient(v) for v in ma],color=colors[0]);ax.axvspan(1.25,3,color='.93',label='端点常值延拓')
ax.scatter([.6,1.1,1.25],[.6,.7375,.692],s=22,color=colors[1],zorder=4)
ax.set(xlabel='马赫数 M',ylabel='阻力系数',xlim=(0,3),ylim=(.5,.8));ax.grid(alpha=.2);ax.legend()
save(fig,'process_q1_drag')
# process 2: multiple trajectories, local sound/terminal velocities.
fig,axs=plt.subplots(1,2,figsize=(6.5,3.7))
for i,label in enumerate(['40km','50km','60km','70km','80km']):
    tr=load('trajectory_'+label+'.csv');axs[0].plot(tr.speed_mps,tr.height_m/1000,color=colors[i],ls=['-','--','-.',':','-'][i],label=label.replace('km',' km'))
tr=load('trajectory_limit.csv');sel=tr.iloc[::5]
vt=np.array([m.terminal_speed(v,190,a) for v in sel.height_m])
axs[1].plot(tr.speed_mps,tr.height_m/1000,label='实际速度');axs[1].plot(tr.sound_mps,tr.height_m/1000,'--',label='当地声速')
axs[1].plot(vt,sel.height_m/1000,':',label='局部收尾速度')
for ax in axs:ax.set(xlabel='速度 (m/s)',ylabel='几何高度 (km)',xlim=(0,1000));ax.legend(fontsize=7);ax.grid(alpha=.2)
save(fig,'process_q1_trajectory')
# process 3: acceleration, recovery and stagnation temperatures.
fig,axs=plt.subplots(1,2,figsize=(6.5,3.2));tr=load('trajectory_limit.csv')
axs[0].plot(tr.time_s,tr.drag_g,color=colors[0]);axs[0].axhline(5,color='.4',ls='--',label='设计值 5 g')
axs[0].set(xlabel='出舱后时间 (s)',ylabel='阻力过载 D/(mg0)');axs[0].legend()
axs[1].plot(tr.time_s,tr.recovery_k,label='恢复温度');axs[1].plot(tr.time_s,tr.stagnation_k,'--',label='滞止温度')
axs[1].plot(tr.time_s,tr.temperature_k,':',label='环境温度');axs[1].axhline(400,color='.4',lw=.8)
axs[1].set(xlabel='出舱后时间 (s)',ylabel='温度 (K)');axs[1].legend()
for ax in axs:ax.grid(alpha=.2)
save(fig,'process_q1_heating')
# process 4: solved angular dynamics and its imposed sensitivity.
fig,axs=plt.subplots(1,2,figsize=(6.5,3.3))
for i,deg in enumerate([5,15,30]):
    p=load(f'pitch_{deg}.csv');axs[0].plot(p.time_s,p.theta_deg,ls=['-','--',':'][i],label=f'初始 {deg}°')
for i,cn in enumerate([.002,.01,.03]):
    p=load(f'attitude_cn_{cn:g}.csv');axs[1].plot(p.time_s,p.spin_g,ls=['-','--',':'][i],label=f'Cn={cn:g}')
axs[1].axhline(2,color='.4',lw=.8,ls='-.');axs[0].set(xlabel='时间 (s)',ylabel='俯仰角 (°)')
axs[1].set(xlabel='时间 (s)',ylabel='头部径向旋转载荷 (g0)')
for ax in axs:ax.legend();ax.grid(alpha=.2)
save(fig,'process_q1_attitude')
# process 5: opening load and speed.
fig,axs=plt.subplots(1,2,figsize=(6.5,3.2))
for tau,name,ls in [(4,'canopy_limit_4s.csv','--'),(8,'canopy_limit.csv','-')]:
    p=load(name);sub=p[p.time_s<=15]
    axs[0].plot(sub.time_s,sub.drag_g,ls,label=f'{tau} s 充气');axs[1].plot(sub.time_s,sub.speed_mps,ls,label=f'{tau} s 充气')
axs[0].axhline(5,color='.4',lw=.8,ls=':');axs[0].set(xlabel='拉伞后时间 (s)',ylabel='总气动过载 (g0)')
axs[1].set(xlabel='拉伞后时间 (s)',ylabel='下降速度 (m/s)')
for ax in axs:ax.legend();ax.grid(alpha=.2)
save(fig,'process_q1_opening')
# result 1: actual v-h validation plus residuals.
e=load('october_curve_comparison.csv');fig,axs=plt.subplots(1,2,figsize=(6.5,3.4),gridspec_kw={'width_ratios':[1.25,1]})
axs[0].plot(e.observed_mps,e.height_m/1000,'k-',label='图线数字化观测')
axs[0].plot(e.modeled_mps,e.height_m/1000,'--',color=colors[0],label='面积标定后预测')
axs[0].set(xlabel='速度 (m/s)',ylabel='高度 (km)');axs[0].legend()
axs[1].plot(e.height_m/1000,e.residual_mps,color=colors[1]);axs[1].axhline(0,color='.5',lw=.8)
axs[1].set(xlabel='高度 (km)',ylabel='模型减观测 (m/s)')
for ax in axs:ax.grid(alpha=.2)
save(fig,'result_q1_validation')
# result 2: risk heights, all bounds normalized by design values.
sw=load('height_sweep.csv');fig,ax=plt.subplots(figsize=(6.5,3.6))
for i,(key,label) in enumerate([('free_g','自由段过载'),('recovery_k','恢复温度'),('spin_g','旋转载荷'),('opening_g','开伞过载'),('landing_mps','着陆速度')]):
    ax.plot(sw.height_m/1000,sw[key]/r['limits'][key],label=label,color=colors[i],ls=['-','--','-.',':','-'][i])
ax.axhline(1,color='black',lw=.8);ax.axvline(r['maximum_height_m']/1000,color='.4',ls=':')
ax.set(xlabel='出舱高度 (km)',ylabel='指标 / 相应设计值',xlim=(0,86),ylim=(0,1.5));ax.legend(ncol=3,loc='lower left',bbox_to_anchor=(0,1.01));ax.grid(alpha=.2)
save(fig,'result_q1_limits')
# result 3: torque sensitivity, not a confidence interval.
d=load('attitude_sensitivity.csv');matrix=d.pivot(index='yaw_bias',columns='damping',values='spin_g')
fig,ax=plt.subplots(figsize=(6.5,3.2));im=ax.pcolormesh(np.arange(4)-.5,np.arange(5)-.5,matrix.values,cmap='cividis',rasterized=False)
ax.set_xticks(range(3),[str(x) for x in matrix.columns]);ax.set_yticks(range(4),[str(x) for x in matrix.index])
ax.set(xlabel='旋转阻尼系数',ylabel='偏置力矩系数 Cn')
for i in range(4):
    for j in range(3):ax.text(j,i,f'{matrix.iloc[i,j]:.2f}',ha='center',va='center',color='white' if matrix.iloc[i,j]<matrix.values.max()*.5 else 'black')
cb=fig.colorbar(im,ax=ax,label='峰值头部径向载荷 (g0)');cb.solids.set_rasterized(False)
save(fig,'result_q1_spin_sensitivity')
# result 4: measured atmosphere distribution and weather-input response.
fig,axs=plt.subplots(1,2,figsize=(6.5,3.2));av=load('atmosphere_validation.csv');wv=load('weather_sensitivity.csv')
axs[0].boxplot([av[av.z_m==z].density_error_pct for z in [10000,20000,25000]],tick_labels=['10','20','25'],showfliers=True)
axs[0].axhline(0,color='.6',lw=.8);axs[0].set(xlabel='位势高度 (km)',ylabel='实测密度偏差 (%)')
labels=['72364\n14日12Z','72364\n15日00Z','72365\n14日12Z','72365\n15日00Z']
axs[1].bar(np.arange(4),wv.delta_speed_mps,color=colors[0]);axs[1].set_xticks(np.arange(4),labels,fontsize=7)
axs[1].set_ylabel('峰速相对标准大气变化 (m/s)');axs[1].axhline(0,color='.6',lw=.8)
save(fig,'result_q1_weather')
(ROOT/'build'/'revision'/'figure_layout.json').write_text(json.dumps(reviews,ensure_ascii=False,indent=2),encoding='utf-8')
print('Generated',len(reviews),'logical figures')


