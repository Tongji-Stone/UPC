"""Rebuild calculation tables and every simulation consumed by the paper."""
from pathlib import Path
from functools import lru_cache
from dataclasses import replace
from collections import defaultdict
import json,csv,time
import numpy as np
from scipy.optimize import brentq
import model as m

ROOT=m.ROOT;OUT=ROOT/'results'
def read_csv(path):return list(csv.DictReader(Path(path).open(encoding='utf-8-sig')))

def atmosphere_check():
    bins=defaultdict(list)
    for r in read_csv(ROOT/'data'/'wyoming_soundings_selected.csv'):
        bins[(r['station'],r['nominal_time_utc'])].append(r)
    out=[];profiles=[]
    for (station,date),rows in sorted(bins.items()):
        lev=defaultdict(list)
        for r in rows:lev[float(r['geopotential_height_m'])].append(r)
        z=np.array(sorted(lev));p=np.array([np.mean([float(r['pressure_pa']) for r in lev[x]]) for x in z])
        t=np.array([np.mean([float(r['temperature_k']) for r in lev[x]]) for x in z])
        for height in [10000,20000,25000]:
            if z[0]<=height<=z[-1]:
                ti=np.interp(height,z,t);pi=np.exp(np.interp(height,z,np.log(p)))
                ts,ps,rs,_=m.atmosphere(m.R_E*height/(m.R_E-height));rho=pi/(m.R_AIR*ti)
                out.append(dict(station=station,date=date,z_m=height,observed_temperature_k=ti,
                                standard_temperature_k=ts,temperature_delta_k=ti-ts,
                                observed_density_kg_m3=rho,standard_density_kg_m3=rs,density_error_pct=100*(rho/rs-1)))
        for zi,pi,ti in zip(z,p,t):
            h=m.R_E*zi/(m.R_E-zi);ts,ps,rs,_=m.atmosphere(h)
            profiles.append(dict(station=station,date=date,z_m=zi,height_m=h,temperature_k=ti,
                                 pressure_pa=pi,density_kg_m3=pi/(m.R_AIR*ti),
                                 density_error_pct=100*(pi/(m.R_AIR*ti)/rs-1)))
    m.write_csv(OUT/'atmosphere_validation.csv',out);m.write_csv(OUT/'sounding_profiles.csv',profiles)
    m.write_csv(ROOT/'data'/'soundings_vs_standard.csv',out)
    return out

def weather_check(area):
    rows=read_csv(OUT/'sounding_profiles.csv');out=[]
    base=m.summarize(m.freefall(38969.4,121.2,area,stop_height=2566.8))
    for station in ['72364','72365']:
        for date in ['2012-10-14T12:00:00Z','2012-10-15T00:00:00Z']:
            rr=[r for r in rows if r['station']==station and r['date']==date]
            z=np.array([float(r['z_m']) for r in rr]);p=np.array([float(r['pressure_pa']) for r in rr]);t=np.array([float(r['temperature_k']) for r in rr])
            def hybrid(h):
                ts,ps,rs,cs=m.atmosphere(h);zz=m.R_E*h/(m.R_E+h)
                if zz<=z[0] or zz>=z[-1]:return ts,ps,rs,cs
                x=np.clip((zz-z[0])/2000,0,1);y=np.clip((z[-1]-zz)/2000,0,1);w=x*x*(3-2*x)*y*y*(3-2*y)
                temp=(1-w)*ts+w*np.interp(zz,z,t)
                pres=np.exp((1-w)*np.log(ps)+w*np.interp(zz,z,np.log(p)))
                return temp,pres,pres/(m.R_AIR*temp),np.sqrt(m.GAMMA*m.R_AIR*temp)
            sol=m.freefall(38969.4,121.2,area,stop_height=2566.8,atmosphere_fn=hybrid)
            sm=m.summarize(sol);label=station+'_'+date[:10]
            m.save_trajectory(OUT/f'weather_{label}.csv',m.trajectory(sol))
            out.append(dict(station=station,date=date,top_z_m=float(z[-1]),speed_mps=sm['speed_mps_peak'],
                            time_s=sm['speed_mps_peak_time_s'],delta_speed_mps=sm['speed_mps_peak']-base['speed_mps_peak']))
    m.write_csv(OUT/'weather_sensitivity.csv',out);m.write_csv(ROOT/'data'/'weather_trajectory_sensitivity.csv',out)
    return out

def feasible_intervals(grid,metric,limit,fun):
    vals=[fun(float(h))[metric]-limit for h in grid]
    roots=[]
    for i in range(len(grid)-1):
        if vals[i]*vals[i+1]<0:
            roots.append(float(brentq(lambda h:fun(h)[metric]-limit,grid[i],grid[i+1],xtol=.05)))
    edges=[float(grid[0])]+roots+[float(grid[-1])];intervals=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        if fun((lo+hi)/2)[metric]<=limit:intervals.append([lo,hi])
    return {'metric':metric,'limit':limit,'crossings_m':roots,'feasible_intervals_m':intervals,
            'upper_censored':bool(intervals and intervals[-1][1]==float(grid[-1]))}

def main():
    start=time.time();OUT.mkdir(exist_ok=True)
    cal=m.calibrate();area=cal['area_m2'];print('Area calibration',area,flush=True)
    validation=[]
    for row in read_csv(ROOT/'data'/'stratos_summary.csv'):
        october=row['flight']=='2012-10-14';law='mach' if october else 'constant'
        sol=m.freefall(float(row['exit_m']),121.2,area,stop_height=2566.8,law=law)
        sm=m.summarize(sol);m.save_trajectory(OUT/('flight_'+row['flight']+'.csv'),m.trajectory(sol))
        validation.append(dict(flight=row['flight'],law=law,height_m=float(row['exit_m']),
            observed_speed_mps=float(row['peak_mps']),predicted_speed_mps=sm['speed_mps_peak'],
            error_pct=100*(sm['speed_mps_peak']/float(row['peak_mps'])-1),observed_low_g_s=float(row['low_g_s']),
            predicted_low_g_s=sm['low_g_s'],predicted_peak_time_s=sm['speed_mps_peak_time_s'],
            predicted_peak_height_m=sm['speed_mps_peak_height_m']))
        if october:octsol=sol;octsummary=sm
    m.write_csv(OUT/'flight_validation.csv',validation)
    obs=read_csv(ROOT/'data'/'stratos_october_digitized.csv');tr=m.trajectory(octsol,step=.02)
    errors=[]
    for r in obs:
        h=float(r['height_m']);v=float(r['speed_mps']);pred=float(np.interp(h,tr['height_m'][::-1],tr['speed_mps'][::-1]))
        errors.append(dict(height_m=h,observed_mps=v,modeled_mps=pred,residual_mps=pred-v,
                           reading_halfwidth_mps=float(r['speed_reading_halfwidth_mps'])))
    m.write_csv(OUT/'october_curve_comparison.csv',errors)
    resid=np.array([r['residual_mps'] for r in errors])
    curve={'n':len(errors),'rmse_mps':float(np.sqrt(np.mean(resid**2))),'mae_mps':float(np.mean(abs(resid))),
           'max_abs_error_mps':float(max(abs(resid))),'bias_mps':float(np.mean(resid))}
    # Plot-derived readings are correlated; sample count is not independent n.
    atmosphere_check();weather=weather_check(area)
    @lru_cache(maxsize=2000)
    def evaluate(h):return m.descent(float(h),area)[3]
    grid=np.r_[100.,np.arange(1000.,86001.,1000.)]
    scan=[]
    for k,h in enumerate(grid):
        scan.append(evaluate(float(h)))
        if k%15==0:print('Height sweep',int(h/1000),'km',flush=True)
    m.write_csv(OUT/'height_sweep.csv',scan)
    limits=m.CONFIG['limits']
    bounds={key:feasible_intervals(grid,key,val,evaluate) for key,val in limits.items()}
    def joint(h):return max(evaluate(float(h))[key]/val-1 for key,val in limits.items())
    def jfun(h):return {'joint':joint(h)}
    joint_bound=feasible_intervals(grid,'joint',0.,jfun)
    if not joint_bound['feasible_intervals_m']:raise RuntimeError('No feasible height; inspect design')
    hmax=joint_bound['feasible_intervals_m'][-1][1]
    if joint_bound['upper_censored']:raise RuntimeError('Feasible domain reaches 86km: high-atmosphere decision required')
    cases={};case_summaries={}
    for label,h in [('40km',40000),('50km',50000),('60km',60000),('70km',70000),('80km',80000),('limit',hmax)]:
        pre,post,fs,metrics=m.descent(h,area)
        tr=m.trajectory(pre);m.save_trajectory(OUT/f'trajectory_{label}.csv',tr);m.save_trajectory(OUT/f'canopy_{label}.csv',post)
        ex={key:m.exceedance(tr['time_s'],tr[key],val) for key,val in [('drag_g',5.),('recovery_k',400.),('spin_g',2.)]}
        cases[label]=metrics;case_summaries[label]={'free':fs,'canopy':post['summary'],'exposures':ex}
    other=m.chute(2566.8,case_summaries['limit']['canopy']['opening_speed_mps'],190,area,inflation_s=4.)
    m.save_trajectory(OUT/'canopy_limit_4s.csv',other)
    # Select minimum inflation time satisfying the explicit 5g design screen.
    vo=case_summaries['limit']['canopy']['opening_speed_mps']
    inflation_needed=brentq(lambda tau:m.chute(2566.8,vo,190,area,inflation_s=tau)['summary']['opening_g_peak']-5,4,12,xtol=1e-4)
    sensitivity=[]
    for key,values in [('free_g',[4,5,6]),('opening_g',[4,5,6]),('recovery_k',[380,400,420]),('spin_g',[1,2,3]),('landing_mps',[5.5,6,6.5])]:
        for val in values:
            changed=dict(limits);changed[key]=val
            jf=lambda h:{'joint':max(evaluate(float(h))[a]/b-1 for a,b in changed.items())}
            b=feasible_intervals(grid,'joint',0.,jf)
            sensitivity.append(dict(parameter=key,value=val,max_height_m=b['feasible_intervals_m'][-1][1] if b['feasible_intervals_m'] else None,
                                    upper_censored=b['upper_censored']))
    m.write_csv(OUT/'threshold_sensitivity.csv',sensitivity)
    # Thermal limit by two calibration flights: equipment is re-sized consistently.
    area_sensitivity=[]
    for label,a in [('March',cal['single_areas_m2']['2012-03-15']),('joint',area),('July',cal['single_areas_m2']['2012-07-25'])]:
        def heat(h):return m.summarize(m.freefall(h,190,a,stop_height=2566.8))['recovery_k_peak']-400
        root=brentq(heat,40000,80000,xtol=.05)
        metrics=m.descent(root,a)[3]
        area_sensitivity.append(dict(case=label,area_m2=a,thermal_height_m=root,**{k:metrics[k] for k in limits}))
    m.write_csv(OUT/'area_sensitivity.csv',area_sensitivity)
    attitude=[]
    for cn in [0.,.002,.01,.03]:
        for damp in [.25,.5,1.]:
            att=replace(m.DEFAULT_ATTITUDE,yaw_bias=cn,damping=damp)
            sol=m.freefall(hmax,190,area,stop_height=2566.8,attitude=att);sm=m.summarize(sol)
            attitude.append(dict(yaw_bias=cn,damping=damp,spin_g=sm['spin_g_peak'],
                                max_rpm=np.sqrt(sm['spin_g_peak']*m.G0/.5)*60/(2*np.pi)))
            if damp==.5:m.save_trajectory(OUT/f'attitude_cn_{cn:g}.csv',m.trajectory(sol))
    m.write_csv(OUT/'attitude_sensitivity.csv',attitude)
    pitch=[]
    for theta in [5,15,30]:
        sol=m.freefall(hmax,190,area,stop_height=2566.8,attitude=replace(m.DEFAULT_ATTITUDE,theta0_deg=theta))
        sm=m.summarize(sol);m.save_trajectory(OUT/f'pitch_{theta}.csv',m.trajectory(sol))
        pitch.append(dict(theta0_deg=theta,speed_mps=sm['speed_mps_peak'],recovery_k=sm['recovery_k_peak'],drag_g=sm['drag_g_peak']))
    m.write_csv(OUT/'pitch_sensitivity.csv',pitch)
    refined=m.descent(hmax,area,max_step=.05,rtol=1e-10)[3];central=cases['limit']
    convergence={key:float(abs(refined[key]/central[key]-1)) for key in limits}
    refined_h=brentq(lambda h:m.summarize(m.freefall(h,190,area,stop_height=2566.8,max_step=.1,rtol=1e-10))['recovery_k_peak']-400,hmax-100,hmax+100,xtol=.01)
    convergence['height_difference_m']=abs(refined_h-hmax)
    assert max(convergence[k] for k in limits)<.001 and convergence['height_difference_m']<50
    layer_rows=[dict(layer=i,z_base_m=m.HB[i],z_top_m=m.HB[i+1],temperature_k=m.TB[i],lapse_k_per_m=m.L[i],pressure_pa=m.PB[i],density_kg_m3=m.RHOB[i]) for i in range(7)]
    m.write_csv(OUT/'atmosphere_layers.csv',layer_rows)
    out={'calibration':cal,'validation':validation,'october':octsummary,'october_curve_errors':curve,
         'limits':limits,'individual_boundaries':bounds,'joint_boundary':joint_bound,'maximum_height_m':hmax,
         'cases':cases,'case_summaries':case_summaries,'canopy_4s':other['summary'],'minimum_inflation_s':inflation_needed,
         'threshold_sensitivity':sensitivity,'area_sensitivity':area_sensitivity,'attitude_sensitivity':attitude,
         'pitch_sensitivity':pitch,'weather_sensitivity':weather,'convergence':convergence,'runtime_s':time.time()-start}
    (ROOT/'results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Maximum height km:',hmax/1000,'runtime',out['runtime_s'],flush=True)
    return out

if __name__=='__main__':main()
