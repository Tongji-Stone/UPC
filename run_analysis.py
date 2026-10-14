"""Rebuild calculation tables and every simulation consumed by the paper."""
from pathlib import Path
from functools import lru_cache
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
    roots=[float(h) for h,value in zip(grid,vals) if value==0.]
    for i in range(len(grid)-1):
        if vals[i]*vals[i+1]<0:
            roots.append(float(brentq(lambda h:fun(h)[metric]-limit,grid[i],grid[i+1],xtol=.05)))
    roots=sorted(set(roots))
    edges=sorted(set([float(grid[0]),*roots,float(grid[-1])]))
    intervals=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        if fun((lo+hi)/2)[metric]<=limit:
            if intervals and intervals[-1][1]==lo:intervals[-1][1]=hi
            else:intervals.append([lo,hi])
    for root in roots:
        if not any(lo<=root<=hi for lo,hi in intervals):intervals.append([root,root])
    intervals.sort()
    return {'metric':metric,'limit':limit,'crossings_m':roots,'feasible_intervals_m':intervals,
            'upper_censored':bool(intervals and intervals[-1][1]==float(grid[-1]))}

def main():
    start=time.time();OUT.mkdir(exist_ok=True)
    cal=m.calibrate();area=cal['area_m2'];print('Area calibration',area,flush=True)
    validation=[]
    for row in read_csv(ROOT/'data'/'stratos_summary.csv'):
        october=row['flight']=='2012-10-14';law='mach'
        sol=m.freefall(float(row['exit_m']),121.2,area,stop_height=2566.8,law=law)
        sm=m.summarize(sol);m.save_trajectory(OUT/('flight_'+row['flight']+'.csv'),m.trajectory(sol))
        validation.append(dict(flight=row['flight'],law=law,height_m=float(row['exit_m']),
            observed_speed_mps=float(row['peak_mps']),predicted_speed_mps=sm['speed_mps_peak'],
            error_pct=100*(sm['speed_mps_peak']/float(row['peak_mps'])-1),predicted_peak_time_s=sm['speed_mps_peak_time_s'],
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
    curve={'n':len(errors),'rmse_mps':float(np.sqrt(np.mean(resid**2)))}
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
    # Isolated chapter 6.1 calculation: never feed this atmosphere into evaluate/joint.
    load_scan=[m.freefall_load_extrapolated(h,190,area) for h in np.arange(86000.,110001.,1000.)]
    load_fun=lambda h:m.freefall_load_extrapolated(h,190,area)['free_g']-limits['free_g']
    load_root=brentq(load_fun,86000.,110000.,xtol=.05)
    load_peak=m.freefall_load_extrapolated(load_root,190,area)
    refined_load=lambda h:m.freefall_load_extrapolated(h,190,area,max_step=.1,rtol=1e-10)['free_g']-limits['free_g']
    refined_load_root=brentq(refined_load,load_root-100,load_root+100,xtol=.01)
    load_extension={'purpose':'chapter_6_1_freefall_load_only','standard_top_m':86000.,'extrapolation_top_m':110000.,
        'maximum_height_m':float(load_root),'limit_g':limits['free_g'],'boundary':load_peak,
        'below_boundary_g':load_fun(load_root-10)+limits['free_g'],
        'above_boundary_g':load_fun(load_root+10)+limits['free_g'],
        'height_convergence_m':float(abs(refined_load_root-load_root))}
    assert load_extension['below_boundary_g']<limits['free_g']<load_extension['above_boundary_g']
    assert load_extension['height_convergence_m']<1
    m.write_csv(OUT/'free_g_extrapolation_sweep.csv',load_scan)
    m.write_csv(OUT/'free_g_extrapolation_boundary.csv',[load_peak])
    bounds={key:feasible_intervals(grid,key,val,evaluate) for key,val in limits.items()}
    def joint(h):return max(evaluate(float(h))[key]/val-1 for key,val in limits.items())
    def jfun(h):return {'joint':joint(h)}
    joint_bound=feasible_intervals(grid,'joint',0.,jfun)
    if not joint_bound['feasible_intervals_m']:raise RuntimeError('No feasible height; inspect design')
    hmax=joint_bound['feasible_intervals_m'][-1][1]
    # A domain edge is reported as censored, never as a physical maximum.
    cases={};case_summaries={}
    for label,h in [('40km',40000),('50km',50000),('60km',60000),('70km',70000),('80km',80000),('limit',hmax)]:
        pre,post,fs,metrics=m.descent(h,area)
        tr=m.trajectory(pre);m.save_trajectory(OUT/f'trajectory_{label}.csv',tr);m.save_trajectory(OUT/f'canopy_{label}.csv',post)
        ex={key:m.exceedance(tr['time_s'],tr[key],val) for key,val in [('drag_g',limits['free_g']),('recovery_k',limits['recovery_k'])]}
        cases[label]=metrics;case_summaries[label]={'free':fs,'canopy':post['summary'],'exposures':ex}
    other=m.chute(2566.8,case_summaries['limit']['canopy']['opening_speed_mps'],190,area,inflation_s=4.)
    m.save_trajectory(OUT/'canopy_limit_4s.csv',other)
    # Select minimum inflation time satisfying the explicit 5g design screen.
    vo=case_summaries['limit']['canopy']['opening_speed_mps']
    inflation_needed=brentq(lambda tau:m.chute(2566.8,vo,190,area,inflation_s=tau)['summary']['opening_g_peak']-limits['opening_g'],4,12,xtol=1e-4)
    sensitivity=[]
    for key,values in [('free_g',[4,5,6]),('opening_g',[4,5,6]),('recovery_k',[limits['recovery_k']-20,limits['recovery_k'],limits['recovery_k']+20]),('spin_g',[4,5,6]),('landing_mps',[5.5,6,6.5])]:
        for val in values:
            changed=dict(limits);changed[key]=val
            jf=lambda h:{'joint':max(evaluate(float(h))[a]/b-1 for a,b in changed.items())}
            b=feasible_intervals(grid,'joint',0.,jf)
            sensitivity.append(dict(parameter=key,value=val,max_height_m=b['feasible_intervals_m'][-1][1] if b['feasible_intervals_m'] else None,
                                    upper_censored=b['upper_censored']))
    m.write_csv(OUT/'threshold_sensitivity.csv',sensitivity)
    # Thermal sensitivity by three calibration flights: equipment is re-sized consistently.
    area_sensitivity=[]
    for label,a in [('March',cal['single_areas_m2']['2012-03-15']),('joint',area),('July',cal['single_areas_m2']['2012-07-25']),('October',cal['single_areas_m2']['2012-10-14'])]:
        def heat(h):return m.summarize(m.freefall(h,190,a,stop_height=2566.8))['recovery_k_peak']-limits['recovery_k']
        root=brentq(heat,40000,80000,xtol=.05)
        metrics=m.descent(root,a)[3]
        area_sensitivity.append(dict(case=label,area_m2=a,thermal_height_m=root,**{k:metrics[k] for k in limits}))
    m.write_csv(OUT/'area_sensitivity.csv',area_sensitivity)
    spin_root=bounds['spin_g']['feasible_intervals_m'][-1][1]
    if bounds['spin_g']['upper_censored']:raise RuntimeError('Spin boundary not bracketed in standard atmosphere')
    def spin_case(h,step=.5,rtol=1e-8,sample_step=.25):
        sol=m.freefall(h,190,area,stop_height=2566.8,spin_response=m.DEFAULT_SPIN,max_step=step,rtol=rtol)
        return sol,m.spin_peak_details(sol,sample_step)
    spin_sol,spin_peak=spin_case(spin_root)
    refined_spin_root=brentq(lambda h:spin_case(h,.1,1e-10)[1]['radial_g']-limits['spin_g'],spin_root-100,spin_root+100,xtol=.01)
    spin_result={'interpretation':'maximum local response after 1 s from zero spin, not accumulated flight rotation',
        'parameters':vars(m.DEFAULT_SPIN),'maximum_height_m':spin_root,'limit_g':limits['spin_g'],
        'boundary':spin_peak,'below_boundary_g':spin_case(spin_root-10)[1]['radial_g'],
        'above_boundary_g':spin_case(spin_root+10)[1]['radial_g'],
        'height_convergence_m':abs(refined_spin_root-spin_root),
        'peak_sampling_difference_g':abs(spin_case(spin_root,.5,1e-8,.05)[1]['radial_g']-spin_peak['radial_g'])}
    assert spin_result['below_boundary_g']<limits['spin_g']<spin_result['above_boundary_g']
    assert spin_result['height_convergence_m']<1 and spin_result['peak_sampling_difference_g']<1e-6
    m.write_csv(OUT/'spin_boundary.csv',[spin_peak])
    m.save_trajectory(OUT/'spin_envelope_boundary.csv',m.trajectory(spin_sol))
    frozen=m.FrozenSpin(spin_peak['onset_density_kg_m3'],spin_peak['onset_speed_mps'],area,m.DEFAULT_SPIN.length_m,
        m.transverse_inertia(190,m.DEFAULT_SPIN),m.DEFAULT_SPIN.c_nomega,m.DEFAULT_SPIN.c_n)
    tau=np.linspace(0,m.DEFAULT_SPIN.response_s,201);omega=frozen.solution(tau)
    m.save_trajectory(OUT/'spin_response_boundary.csv',{'elapsed_s':tau,'omega_rad_s':omega,'radial_g':m.DEFAULT_SPIN.head_radius_m*omega**2/m.G0})
    refined=m.descent(hmax,area,max_step=.05,rtol=1e-10)[3];central=cases['limit']
    convergence={key:float(abs(refined[key]/central[key]-1)) for key in limits}
    if joint_bound['upper_censored']:
        convergence['height_difference_m']=0.
    else:
        def refined_joint(h):
            mm=m.descent(h,area,max_step=.1,rtol=1e-10)[3]
            return max(mm[key]/val-1 for key,val in limits.items())
        refined_h=brentq(refined_joint,max(100.,hmax-100),min(86000.,hmax+100),xtol=.01)
        convergence['height_difference_m']=abs(refined_h-hmax)
    assert max(convergence[k] for k in limits)<.001 and convergence['height_difference_m']<50
    boundary_rows=[]
    for key,bound in bounds.items():
        hi=bound['feasible_intervals_m'][-1][1] if bound['feasible_intervals_m'] else None
        boundary_rows.append(dict(metric=key,limit=limits[key],maximum_height_m=hi,
            upper_censored=bound['upper_censored'],active=bool(hi is not None and not bound['upper_censored'] and abs(hi-hmax)<1)))
    m.write_csv(OUT/'constraint_boundaries.csv',boundary_rows)
    terminal={str(h):m.terminal_speed(h,190,area) for h in [0.,2566.8,5000.]}
    family_canopy=[case_summaries[k]['canopy'] for k in ['40km','50km','60km','70km','80km']]
    low_altitude={
        'exit_heights_m':[40000,50000,60000,70000,80000],
        'opening_speed_spread_mps':float(np.ptp([x['opening_speed_mps'] for x in family_canopy])),
        'opening_load_spread_g':float(np.ptp([x['opening_g_peak'] for x in family_canopy])),
        'landing_speed_spread_mps':float(np.ptp([x['landing_mps'] for x in family_canopy])),
        'opening_load_max_g':max(x['opening_g_peak'] for x in family_canopy),
        'opening_load_margin_g':limits['opening_g']-max(x['opening_g_peak'] for x in family_canopy),
        'landing_speed_margin_mps':limits['landing_mps']-max(x['landing_mps'] for x in family_canopy)}

    layer_rows=[dict(layer=i,z_base_m=m.HB[i],z_top_m=m.HB[i+1],temperature_k=m.TB[i],lapse_k_per_m=m.L[i],pressure_pa=m.PB[i],density_kg_m3=m.RHOB[i]) for i in range(7)]
    m.write_csv(OUT/'atmosphere_layers.csv',layer_rows)
    out={'calibration':cal,'validation':validation,'october':octsummary,'october_curve_errors':curve,'free_g_extrapolation':load_extension,
         'terminal_speeds_mps':terminal,'low_altitude_convergence':low_altitude,'boundary_table':boundary_rows,'maximum_height_status':'domain_limited' if joint_bound['upper_censored'] else 'finite','limits':limits,'individual_boundaries':bounds,'joint_boundary':joint_bound,'maximum_height_m':hmax,
         'cases':cases,'case_summaries':case_summaries,'canopy_4s':other['summary'],'minimum_inflation_s':inflation_needed,
         'threshold_sensitivity':sensitivity,'area_sensitivity':area_sensitivity,'spin_response':spin_result,'weather_sensitivity':weather,'convergence':convergence,'runtime_s':time.time()-start}
    (ROOT/'results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Maximum height km:',hmax/1000,'runtime',out['runtime_s'],flush=True)
    return out

if __name__=='__main__':main()
