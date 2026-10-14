"""Explicit-area, Mach-dependent space diving model. SI units throughout."""
from __future__ import annotations
from spin_model import SpinResponse,FrozenSpin,transverse_inertia,local_response
from pathlib import Path
import math,csv,json,argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq,minimize_scalar
ROOT=Path(__file__).resolve().parent
CONFIG=json.loads((ROOT/'config.json').read_text(encoding='utf-8'))
R_E=6356766.;G0=9.80665;R_AIR=287.05287;GAMMA=1.4
CP=GAMMA*R_AIR/(GAMMA-1);RECOVERY=.71**(1/3)
HB=np.array([0.,11000.,20000.,32000.,47000.,51000.,71000.,84852.])
L=np.array([-.0065,0.,.001,.0028,0.,-.0028,-.002])
TB=np.empty(8);PB=np.empty(8);TB[0]=288.15;PB[0]=101325.
for i in range(7):
    dz=HB[i+1]-HB[i];TB[i+1]=TB[i]+L[i]*dz
    PB[i+1]=PB[i]*(math.exp(-G0*dz/(R_AIR*TB[i])) if L[i]==0 else (TB[i+1]/TB[i])**(-G0/(R_AIR*L[i])))
RHOB=PB/(R_AIR*TB)
def atmosphere(h):
    if h>86000.01 or not math.isfinite(h):raise ValueError('Atmosphere ends at 86 km')
    return _layer_atmosphere(h)

def _layer_atmosphere(h):
    """Seven-layer formula, shared with the explicitly isolated load extrapolation."""
    h=max(float(h),0.);z=R_E*h/(R_E+h)
    i=min(6,int(np.searchsorted(HB,z,side='right')-1));dz=z-HB[i];t=TB[i]+L[i]*dz
    p=PB[i]*(math.exp(-G0*dz/(R_AIR*TB[i])) if L[i]==0 else (t/TB[i])**(-G0/(R_AIR*L[i])))
    return t,p,p/(R_AIR*t),math.sqrt(GAMMA*R_AIR*t)
def gravity(h):return G0*(R_E/(R_E+max(h,0.)))**2
def drag_coefficient(mach,law='mach'):
    if law=='constant':return .6
    if law!='mach':raise ValueError('Unknown drag law')
    x=abs(float(mach))
    if x<=.6:return .6
    if x<=1.1:return .6+.55*(x-.6)**2
    if x<=1.25:return .74-.32*(x-1.1)
    return .692
DEFAULT_SPIN=SpinResponse(**CONFIG['spin_response'])

def freefall(h0,mass,area,*,stop_height=0.,law='mach',spin_response=None,atmosphere_fn=atmosphere,max_step=.5,rtol=1e-8):
    if not (0<=stop_height<h0<=86000 and mass>0 and area>0):raise ValueError('Invalid freefall parameters')
    return _integrate_freefall(h0,mass,area,stop_height=stop_height,law=law,spin_response=spin_response,
                              atmosphere_fn=atmosphere_fn,max_step=max_step,rtol=rtol)

def _integrate_freefall(h0,mass,area,*,stop_height,law,spin_response,atmosphere_fn,max_step,rtol):
    y0=[h0,0.]
    def rhs(t,y):
        h,v=y[:2];_,_,rho,sound=atmosphere_fn(h)
        drag=.5*rho*drag_coefficient(v/sound,law)*area*v*abs(v)
        out=[-v,gravity(h)-drag/mass]
        return out
    def stop(t,y):return y[0]-stop_height
    stop.terminal=True;stop.direction=-1
    sol=solve_ivp(rhs,(0,3000),y0,events=stop,dense_output=True,max_step=max_step,rtol=rtol,atol=rtol*.1)
    if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('Freefall event not reached')
    sol.mass=mass;sol.area=area;sol.law=law;sol.spin_response=spin_response;sol.atmosphere_fn=atmosphere_fn
    return sol

def _load_extrapolation_atmosphere(h):
    """Chapter 6.1 only: continue the final lapse-rate layer up to 110 km."""
    if not math.isfinite(h) or h>110000.01:raise ValueError('Load extrapolation ends at 110 km')
    return _layer_atmosphere(h)

def freefall_load_extrapolated(h0,mass,area,*,stop_height=2566.8,max_step=.5,rtol=1e-8):
    """Only load metrics are exposed; no extrapolated heat, spin or canopy result."""
    if not (0<=stop_height<h0<=110000 and mass>0 and area>0):raise ValueError('Invalid load extrapolation parameters')
    sol=_integrate_freefall(h0,mass,area,stop_height=stop_height,law='mach',spin_response=None,
                           atmosphere_fn=_load_extrapolation_atmosphere,max_step=max_step,rtol=rtol)
    def load(t):
        h,v=sol.sol(t);_,_,rho,sound=_load_extrapolation_atmosphere(h)
        return .5*rho*drag_coefficient(v/sound)*area*v*v/(mass*G0)
    times=np.linspace(0,sol.t[-1],max(501,int(sol.t[-1]/.05)+1))
    values=np.array([load(t) for t in times]);k=int(np.argmax(values))
    opt=minimize_scalar(lambda t:-load(t),bounds=(times[max(0,k-1)],times[min(len(times)-1,k+1)]),method='bounded')
    peak,t=max((float(values[k]),float(times[k])),(-float(opt.fun),float(opt.x)))
    out={'height_m':float(h0),'free_g':peak,'peak_time_s':t,'peak_height_m':float(sol.sol(t)[0]),
         'freefall_duration_s':float(sol.t[-1]),'uses_extrapolation':bool(h0>86000)}
    for threshold in [3.,4.,5.]:out[f'duration_above_{threshold:g}g_s']=exceedance(times,values,threshold)['total_s']
    return out
def trajectory(sol,step=.15):
    t=np.linspace(0,sol.t[-1],max(501,int(sol.t[-1]/step)+1));y=sol.sol(t);h,v=y[:2]
    temp,p,rho,sound=np.array([sol.atmosphere_fn(float(x)) for x in h]).T
    cd=np.array([drag_coefficient(x,sol.law) for x in v/sound])
    ae=np.full(len(t),sol.area)
    nd=.5*rho*cd*ae*v*v/(sol.mass*G0)
    omega=np.zeros(len(t)) if sol.spin_response is None else local_response(rho,v,sol.area,sol.mass,sol.spin_response)
    spin=np.zeros(len(t)) if sol.spin_response is None else sol.spin_response.head_radius_m*omega**2/G0
    return dict(time_s=t,height_m=h,speed_mps=v,temperature_k=temp,pressure_pa=p,density_kg_m3=rho,
        sound_mps=sound,mach=v/sound,cd=cd,area_m2=ae,drag_g=nd,drag_force_n=nd*sol.mass*G0,spin_g=spin,
        recovery_k=temp+RECOVERY*v*v/(2*CP),stagnation_k=temp+v*v/(2*CP),
        spin_response_rad_s=omega,spin_response_rpm=omega*60/(2*np.pi))
def scalar_metrics(sol,t):
    y=sol.sol(t);h,v=y[:2];temp,p,rho,sound=sol.atmosphere_fn(float(h))
    ae=sol.area
    return {'speed_mps':float(v),'mach':float(v/sound),
        'drag_g':.5*rho*drag_coefficient(v/sound,sol.law)*ae*v*v/(sol.mass*G0),
        'recovery_k':temp+RECOVERY*v*v/(2*CP),'stagnation_k':temp+v*v/(2*CP),
        'spin_g':0. if sol.spin_response is None else float(sol.spin_response.head_radius_m*local_response(rho,v,sol.area,sol.mass,sol.spin_response)**2/G0)}
def summarize(sol):
    tr=trajectory(sol,step=.25);out={}
    for key in ['speed_mps','mach','drag_g','recovery_k','stagnation_k','spin_g']:
        values=tr[key];times=tr['time_s']
        candidates=np.flatnonzero((values[1:-1]>=values[:-2])&(values[1:-1]>=values[2:]))+1
        if np.ptp(values)==0:candidates=[]
        best=[(float(values[0]),float(times[0])),(float(values[-1]),float(times[-1]))]
        for k in candidates:
            opt=minimize_scalar(lambda t:-scalar_metrics(sol,t)[key],bounds=(times[k-1],times[k+1]),method='bounded')
            best.append((-float(opt.fun),float(opt.x)))
        value,t=max(best)
        out[key+'_peak']=value;out[key+'_peak_time_s']=t;out[key+'_peak_height_m']=float(sol.sol(t)[0])
    sonic=[]
    for i in np.flatnonzero(np.diff(np.sign(tr['mach']-1))):
        sonic.append(float(brentq(lambda t:scalar_metrics(sol,t)['mach']-1,tr['time_s'][i],tr['time_s'][i+1])))
    out['sonic_times_s']=sonic
    out['duration_s']=float(sol.t[-1]);return out
def spin_peak_details(sol,sample_step=.25):
    """Most dangerous onset location along freefall, NOT accumulated spin."""
    if sol.spin_response is None:raise ValueError('Spin response parameters required')
    tr=trajectory(sol,step=sample_step);values=tr['spin_g'];times=tr['time_s']
    ks=np.flatnonzero((values[1:-1]>=values[:-2])&(values[1:-1]>=values[2:]))+1
    best=[(float(values[0]),float(times[0])),(float(values[-1]),float(times[-1]))]
    if np.ptp(values)>0:
        for k in ks:
            opt=minimize_scalar(lambda t:-scalar_metrics(sol,t)['spin_g'],bounds=(times[k-1],times[k+1]),method='bounded')
            best.append((-float(opt.fun),float(opt.x)))
    load,t=max(best);h,v=sol.sol(t);rho=sol.atmosphere_fn(float(h))[2];p=sol.spin_response
    omega=float(local_response(rho,v,sol.area,sol.mass,p))
    next_t=min(t+p.response_s,sol.t[-1]);hn,vn=sol.sol(next_t);rn=sol.atmosphere_fn(float(hn))[2]
    return dict(exit_height_m=float(sol.y[0,0]),onset_time_s=t,onset_height_m=float(h),
        onset_speed_mps=float(v),onset_density_kg_m3=float(rho),response_s=p.response_s,
        omega0_rad_s=0.,omega_response_rad_s=omega,radial_g=load,
        reduced_rate=abs(omega)*p.length_m/(2*v) if v>0 else 0.,
        frozen_density_change_fraction=float(rn/rho-1),frozen_speed_change_fraction=float(vn/v-1) if v>0 else 0.)

def exceedance(t,values,limit):
    total=longest=run=0.
    for i in range(len(t)-1):
        a,b=values[i]-limit,values[i+1]-limit;dt=t[i+1]-t[i]
        if a>0 and b>0:part=dt
        elif a<=0 and b<=0:part=0.
        else:part=dt*(a/(a-b) if a>0 else b/(b-a))
        total+=part
        if a<=0:run=0.
        run+=part;longest=max(longest,run)
        if b<=0:run=0.
    return {'total_s':float(total),'longest_s':float(longest)}
def calibrate():
    rows=list(csv.DictReader((ROOT/'data'/'stratos_summary.csv').read_text(encoding='utf-8').splitlines()))
    def peak(a,row):
        sol=freefall(float(row['exit_m']),CONFIG['reference_mass_kg'],a,stop_height=CONFIG['opening_height_m'],law='mach')
        k=np.argmax(sol.y[1]);lo=sol.t[max(0,k-1)];hi=sol.t[min(len(sol.t)-1,k+1)]
        return -minimize_scalar(lambda t:-sol.sol(t)[1],bounds=(lo,hi),method='bounded').fun
    training=[r for r in rows if r['flight'] in ('2012-03-15','2012-07-25','2012-10-14')]
    if len(training)!=3 or len({r['flight'] for r in training})!=3:raise ValueError('Require all three flights exactly once')
    def loss(a):return sum((peak(a,r)/float(r['peak_mps'])-1)**2 for r in training)
    fit=minimize_scalar(loss,bounds=(.3,3.),method='bounded',options={'xatol':1e-9})
    if not fit.success:raise RuntimeError('Calibration failed')
    singles={r['flight']:float(brentq(lambda a:peak(a,r)-float(r['peak_mps']),.3,3.)) for r in training}
    return {'area_m2':float(fit.x),'single_areas_m2':singles,'relative_sse':float(fit.fun),
        'training_flights':[r['flight'] for r in training],'drag_law':'mach','reference_mass_kg':CONFIG['reference_mass_kg']}
def terminal_speed(h,mass,area,law='mach'):
    _,_,rho,sound=atmosphere(h)
    fun=lambda v:.5*rho*drag_coefficient(v/sound,law)*area*v*v-mass*gravity(h)
    speed=brentq(fun,0.,max(50000.,math.sqrt(2*mass*gravity(h)/(.6*rho*area))*2))
    # Rounded source coefficients have a jump at M=1.1. Leave a plot gap if
    # balance lies in the jump; do not present a discontinuity as an exact root.
    return speed if abs(fun(speed))/(mass*gravity(h))<1e-8 else math.nan
def canopy_area(mass,area,target=5.5):return 2*mass*G0/(atmosphere(0)[2]*target**2)-.6*area
def chute(h_open,v_open,mass,area,*,inflation_s=8.,bc=None,max_step=.25,rtol=1e-8):
    if inflation_s<=0 or mass<=0 or area<=0 or h_open<=0 or v_open<0:raise ValueError('Invalid canopy parameters')
    bc=canopy_area(mass,area,CONFIG['landing_design_speed_mps']) if bc is None else bc
    if bc<=0:raise ValueError('Canopy drag area must be positive')
    def values(t,h,v):
        _,_,rho,sound=atmosphere(h);x=min(1.,max(0.,t/inflation_s));s=x*x*(3-2*x)
        return .5*rho*drag_coefficient(v/sound)*area*v*v,.5*rho*bc*s*v*v
    def rhs(t,y):return [-y[1],gravity(y[0])-sum(values(t,*y))/mass]
    def ground(t,y):return y[0]
    ground.terminal=True;ground.direction=-1
    sol=solve_ivp(rhs,(0,3000),(h_open,v_open),events=ground,dense_output=True,rtol=rtol,atol=rtol*.1,max_step=max_step)
    if not sol.success or not len(sol.t_events[0]):raise RuntimeError('Canopy ground event failed')
    t=np.unique(np.r_[np.linspace(0,min(30,sol.t[-1]),1001),np.linspace(0,sol.t[-1],1001)])
    h,v=sol.sol(t);forces=np.array([values(ti,hi,vi) for ti,hi,vi in zip(t,h,v)]);ng=forces.sum(axis=1)/(mass*G0)
    k=int(np.argmax(ng));lo=t[max(0,k-1)];hi=t[min(len(t)-1,k+1)]
    opt=minimize_scalar(lambda s:-sum(values(s,*sol.sol(s)))/(mass*G0),bounds=(lo,hi),method='bounded')
    return {'time_s':t,'height_m':h,'speed_mps':v,'drag_g':ng,'drag_force_n':forces.sum(axis=1),'body_force_n':forces[:,0],'canopy_force_n':forces[:,1],
        'summary':{'opening_speed_mps':v_open,'opening_g_peak':max(float(ng[k]),-float(opt.fun)),
            'peak_time_after_open_s':float(opt.x),'landing_mps':float(v[-1]),'duration_s':float(sol.t[-1]),
            'canopy_drag_area_m2':bc,'canopy_reference_area_m2':bc/CONFIG['canopy_cd'],
            'opening_exposure':exceedance(t,ng,CONFIG['limits']['opening_g']),'inflation_s':inflation_s,
            'inflation_completed':bool(sol.t[-1]>=inflation_s),'height_at_full_inflation_m':float(sol.sol(inflation_s)[0]) if sol.t[-1]>=inflation_s else None}}
def descent(h0,area,*,mass=190.,spin_response=DEFAULT_SPIN,inflation_s=8.,max_step=.5,rtol=1e-8):
    if not 0<h0<=86000:raise ValueError('Exit height must be in (0,86000]')
    hop=min(h0,CONFIG['opening_height_m'])
    if h0>hop:
        pre=freefall(h0,mass,area,stop_height=hop,spin_response=spin_response,max_step=max_step,rtol=rtol)
        fs=summarize(pre);vo=float(pre.y[1,-1])
    else:
        pre=None;vo=0.;fs={k+'_peak':(atmosphere(h0)[0] if k.endswith('_k') else 0.) for k in ['speed_mps','mach','drag_g','recovery_k','stagnation_k','spin_g']};fs['duration_s']=0.
    post=chute(hop,vo,mass,area,inflation_s=inflation_s,max_step=min(.25,max_step),rtol=rtol)
    metrics={'height_m':h0,'free_g':fs['drag_g_peak'],'recovery_k':fs['recovery_k_peak'],'stagnation_k':fs['stagnation_k_peak'],
        'spin_g':fs['spin_g_peak'],'opening_g':post['summary']['opening_g_peak'],'landing_mps':post['summary']['landing_mps'],
        'speed_mps':fs['speed_mps_peak'],'mach':fs['mach_peak'],'total_time_s':fs['duration_s']+post['summary']['duration_s']}
    return pre,post,fs,metrics
def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
def save_trajectory(path,tr):
    keys=[k for k,v in tr.items() if isinstance(v,np.ndarray)]
    write_csv(path,[{k:float(tr[k][i]) for k in keys} for i in range(len(tr[keys[0]]))])
def minimum_run():
    cal=calibrate();area=cal['area_m2'];oct_sol=freefall(38969.4,121.2,area,stop_height=2566.8)
    pre,post,fs,metrics=descent(60000,area)
    out={'calibration':cal,'october':summarize(oct_sol),'case_60km':metrics,'canopy':post['summary']}
    (ROOT/'results').mkdir(exist_ok=True);(ROOT/'results'/'minimum_run.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2));return out
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--minimum',action='store_true');parser.parse_args();minimum_run()
