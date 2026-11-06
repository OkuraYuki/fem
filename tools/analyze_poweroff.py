from pathlib import Path
import csv,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'output/poweroff_20261026'
COLS=['r100_V','r500_V','r1200_V','r2000_V','r2400_V']
CASES=['base','dt05','dt025','delay','rise5','fall5','fall20','high10','high50']
LABELS={'base':'Immediate ON (dt=1 ns)','dt05':'dt=0.5 ns','dt025':'dt=0.25 ns','delay':'Delay=5 ns','rise5':'Rise=5 ns','fall5':'Fall=5 ns','fall20':'Fall=20 ns','high10':'High=10 ns','high50':'High=50 ns'}
def fields(path):
    raw=path.read_bytes().replace(b'\n',b'').replace(b'\r',b'')
    return np.frombuffer(raw,dtype='S12').astype(float)
def table(path):
    return np.genfromtxt(path,delimiter=',',names=True,encoding='utf8')
def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    data={};metrics=[];status=[];global_metrics=[]
    for name in CASES:
        p=OUT/name; log=(p/'run.log').read_text() if (p/'run.log').exists() else ''
        steps=re.findall(r'Step (\d+)/(\d+): t=([\deE.+-]+) s, Vmin=([\deE.+-]+), Vmax=([\deE.+-]+), max_diff=([\deE.+-]+)',log)
        status.append({'case':name,'completed':'Analysis completed.' in log,'last_step':int(steps[-1][0]) if steps else 0,'expected_steps':int(steps[-1][1]) if steps else 0,'warnings':log.count('WARNING')})
        if not (p/'timeseries.csv').exists():continue
        d=data[name]=table(p/'timeseries.csv')
        rel=[float(x) for x in re.findall(r'Solved linear system with .*?rel_residual=([\deE.+-]+)',log)]
        off_indices=np.flatnonzero((d['input_V']==0)&(d['time_ns']> (30 if name in ['delay','rise5'] else 10 if name=='high10' else 50 if name=='high50' else 25)-1e-7))
        off=int(off_indices[0]) if len(off_indices) else -1
        for col in COLS:
            # Sample post-OFF response at fixed offsets; report sign before fitting a decay.
            minimum=int(np.argmin(d[col])); vmax=float(np.max(np.abs(d[col][off:]))) if off>=0 else float('nan')
            tau=float('nan')
            if off>=0:
                mask=(d['time_ns']>=d['time_ns'][off])&(d['time_ns']<=d['time_ns'][off]+50)&(np.abs(d[col])>1e-13)
                vals=d[col][mask]
                if len(vals)>2 and (np.all(vals>0) or np.all(vals<0)):
                    slope=np.polyfit(d['time_ns'][mask],np.log(np.abs(vals)),1)[0]
                    if slope<0:tau=-1/slope
            metrics.append({'case':name,'probe':col,'dt_ns':float(d['time_ns'][1]-d['time_ns'][0]),'off_ns':float(d['time_ns'][off]) if off>=0 else '',
                'V_before_off':float(d[col][off-1]) if off>0 else '', 'V_first_off':float(d[col][off]) if off>=0 else '',
                'V_min':float(d[col][minimum]),'t_min_ns':float(d['time_ns'][minimum]),'post_off_max_abs_V':vmax,'final_V':float(d[col][-1]),
                'tau_fit_ns':tau,'max_logged_solver_relative_residual':max(rel) if rel else ''})
        if off>=0 and steps:
            times=np.array([float(row[2])*1e9 for row in steps]); minima=np.array([float(row[3]) for row in steps]);maxima=np.array([float(row[4]) for row in steps])
            mask=times>=float(d['time_ns'][off])-1e-7
            local=np.flatnonzero(mask);k=local[int(np.argmin(minima[mask]))]
            fit=(times>=times[k])&(times<=times[k]+50)&(minima<0)
            tau=float('nan')
            if np.count_nonzero(fit)>2:
                slope=np.polyfit(times[fit],np.log(-minima[fit]),1)[0]
                if slope<0:tau=-1/slope
            global_metrics.append({'case':name,'off_ns':float(d['time_ns'][off]),'post_off_global_min_V':float(minima[k]),'minimum_time_ns':float(times[k]),'post_off_global_max_V':float(np.max(maxima[mask])),'global_min_tau_fit_ns':tau,'last_global_min_V':float(minima[-1]),'last_global_max_V':float(maxima[-1])})
        write_csv(p/'extrema_timeseries.csv',[{'step':int(s),'time_ns':float(t)*1e9,'global_min_V':float(lo),'global_max_V':float(hi),'max_step_change_V':float(di)} for s,_,t,lo,hi,di in steps])
    write_csv(OUT/'case_status.csv',status);write_csv(OUT/'metrics.csv',metrics);write_csv(OUT/'global_metrics.csv',global_metrics)
    # Existing timeseries is used only after checking its representative values against raw nodal files.
    oldfile=ROOT/'output/model_entyu_1/single_pulse_nodes_timeseries.csv'
    oldraw=np.loadtxt(oldfile,delimiter=',',skiprows=1)
    checks=[]
    for step in [1,24,25,50,75,150]:
        raw=fields(ROOT/f'output/model_entyu_1/tempa.dat{step:03}')
        csv_v=oldraw[step-1,3:8];raw_v=raw[np.array([181500,181254,180954,180733,180684])-1]
        checks.append({'step':step,'max_abs_csv_raw_difference_V':float(np.max(np.abs(csv_v-raw_v)))})
    write_csv(OUT/'old_csv_verification.csv',checks)
    if 'base' in data:
        base=data['base'];rows=[]
        for step in [1,24,25,50,75,100,150]:
            if step>len(base):continue
            for i,col in enumerate(COLS):
                before=float(oldraw[step-1,3+i]);after=float(base[col][step-1])
                rows.append({'time_ns':step,'probe':col,'old_V':before,'fixed_V':after,'change_V':after-before})
        write_csv(OUT/'before_after.csv',rows)
        fig,axes=plt.subplots(1,2,figsize=(12,4))
        for i,col in enumerate(COLS[:3]):
            axes[0].plot(oldraw[:,1],oldraw[:,3+i],label=f'Old {col[:-2]}',ls='--')
            axes[0].plot(base['time_ns'],base[col],label=f'Fixed {col[:-2]}')
        axes[0].axvline(25,color='gray',lw=.8);axes[0].set(xlabel='Time (ns)',ylabel='Potential (V)',title='Before / after correction');axes[0].legend(fontsize=7)
        for col in COLS:
            axes[1].plot(base['time_ns'],base[col]*1e9,label=col[:-2])
        axes[1].set(xlim=(25,100),ylim=(-30,5),xlabel='Time (ns)',ylabel='Potential (nV)',title='Fixed central-plane probes after OFF');axes[1].legend(fontsize=7)
        fig.tight_layout();fig.savefig(OUT/'before_after.png',dpi=180);fig.savefig(OUT/'before_after.svg');plt.close(fig)
        if all(x in data for x in ['dt05','dt025']):
            rows=[];fine=data['dt025']
            for name in ['base','dt05']:
                coarse=data[name];mask=coarse['time_ns']<=min(coarse['time_ns'][-1],fine['time_ns'][-1])
                # Exact matching of common sample times avoids interpolation across square-wave edges.
                indices=np.rint(coarse['time_ns'][mask]/.25).astype(int)-1
                for col in COLS:
                    diff=np.abs(coarse[col][mask]-fine[col][indices]);ts=coarse['time_ns'][mask]
                    k=int(np.argmax(diff));post=ts>=25
                    rows.append({'case':name,'reference':'dt025','probe':col,'window_end_ns':float(ts[-1]),'max_abs_error_V':float(diff[k]),'error_max_time_ns':float(ts[k]),'post_off_max_abs_error_V':float(np.max(diff[post]))})
            write_csv(OUT/'time_convergence.csv',rows)
            fig,axes=plt.subplots(1,2,figsize=(12,4))
            for name in ['base','dt05','dt025']:
                d=data[name];ext=table(OUT/name/'extrema_timeseries.csv');axes[0].plot(ext['time_ns'],ext['global_min_V']*1e3,label=LABELS[name])
                axes[1].plot(d['time_ns'],d['r1200_V']*1e9,label=LABELS[name])
            axes[0].set(xlim=(25,75),xlabel='Time (ns)',ylabel='Global minimum (mV)',title='Time-step comparison')
            axes[1].set(xlim=(25,75),ylim=(-30,5),xlabel='Time (ns)',ylabel='r1200 potential (nV)',title='Central-plane response')
            for ax in axes:ax.legend()
            fig.tight_layout();fig.savefig(OUT/'time_convergence.png',dpi=180);fig.savefig(OUT/'time_convergence.svg');plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(12,4))
        for name in ['base','delay','rise5','fall5','fall20','high10','high50']:
            if name not in data:continue
            d=data[name];ext=table(OUT/name/'extrema_timeseries.csv')
            axes[0].plot(d['time_ns'],d['r100_V'],label=LABELS[name])
            axes[1].plot(ext['time_ns'],ext['global_min_V']*1e3,label=LABELS[name])
        axes[0].set(xlabel='Time (ns)',ylabel='r100 potential (V)',title='Initial state / waveform comparison')
        axes[1].set(xlabel='Time (ns)',ylabel='Global minimum (mV)',title='Residual response over the whole mesh')
        for ax in axes:ax.legend(fontsize=7)
        fig.tight_layout();fig.savefig(OUT/'waveforms.png',dpi=180);fig.savefig(OUT/'waveforms.svg');plt.close(fig)
    print(json.dumps(status,indent=2));print('Metrics and available plots updated.')
if __name__=='__main__':main()
