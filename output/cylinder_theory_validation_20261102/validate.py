from pathlib import Path
import numpy as np
import csv,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent
src=out.parent/'old_models_comparison_20261031'
def save(name,header,rows):
 with (out/name).open('w',newline='') as f:
  w=csv.writer(f);w.writerow(header);w.writerows(rows)
def read(p):
 v=np.frombuffer(p.read_bytes().replace(b'\n',b'').replace(b'\r',b''),dtype='S12').astype(float)
 assert len(v)==363759 and np.isfinite(v).all()
 return v
p=src/'models/model_entyu/in.dat'
with p.open() as f:
 h=f.readline().split();n,e=map(int,h[:2]); data=np.fromstring(f.read(),sep=' ')
conn=data[:4*e].astype(int).reshape(e,4)-1
xyz=data[4*e:].reshape(n,3);r=np.hypot(xyz[:,0],xyz[:,1]);a,b=25.,2500.
# Independently determine material membership from the saved element blocks.
sin=(src/'models/model_entyu/sin.dat').read_text().splitlines();nb,nm=map(int,sin[0].split()[:2])
mat=np.zeros(e,dtype=int)
for line in sin[1+nb:1+nb+nm]:
 s,t,m=map(int,line.split());mat[s-1:t]=m
inner=np.unique(conn[mat==1]);outer=np.unique(conn[mat==3]);soil=np.unique(conn[mat==2]);free=np.setdiff1d(soil,np.union1d(inner,outer))
assert np.isclose(r[inner].max(),a) and np.isclose(r[outer].min(),b)
alpha=np.clip(np.log(b/np.maximum(r,1e-100))/np.log(b/a),0,1)
steady=read(src/'output/model_entyu_teijo/tempa.dat001')
probes=np.array([181500,181254,180954,180733,180684])-1
assert np.all(np.isin(probes,free))
err=steady[free]-alpha[free]
center=free[np.isclose(xyz[free,2],750)]
save('steady_all_free_nodes.csv',['node','x','y','z','radius','fem_V','theory_V','error_V'],[[i+1,*xyz[i],r[i],steady[i],alpha[i],steady[i]-alpha[i]] for i in free])
rows=[];ts=[];hashes={};fig,axes=plt.subplots(2,1,figsize=(10,9))
axes[0].scatter(r[center],steady[center],s=4,label='FEM, central plane')
rr=np.geomspace(a,b,400);axes[0].plot(rr,np.log(b/rr)/np.log(b/a),color='black',label='Coaxial theory');axes[0].set(xlabel='Radius (mesh units)',ylabel='Potential (V)');axes[0].legend();axes[0].grid(alpha=.3)
for name in (src/'models.dat').read_text().split():
 if name.endswith('teijo'):continue
 maxerr=0.;sq=0.;ab=0.;cnt=0;shape=0.;off=0.;boundary=0.;zspread=0.;y=1.;rcerr=0.;probeerr=np.zeros(5);probemax=np.zeros(5)
 pulse=[l for l in (src/'models'/name/'pulse.dat').read_text().splitlines() if l.strip() and not l.startswith('#')][0].split()
 _,lo,hi,rise,high,fall,low,phase=map(float,pulse)
 assert rise==0 and lo==0 and hi==1
 for step in range(1,151):
  t=step*1e-9;u=(t+phase)%(high+fall+low)
  # Stabilize exact switching instants in decimal ns.
  u=round(u*1e9,8)*1e-9
  g=1. if u<high-1e-18 else (max(0.,1-(u-high)/fall) if fall else 0.)
  path=src/'output'/name/f'tempa.dat{step:03d}';v=read(path);hashes[str(path.relative_to(src))]=hashlib.sha256(path.read_bytes()).hexdigest()
  d=v[free]-g*alpha[free];maxerr=max(maxerr,float(np.abs(d).max()));sq+=float(d@d);ab+=float(np.abs(d).sum());cnt+=len(d)
  shape=max(shape,float(np.abs(v-g*steady).max()));boundary=max(boundary,float(np.abs(v[inner]-g).max()),float(np.abs(v[outer]).max()))
  zspread=max(zspread,float(np.ptp(v.reshape(151,2409),axis=0).max()))
  if g==0:off=max(off,float(np.abs(v).max()))
  y=(y+.1*g)/1.1;rcerr+=float(np.abs(v[probes]-alpha[probes]*y).sum())
  de=np.abs(v[probes]-g*alpha[probes]);probeerr+=de;probemax=np.maximum(probemax,de)
  ts.append([name,step,t*1e9,g,*v[probes],*list(g*alpha[probes]),*list(v[probes]-g*alpha[probes])])
 rows.append([name,ab/cnt,np.sqrt(sq/cnt),maxerr,shape,off,boundary,zspread,float(probeerr.mean()/150),rcerr/750])
 values=np.array([row[1:] for row in ts if row[0]==name]);axes[1].plot(values[:,1],values[:,4],label=name+' FEM r=500')
save('summary.csv',['model','free_MAE_V','free_RMSE_V','free_max_error_V','max_vs_g_times_steady_V','max_OFF_abs_V','max_boundary_error_V','max_z_spread_V','probe_MAE_instant_V','probe_MAE_RC_V'],rows)
save('timeseries.csv',['model','step','time_ns','input_V']+[f'node{i+1}_{s}_V' for s in ('fem','theory','error') for i in probes],ts)
axes[1].set(xlabel='Time (ns)',ylabel='Potential (V)');axes[1].legend();axes[1].grid(alpha=.3);fig.tight_layout();fig.savefig(out/'comparison.png',dpi=160)
metrics={'free_nodes':len(free),'center_nodes':len(center),'a':a,'b':b,'steady_MAE_V':float(np.abs(err).mean()),'steady_RMSE_V':float(np.sqrt(np.mean(err**2))),'steady_max_error_V':float(np.abs(err).max()),'inner_radii':np.unique(np.round(r[inner],6)).tolist(),'outer_min_radius':float(r[outer].min()),'probes':[[int(i+1),float(r[i]),float(alpha[i]),float(steady[i]),float(steady[i]-alpha[i]),float(100*(steady[i]/alpha[i]-1))] for i in probes]}
(out/'metrics.json').write_text(json.dumps(metrics,indent=2));(out/'output_sha256.json').write_text(json.dumps(hashes,indent=2));print(json.dumps(metrics,indent=2));print(rows)
