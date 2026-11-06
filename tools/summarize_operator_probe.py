from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1];out=root/'output/poweroff_20261029'
rows=[]
for mode,spacing in [('coarse',20),('full',10),('fine',5),('iso',10),('collapsed',0)]:
 p=out/mode;d=np.genfromtxt(str(p/'history.csv'),delimiter=',',names=True);assert len(d)==26
 log=(p/'run.log').read_text();assert 'WARNING' not in log
 vals={}
 for line in (p/'initial_metrics.txt').read_text().splitlines():
  k,v=line.split();vals[k]=float(v)
 row={'case':mode,'axial_spacing_mesh':spacing,'nodes':int(vals['nodes']),'initial_z_difference_V':vals['max_difference_from_bottom_layer'],'initial_conductive_dielectric_difference_V':vals['max_conductive_dielectric_field_difference'],'off25_min_V':float(d['min_V'][24]),'off25_max_V':float(d['max_V'][24]),'warnings':log.count('WARNING')}
 rows.append(row)
with (out/'cause_summary.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
a=np.fromfile(str(out/'full/phi0.bin'),dtype=np.float64);b=np.fromfile(str(out/'full/phi24.bin'),dtype=np.float64);c=np.fromfile(str(out/'full/phi25.bin'),dtype=np.float64)
# Compare the independently stepped, original-order solve to the existing RCM solve.
old=np.frombuffer((root/'output/poweroff_20261026/base/tempa.dat025').read_bytes().replace(b'\n',b''),dtype='S12').astype(float)
check={'max_abs_difference_from_production_rounded_snapshot_V':float(np.max(np.abs(c-old))),'trace_node884':{'initial_dielectric_V':float(a[883]),'ON24_minus_initial_V':float(b[883]-a[883]),'OFF25_V':float(c[883])},'all_modes_26_steps_and_zero_warnings':True}
(out/'verification.json').write_text(json.dumps(check,indent=2))
assert check['max_abs_difference_from_production_rounded_snapshot_V']<5.1e-9
fig,ax=plt.subplots(1,2,figsize=(12,4))
for mode,spacing in [('coarse',20),('full',10),('fine',5)]:
 field=np.fromfile(str(out/mode/'phi0.bin'),dtype=np.float64).reshape(-1,2409)
 ref=np.fromfile(str(out/'collapsed/phi0.bin'),dtype=np.float64)
 z=np.linspace(0,1500,len(field));ax[0].plot(z,(field[:,883]-ref[883])*1e3,label='dz='+str(spacing))
xs=np.array([20,10,5]);ys=np.array([-r['off25_min_V']*1e3 for r in rows[:3]])
ax[1].plot(xs,ys,'o-',label='3D tetrahedral FEM');ax[1].plot([0,20],[0,ys[0]],'--',color='gray',label='Linear reference')
ax[0].set(xlabel='z (mesh units)',ylabel='Initial potential - 2D reference (mV)',title='Artificial axial variation near end faces')
ax[1].set(xlabel='Axial spacing (mesh units)',ylabel='Magnitude of negative OFF peak (mV)',title='Spatial refinement reduces undershoot',xlim=(0,21),ylim=(0,.34))
for axy in ax:axy.grid(alpha=.2);axy.legend(fontsize=8)
fig.tight_layout();fig.savefig(out/'spatial_cause.png',dpi=180);fig.savefig(out/'spatial_cause.svg');plt.close(fig)
manifest={'date':'2026-10-29','files_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'src/electrostatic_analyzer.cpp',root/'src/electrostatic_analyzer.h',root/'tools/poweroff_operator_probe.cpp',root/'tools/fortran_full_residual.f90',root/'tools/compare_fortran_elements.cpp']},'fortran_reference_root':'/home/ohkura/program/femeemt.ver5.4.2+ohkura_2'}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(rows,indent=2));print(json.dumps(check,indent=2))
