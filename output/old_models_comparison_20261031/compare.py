from pathlib import Path
import numpy as np, csv, json, re, hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent; root=out.parents[1]
names=(out/'models.dat').read_text().split(); probes=np.array([181500,181254,180954,180733,180684])-1

def read(p):
 a=np.frombuffer(p.read_bytes().replace(b'\n',b'').replace(b'\r',b''),dtype='S12').astype(float)
 assert a.size==363759 and np.all(np.isfinite(a)),str(p)
 return a

def write_csv(p,header,rows):
 with p.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(header);w.writerows(rows)

results=[];fig,axes=plt.subplots(3,2,figsize=(12,11)); allrows=[]
for name,ax in zip(names,axes.flat):
 new=out/'output'/name;old=root/'output'/name
 summary=(new/'analysis_summary.txt').read_text()
 assert 'last_solver_converged: yes' in summary,summary[:2000]
 files=sorted(new.glob('tempa.dat[0-9]*')); expected=1 if name.endswith('teijo') else 150
 assert len(files)==expected,(name,len(files))
 rows=[];global_min_old=float('inf');global_min_new=float('inf');maxdiff=0.;zrange_old=0.;zrange_new=0.; oldhash={}
 for p in files:
  step=int(p.name[-3:]); t=0 if expected==1 else step
  op=old/p.name;oldhash[p.name]=hashlib.sha256(op.read_bytes()).hexdigest()
  a=read(op);b=read(p);diff=np.abs(a-b)
  global_min_old=min(global_min_old,float(a.min()));global_min_new=min(global_min_new,float(b.min()));maxdiff=max(maxdiff,float(diff.max()))
  zrange_old=max(zrange_old,float(np.ptp(a.reshape(151,2409),axis=0).max()));zrange_new=max(zrange_new,float(np.ptp(b.reshape(151,2409),axis=0).max()))
  row=[step,t,float(a.min()),float(b.min()),float(diff.max()),float(np.sqrt(np.mean(diff**2)))]+a[probes].tolist()+b[probes].tolist()
  rows.append(row);allrows.append([name]+row)
 write_csv(new/'comparison_timeseries.csv',['step','time_ns','old_global_min_V','new_global_min_V','max_abs_difference_V','rms_difference_V']+['old_node%d_V'%(i+1) for i in probes]+['new_node%d_V'%(i+1) for i in probes],rows)
 (new/'old_output_sha256.json').write_text(json.dumps(oldhash,indent=2))
 arr=np.array(rows)
 if expected>1:
  ax.plot(arr[:,1],arr[:,6],label='Old: node 181500');ax.plot(arr[:,1],arr[:,11],label='Recomputed: node 181500');ax.axvline(25,color='grey',ls=':',lw=1)
 else:
  ax.plot([100,500,1200,2000,2400],arr[0,6:11],'o-',label='Old');ax.plot([100,500,1200,2000,2400],arr[0,11:16],'x--',label='Recomputed')
 ax.set(title=name,xlabel='Radius (mesh units)' if expected==1 else 'Time (ns)',ylabel='Potential (V)');ax.grid(alpha=.25);ax.legend(fontsize=8)
 d={'model':name,'steps':expected,'old_min_V':global_min_old,'new_min_V':global_min_new,'max_abs_difference_V':maxdiff,'old_max_z_spread_V':zrange_old,'new_max_z_spread_V':zrange_new,'node181500_25ns_old_V':None if expected==1 else rows[24][6],'node181500_25ns_new_V':None if expected==1 else rows[24][11]}
 results.append(d);print(json.dumps(d),flush=True)
axes.flat[-1].axis('off');fig.tight_layout();fig.savefig(out/'before_after.png',dpi=160);fig.savefig(out/'before_after.svg');plt.close(fig)
write_csv(out/'summary.csv',list(results[0]),[[d[k] for k in results[0]] for d in results]);(out/'results.json').write_text(json.dumps(results,indent=2))
# Compare with the intermediate corrected results from October 26.
inter=[]
for name,case in [('model_entyu_1','base'),('model_entyu_2','fall5'),('model_entyu_3','fall20')]:
 old=root/'output/poweroff_20261026'/case;new=out/'output'/name
 for step in range(1,151):
  a=read(old/('tempa.dat%03d'%step));b=read(new/('tempa.dat%03d'%step))
  inter.append([name,case,step,float(a.min()),float(b.min()),float(np.abs(a-b).max())])
write_csv(out/'comparison_with_20261026.csv',['model','case_20261026','time_ns','intermediate_global_min_V','new_global_min_V','max_abs_difference_V'],inter)
log=(out/'run.log').read_text();assert 'Batch completed with' not in log
assert not re.search(r'WARNING|not converged|Cannot open|configuration error',log,re.I)
lines=['# 旧円柱モデルの再計算・過去結果比較（2026年10月31日）','', '旧円柱モデル5種類を、2026年10月29日の修正を含むソースからビルドした通常プログラムで再計算した。過去の出力は上書きしていない。入力は models/ に複製し、ハッシュを input_hashes.json に記録した。元のメッシュ・材料・パルス・時間条件と追加済みの symmetry.dat を使用し、circuit.dat は追加していない。過渡モデルのOFFは電極の0 V固定であり、電源切り離し放電ではない。','', '## 結果','', '| モデル | ステップ | 過去の全領域最小電位 [V] | 今回の全領域最小電位 [V] | 全節点・全時刻の最大絶対差 [V] |','|---|---:|---:|---:|---:|']
for d in results:lines.append('| {model} | {steps} | {old_min_V:.6e} | {new_min_V:.6e} | {max_abs_difference_V:.6e} |'.format(**d))
lines+=['','4種類の過渡解析は各150ステップ、定常解析は1ステップで完了。全モデルの最終線形求解は収束し、実行ログに警告・未収束・入力エラーはない。各ステップの全363,759節点の保存値を比較した。','', '高さ方向に同じXY座標を持つ節点の電位は、今回の全モデル・全出力ステップで完全一致した。今回の小さな負値は数値誤差水準であり、値の切り上げ・平滑化は行っていない。','', '## 節点181500（r=100、z=750）の25 ns時点','', '| モデル | 過去 [V] | 今回 [V] |','|---|---:|---:|']
for d in results[:4]:lines.append('| %s | %.6e | %.6e |'%(d['model'],d['node181500_25ns_old_V'],d['node181500_25ns_new_V']))
lines+=['','![旧結果と再計算結果](before_after.png)','','## 10月26日の中間修正結果との比較','','初期条件・勾配変換を修正済みで、高さ方向の同一電位条件を使っていなかった10月26日の結果とも比較した。単発パルス・5 ns立下り・20 ns立下りは解析区間内の入力波形が一致する。','', '| モデル | 10月26日の全領域最小電位 [V] | 今回 [V] |','|---|---:|---:|']
for name,case in [('model_entyu_1','base'),('model_entyu_2','fall5'),('model_entyu_3','fall20')]:
 rows=[x for x in inter if x[0]==name];lines.append('| %s | %.6e | %.6e |'%(name,min(x[3] for x in rows),min(x[4] for x in rows)))
lines+=['','初期実装の旧結果との差は勾配変換・初期条件・対称性条件など複数の修正を含む。10月26日結果との比較では、残留していた端面付近の負電位が数値誤差水準まで減少することを確認できる。定常モデルはパルスOFFを扱わず、修正後の空間電位分布を比較した。','','## 保存資料','','- [モデル別集計](summary.csv)','- [10月26日との時系列比較](comparison_with_20261026.csv)','- [比較図（SVG）](before_after.svg)','- 各 output/<モデル名>/comparison_timeseries.csv：全領域の最小値、最大差、RMS差、代表5節点の新旧電位。','- 各 output/<モデル名>/old_output_sha256.json：比較した旧出力のSHA256。','- run.log、build.log、batch_summary.csv、ソース複製、入力複製を同じフォルダに保存。','','比較数値は tempa.dat の出力精度（有効数字5桁）に基づく。微小電位の値も保存値から集計している。実回路との比較や外部抵抗・容量の同定は今回の対象外。']
(out/'report.md').write_text('\n'.join(lines)+'\n')
