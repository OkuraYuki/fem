from pathlib import Path
import subprocess, shutil, json, hashlib, csv, re, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path(__file__).resolve().parent
R=D.parents[1]
base=R/'models/model_entyu'
cases=[('baseline',1,0,25,0,25),('amplitude_half',.5,0,25,0,25),('amplitude_double',2,0,25,0,25),('width_short',1,0,10,0,40),('width_long',1,0,40,0,10),('period_short',1,0,25,0,10),('period_long',1,0,25,0,75),('edge_5ns',1,5,25,5,15),('edge_10ns',1,10,25,10,5)]
subprocess.run(['g++','-O2','-std=c++17',str(R/'main.cpp'),str(R/'electrostatic_analyzer.cpp'),'-o',str(D/'fem')],check=True)
manifest={'base_model':str(base),'sources':{n:hashlib.sha256((R/n).read_bytes()).hexdigest() for n in ['main.cpp','electrostatic_analyzer.cpp','electrostatic_analyzer.h']},'cases':cases,'dt_s':1e-9,'steps':200}
(D/'manifest.json').write_text(json.dumps(manifest,indent=2))
rows=[]; data={}; probes=[181500,181254,180954,180733]
for name,amp,rise,high,fall,low in cases:
    m=D/'models'/name;m.mkdir(parents=True,exist_ok=True)
    for n in ['in.dat','sin.dat','tmate.dat','symmetry.dat','voltage_difference.dat']:
        shutil.copy2(base/n,m/n)
    (m/'pulse.dat').write_text(f'# material low high rise high_time fall low_time phase\n1 0 {amp} {rise*1e-9:.12g} {high*1e-9:.12g} {fall*1e-9:.12g} {low*1e-9:.12g} 0\n')
    config=(base/'sina.dat.t').read_text().replace(' 150   0',' 200   0')
    (m/'sina.dat.t').write_text(config)
    start=time.time()
    with (D/(name+'.log')).open('w') as f: subprocess.run([str(D/'fem'),name],cwd=D,stdout=f,stderr=subprocess.STDOUT,check=True)
    out=D/'output'/name;log=(D/(name+'.log')).read_text()
    files=sorted(out.glob('tempa.dat*'));assert len(files)==200,(name,len(files))
    values=[]
    for i,p in enumerate(files,1):
        v=[]
        with p.open('rb') as f:
            for node in probes:
                j=node-1;f.seek((j//6)*73+(j%6)*12);v.append(float(f.read(12)))
        values.append([i]+v)
    arr=np.array(values);assert np.isfinite(arr).all()
    data[name]=arr
    with (out/'probe_timeseries.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['time_ns']+[f'node_{n}_V' for n in probes]);w.writerows(arr)
    mins=[float(v) for v in re.findall(r'Vmin=([\deE.+-]+)',log)]
    assert len(mins)==200 and 'WARNING' not in log
    rows.append(dict(case=name,amplitude_V=amp,rise_ns=rise,high_ns=high,fall_ns=fall,low_ns=low,period_ns=rise+high+fall+low,peak_node181500_V=arr[:,1].max(),mean_node181500_V=arr[:,1].mean(),rms_node181500_V=np.sqrt(np.mean(arr[:,1]**2)),min_field_V=min(mins),runtime_s=time.time()-start))
    print(name,'complete',round(time.time()-start,1),'s',flush=True)
with (D/'summary.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
fig,axs=plt.subplots(2,2,figsize=(13,8))
for ax,names,title in zip(axs.flat,[['baseline','amplitude_half','amplitude_double'],['baseline','width_short','width_long'],['baseline','period_short','period_long'],['baseline','edge_5ns','edge_10ns']],['Amplitude','Pulse width (period=50 ns)','Repetition period','Rise / fall time (period=50 ns)']):
    for name in names:ax.plot(data[name][:,0],data[name][:,1],label=name)
    ax.set(xlabel='Time (ns)',ylabel='Potential at node 181500 (V)',title=title);ax.legend(fontsize=8);ax.grid(alpha=.3)
fig.tight_layout();fig.savefig(D/'comparison.png',dpi=180);plt.close(fig)
err=max(float(np.max(np.abs(data[n][:,1:]-data['baseline'][:,1:]*a))) for n,a in [('amplitude_half',.5),('amplitude_double',2)])
text='''# パルス条件比較解析報告書

作成日：2026年11月4日

## 解析対象と目的

作業フォルダの既存資料と入力を確認し、パルス・理論検証に使用されている円柱モデル `models/model_entyu` を選定した。対象の指定がなかったための選定であり、他形状の解析結果ではない。振幅、High保持時間、周期、立上り・立下り時間の影響を、一度に一種類の条件を変更して比較した。

## 共通条件

363,759節点、2,095,200四面体要素。既存の材料、接地境界、電極配置、高さ方向の同一電位拘束を保持。媒質の面内導電率0.01 S/m、軸方向0.008 S/m、誘電率1e-10 F/m。後退Euler法、時間刻み1 ns、200ステップ、解析期間200 ns。初期状態は開始時の電圧に対応した無自由電荷の誘電場であり、矩形波では開始時からHighである。

回路結合は使用せず、Low期間は電極を0 Vへ固定する。電源切り離し放電とは条件が異なる。空間スケールの単位は追加仮定せず、電位[V]と時間[ns]を比較する。絶対電場・電流・エネルギーは評価していない。

## 条件と計算結果

基準は振幅1 V、High 25 ns、Low 25 ns、周期50 nsの矩形波。幅変更では周期50 nsを、周期変更ではHigh 25 nsを保持。緩やかな波形ではHigh保持25 nsを保持し、Low時間を短縮して周期50 nsを維持した。このため傾斜時間比較には電圧時間積分とLow時間の変化も含まれる。パルス幅は本表ではHigh保持時間を意味する。

|条件|振幅 V|立上り ns|High ns|立下り ns|Low ns|周期 ns|観測点最大 V|平均 V|RMS V|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
'''
for r in rows:
    text+='|'+ '|'.join(str(r[k]) if k not in ['peak_node181500_V','mean_node181500_V','rms_node181500_V'] else f'{r[k]:.7f}' for k in ['case','amplitude_V','rise_ns','high_ns','fall_ns','low_ns','period_ns','peak_node181500_V','mean_node181500_V','rms_node181500_V'])+'|\n'
text+='''
平均とRMSはt=1〜200 nsの200サンプルから計算した。初期値は含まない。周期35 nsの条件は観測期間が整数周期でないため、平均値には窓長の影響が含まれる。観測点は節点181500（既存資料のr=100位置）。4節点の履歴を各ケースのCSVへ保存した。

![パルス条件別の電位履歴](comparison.png)

## 考察

'''
b=rows[0]
text+=f'振幅0.5倍・2倍の4観測点で、基準波形を同倍率にした値との差は最大{err:.3e} Vであった。今回の線形な材料・境界条件では、電位応答は振幅に比例する。\n\n'
for pair,label in [(['width_short','width_long'],'High保持時間'),(['period_short','period_long'],'周期'),(['edge_5ns','edge_10ns'],'立上り・立下り時間')]:
    text+=label+'の変更結果：'+ '、'.join(f"{r['case']}の平均は{r['mean_node181500_V']:.6f} V、RMSは{r['rms_node181500_V']:.6f} V" for r in rows if r['case'] in pair)+f"（基準平均{b['mean_node181500_V']:.6f} V、RMS{b['rms_node181500_V']:.6f} V）。\n\n"
text+='''この一様な円柱条件では面内の誘電率／導電率が共通なので、理想電極の電圧変化が内部電位に直接反映される。短パルスで内部が必ず時定数に従ってゆっくり充電される、という回路応答をこの固定電圧条件へ当てはめることはできない。残留放電や電源抵抗による波形の遅れを調べる場合は、既存のcircuit.datを用いた別条件が必要となる。

## 検証と制約

全9条件で200ステップの出力を確認し、観測値の有限性、ログの警告なしを検査した。振幅比例性は4観測点で照合した。これらは実行整合性の検査であり、メッシュ・時間刻みの収束試験や実測との照合ではない。1 ns未満の過渡や絶縁破壊・非線形材料・発熱は本報告の評価対象外。出力の微小負値は加工せずsummary.csvへ領域最小値として記録した。

## 再現と保存先

`python3 output/pulse_study_20261104/run_study.py` で再実行できる。現在のソースを専用バイナリにビルドし、専用フォルダ内で実行する。既存入力・結果・ルートのbatch_summary.csvは変更しない。再実行時は本比較解析の結果を上書きする。

- manifest.json：使用ソースのSHA-256、基準モデル、条件一覧
- models/条件名/：各ケースの実入力
- output/条件名/：全節点200ステップ、analysis_summary.txt、probe_timeseries.csv
- 条件名.log：実行ログ
- summary.csv：条件・指標・実行時間
- comparison.png：比較グラフ
'''
(D/'report.md').write_text(text)
subprocess.run(['python3',str(D/'verify_report.py')],check=True)
print('REPORT COMPLETE',flush=True)
