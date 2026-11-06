from pathlib import Path
import numpy as np,csv,json,re,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1];out=root/'output/poweroff_implementation_20261029'
probes=[181500,181254,180954,180733,180684]
names=['model_entyu_open','model_entyu_fall_sym','model_entyu_clamp_sym'];data={};results=[]
def sample(path,node):
 i=node-1
 with path.open('rb') as f:
  f.seek((i//6)*73+(i%6)*12)
  return float(f.read(12))
for name in names:
 folder=root/'output'/name
 arr=np.array([[i]+[sample(folder/('tempa.dat%03d'%i),n) for n in probes] for i in range(1,151)])
 data[name]=arr
 with (folder/'probe_timeseries.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['time_ns','r100_V','r500_V','r1200_V','r2000_V','r2400_V']);w.writerows(arr)
 log=(out/(name+'.log')).read_text();mins=[float(v) for v in re.findall(r'Vmin=([\deE.+-]+)',log)]
 assert len(mins)==150 and 'WARNING' not in log and min(mins)>-1e-9
 results.append({'case':name,'steps':150,'warnings':0,'minimum_over_steps_V':min(mins),'r100_25ns_V':float(arr[24,1]),'r100_26ns_V':float(arr[25,1]),'r100_50ns_V':float(arr[49,1]),'r100_150ns_V':float(arr[-1,1])})
 for step in [24,25,26,45,150]:
  p=folder/('tempa.dat%03d'%step);v=np.frombuffer(p.read_bytes().replace(b'\n',b''),dtype='S12').astype(float)
  assert len(v)==363759 and np.max(np.ptp(v.reshape(151,2409),axis=0))==0
c=np.genfromtxt(str(root/'output/model_entyu_open/circuit_history.csv'),delimiter=',',names=True)
i=np.arange(len(c));be=np.ones(len(c));be[i>25]=1.1**(-(i[i>25]-25));error=float(np.max(np.abs(c['electrode_V']-be)))
assert error<1e-8 and abs(c['electrode_V'][25]-1)<1e-9 and np.all(np.diff(c['electrode_V'][25:])<=1e-10)
C=float(c['electrode_charge_C'][0]);G=C/10e-9
fig,axes=plt.subplots(1,2,figsize=(12,4))
labels={'model_entyu_open':'Disconnect source; discharge through medium','model_entyu_fall_sym':'Source voltage falls over 20 ns','model_entyu_clamp_sym':'Electrode instantly clamped to 0 V'}
for name in names:axes[0].plot(data[name][:,0],data[name][:,1],label=labels[name])
axes[0].set(xlabel='Time (ns)',ylabel='r=100 potential (V)',title='Implemented OFF boundary conditions',xlim=(15,100));axes[0].legend(fontsize=7)
t=c['time_s']*1e9;exact=np.exp(-np.maximum(t-25,0)/10)
axes[1].plot(t,c['electrode_V'],label='Coupled FEM + circuit (dt=1 ns)');axes[1].plot(t,exact,'--',label='Continuous RC: tau=10 ns')
axes[1].set(xlabel='Time (ns)',ylabel='Electrode potential (V)',title='Positive discharge after source disconnection',xlim=(20,100));axes[1].legend(fontsize=8)
for ax in axes:ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(out/'implemented_discharge.png',dpi=180);fig.savefig(out/'implemented_discharge.svg');plt.close(fig)
info={'cases':results,'open_BE_max_abs_electrode_voltage_error':error,'open_voltage_at_disconnect_V':float(c['electrode_V'][25]),'open_voltage_at_150ns_V':float(c['electrode_V'][-1]),'effective_capacitance_F':C,'effective_medium_resistance_ohm':1/G,'reference_time_constant_s':10e-9,'raw_values_were_not_clipped':True}
(out/'implemented_results.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2))
report='''# 電源OFFモデルの実装と検証（2026-10-29）

高さ方向の同一電位と、電源・放電回路との結合を通常の解析プログラムに実装した。電源を切り離す場合と、接続を保って電圧を低下させる場合を入力で使い分けられる。

## 結果

- 切り離し放電：25 nsの電極電位は1 Vを保持し、26 nsから正のまま減衰した。全150ステップの領域最小値は0 V（接地電極を含む）。
- 電圧立下り：既存モデルの25〜45 nsの立下りに合わせて内部電位も低下した。残る負値は約2.3e-11 Vの数値誤差水準。
- 瞬時0 V固定：従来のように直ちに0 Vへ変わる比較条件。残る負値は約4.8e-11 Vの数値誤差水準。
- 高さ方向に同じXY位置の出力電位は、確認した各時刻で完全に一致した。3ケースとも363,759節点・2,095,200要素、150ステップを完走し、収束警告0件。

![更新後の放電応答](implemented_discharge.png)

## 切り離し放電の参照条件

外部放電抵抗の値が未提示なので、参照例では外部抵抗・追加容量を置かず、導電性の媒質そのものを通って外側接地電極へ放電する。実測した外部抵抗や出力抵抗の影響は、値を入力してから評価する。

媒質のε/σ_xyから時定数は10 ns。dt=1 nsの後退Euler法の理論値と、電極電位は最大差'''+str(error)+''' Vで一致した。25 nsの開放瞬間には電荷を保持し、次の区間から放電を計算する。

連続時間の指数関数との違いは時間積分誤差である。dt=1 nsでは1ステップごとに1/1.1を掛ける減衰となる。電位を目標曲線へ加工しているわけではない。

長さ換算0.001（Fortranのfact=1000に合わせた参照規約）での実効容量は約'''+str(C)+''' F、媒質の実効抵抗は約'''+str(1/G)+''' Ω。外部抵抗を接続すれば、その並列導電経路と追加容量を含めて時定数が変わる。

## 実装した内容

1. symmetry.dat：同じXY座標の節点を同一電位へまとめ、縮約したFEM行列を解く。材料・境界・押出し形状の条件を検査し、元の全節点順に展開して保存する。
2. circuit.dat：対象電極を一つの未知電位にし、電源抵抗・放電抵抗・追加容量と導電電流／変位電流のKCLを結合する。回路を指定しない場合のpulse.datは従来の固定電位条件。
3. initial_state：無自由電荷の誘電場、または開始前から通電した導電定常場を選択。OFF時に内部電位・電荷をリセットしない。
4. circuit_history.csv：初期状態と各ステップの電極電位・接続先基準電位・電荷・抵抗枝電流を保存。
5. 高導電率電極の自由度をまとめる際の相殺誤差が偽の接地抵抗にならないよう、勾配行列が定数電位を消す性質を保って縮約する。

負電位の切り上げ、波形の平滑化、結果の後処理による指数関数への置換は行っていない。

## 使用方法

```bash
./fem model_entyu_open
./fem model_entyu_fall_sym
./fem model_entyu_clamp_sym
```

それぞれoutput/モデル名/へ出力する。元の解析結果は保持した。元の円柱モデル5種類にも、確認済みのsymmetry.datを追加しているが、元の出力フォルダは再計算していない。

詳しい設定と制約は[使い方](../../electrical_circuit.md)を参照。回路結合は単一電極・一度の切り離し・過渡解析に対応する。切り離し時刻は時間格子上に置く。任意の初期電荷ファイルや途中からの再開は未対応。

## 検証記録

- [円柱モデルの指標](implemented_results.json)
- [既存物理回帰試験](physics_regression.log)
- [回路・電荷保存・切替連続性等の試験](circuit_regression.log)
- [設定ファイルの正常・異常系試験](electrical_cli_checks.json)
- 各モデルのoutputフォルダにanalysis_summary.txt、全節点電位、probe_timeseries.csvを保存。切り離しモデルにはcircuit_history.csvも保存。

実機の抵抗・容量・電源立下り波形への同定はまだ行っていない。参照例の正の減衰を、実機の減衰時間と一致した結果として扱うものではない。
'''
(out/'report.md').write_text(report)
manifest={'source_sha256':{str(p):hashlib.sha256((root/p).read_bytes()).hexdigest() for p in [Path('src/electrostatic_analyzer.cpp'),Path('src/electrostatic_analyzer.h'),Path('src/main.cpp'),Path('tests/circuit_regression.cpp')]},'binary_sha256':hashlib.sha256((out/'fem_final').read_bytes()).hexdigest()}
(out/'build_manifest.json').write_text(json.dumps(manifest,indent=2))
