"""0 V固定パルスの保存済みFEMと同軸理論を比較する。"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = Path(__file__).resolve().parent
source = out.parent / 'cylinder_theory_validation_20261102'
plt.rcParams.update({'font.family': 'Noto Sans CJK JP', 'font.size': 12})
with (source/'timeseries.csv').open() as f:
    rows = [r for r in csv.DictReader(f) if r['model'] == 'model_entyu']
with (source/'steady_all_free_nodes.csv').open() as f:
    radial = [r for r in csv.DictReader(f) if float(r['z']) == 750 and float(r['y']) == 0 and float(r['x']) > 0]
radial.sort(key=lambda r: float(r['radius']))
with (out/'steady_centerline.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(radial[0]))
    writer.writeheader()
    writer.writerows(radial)
with (out/'pulse_timeseries.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
t = np.array([float(r['time_ns']) for r in rows])
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
fig.subplots_adjust(left=.08, right=.98, top=.85, bottom=.13, hspace=.43, wspace=.24)
fig.suptitle('円柱モデル：0 V固定パルスの同軸理論と有限要素法（FEM）の比較', fontsize=20, y=.97)
fig.text(.5, .915, '中心電極：1 Vを25 ns保持 → 0 Vを25 ns保持（周期50 ns）／ 外側電極：常時0 V', ha='center')
fig.text(.5, .88, '保存済みの修正後FEM結果を使用・電源切り離しなし・高さ方向は同電位に拘束', ha='center', fontsize=11)
p = axes[0, 0]
r = np.linspace(25, 2500, 2000)
p.plot(r, np.log(2500/r)/np.log(100), color='#1565a6', lw=2, label='同軸円筒の理論解')
p.plot([float(v['radius']) for v in radial], [float(v['fem_V']) for v in radial], 'o--', color='#c36718', ms=4, lw=1.5, label='有限要素法の解析結果')
p.set(title='定常電位：中心電極1 Vでの半径分布', xlabel='中心軸からの半径 [メッシュ座標単位]', ylabel='電位 [V]', xlim=(0,2500), ylim=(0,1.04))
p.text(.35,.63, '中央断面：高さ750、正のx軸上\n代表点ではFEMが理論より約7%高い', transform=p.transAxes, fontsize=11)
p.legend(fontsize=10)
p.grid(alpha=.25)
for p, node, radius in zip([axes[0,1], axes[1,0], axes[1,1]], [181254,180954,180733], [500,1200,2000]):
    theory = np.array([float(r['node%d_theory_V'%node]) for r in rows])
    fem = np.array([float(r['node%d_fem_V'%node]) for r in rows])
    # 値は各計算時刻の境界値。階段表示で矩形波の切替を示す。
    tt = np.r_[0,t]
    p.step(tt,np.r_[theory[0],theory],where='post',color='#1565a6',lw=2.6,label='同軸円筒の理論解')
    p.step(tt,np.r_[fem[0],fem],where='post',color='#c36718',ls='--',lw=1.8,label='有限要素法の解析結果')
    p.plot(t[::5],fem[::5], 'o',color='#c36718',ms=3)
    p.set(title='パルス応答：半径%d・高さ750（節点%d）'%(radius,node), xlabel='時間 [ns]', ylabel='観測点の電位 [V]', xlim=(0,150),ylim=(-.04*max(fem),1.38*max(fem)))
    p.legend(loc='upper right', fontsize=10)
    high = np.array([float(r['input_V']) for r in rows]) == 1
    p.text(.025,.77,'High中：理論 %.6f V ／ FEM %.6f V'%(np.mean(theory[high]), np.mean(fem[high])),transform=p.transAxes,fontsize=10)
    p.grid(alpha=.25)
fig.text(.08,.045,'共通条件：内側境界半径 a = 25、外側境界半径 b = 2500（メッシュ座標単位）。半径方向 σ = 0.01 S/m、ε = 1e-10 F/m。\n'
         '理論式：φ(r,t) = g(t) ln(b/r) / ln(b/a)。FEMの時間刻み1 ns、解析期間150 ns。時刻0の点は最初の保存値を補って表示。', fontsize=11)
fig.savefig(out/'comparison.png',dpi=170)
fig.savefig(out/'comparison.svg')
print('中央断面の比較点数:', len(radial), '過渡ステップ数:',len(rows))
