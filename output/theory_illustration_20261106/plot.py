"""同軸円筒の理論式を図示する（FEMの再計算は行わない）。"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

plt.rcParams.update({'font.family': 'Noto Sans CJK JP', 'font.size': 11,
                     'axes.titlesize': 13, 'axes.labelsize': 11})
out = Path(__file__).resolve().parent
a, b, tau, off = 25., 2500., 10., 25.
alpha = lambda r: np.log(b / r) / np.log(b / a)
fig, ax = plt.subplots(2, 2, figsize=(15, 10))
fig.subplots_adjust(left=.075, right=.97, bottom=.12, top=.87, hspace=.46, wspace=.25)
fig.suptitle('同軸円筒モデル：電位の空間分布と電源OFF条件による時間応答', fontsize=19, y=.975)
fig.text(.5, .925, '理論式の図示（一様媒質・高さ方向の変化なし・境界に整合した初期状態）／ FEM計算結果ではありません', ha='center')

p = ax[0, 0]
for radius, color in [(1, '#56718b'), (.91, '#e5f0f7'), (.15, '#ed9e41')]:
    p.add_patch(Circle((0, 0), radius, color=color))
p.annotate('中心電極：電圧 g(t)\n半径 a = 25', xy=(0, 0), xytext=(-1.48, .70),
           arrowprops={'arrowstyle': '->'}, fontsize=11)
p.annotate('外側電極：0 V\n内半径 b = 2500', xy=(.95, 0), xytext=(.35, -.98),
           arrowprops={'arrowstyle': '->'}, fontsize=11)
p.plot(.48, .18, 'ko')
p.annotate('媒質内の観測点\n中心軸からの距離 r', xy=(.48, .18), xytext=(.38, .70),
           arrowprops={'arrowstyle': '->'}, fontsize=11)
p.text(-1.4, -1.18, '横断面の模式図（寸法比は実際と異なる）', fontsize=10)
p.set(xlim=(-1.55, 1.55), ylim=(-1.3, 1.2), aspect='equal')
p.axis('off')
p.set_title('① 中心電極と外側接地電極の配置')

p = ax[0, 1]
r = np.linspace(a, b, 1500)
p.plot(r, alpha(r), color='#1565a6', lw=2.5, label='同軸円筒の理論電位（中心電極 1 V）')
rp = np.array([100., 500., 1200., 2000.])
p.plot(rp, alpha(rp), 'o', color='#1565a6')
p.annotate('半径500：0.3495 V', xy=(500, alpha(500)), xytext=(850, .58), arrowprops={'arrowstyle': '->'})
p.text(.38, .78, r'$\phi(r)=\frac{\ln(b/r)}{\ln(b/a)}\times 1\ \mathrm{V}$', transform=p.transAxes, fontsize=16)
p.set(title='② 電位は半径に対して対数的に低下', xlabel='中心軸からの半径 [メッシュ座標単位]', ylabel='電位 [V]', xlim=(0, b), ylim=(0, 1.06))
p.legend(loc='upper right', fontsize=9)
p.grid(alpha=.25)

t = np.linspace(0, 80, 1601)
clamp = np.where(t < off, 1., 0.)
opened = np.exp(-np.maximum(t-off, 0)/tau)
for p, factor, title, ylabel in [
    (ax[1, 0], 1., '③ 電極電圧：0 V固定と電源切り離し', '中心電極の電位 [V]'),
    (ax[1, 1], alpha(500.), '④ 媒質内：半径500の観測点も電極に比例', '半径500の観測点電位 [V]')]:
    p.plot(t, factor*opened, color='#c36718', lw=2.5, label='電源切り離し：媒質を通じて放電')
    p.plot(t, factor*clamp, '--', color='#1565a6', lw=2.2, label='電極を0 V固定：入力に追従')
    p.axvline(off, color='gray', ls=':', lw=1.3)
    p.text(26.5, factor*.88, '25 nsで切替', fontsize=10)
    p.set(title=title, xlabel='時間 [ns]', ylabel=ylabel, xlim=(0, 80), ylim=(-.04*factor, 1.13*factor))
    p.legend(loc='upper right', fontsize=9)
    p.grid(alpha=.25)
ax[1, 0].text(35, .40, r'$V(t)=1\ \mathrm{V}\,e^{-(t-25\,\mathrm{ns})/\tau}$'+'\n時定数 τ = ε / σ = 10 ns', fontsize=12)
ax[1, 1].text(35, .13, r'$\phi(500,t)=0.349485\,V(t)$', fontsize=13)
fig.text(.075, .035, '共通条件：a = 25、b = 2500（同じメッシュ座標単位）、半径方向 σ = 0.01 S/m、ε = 1.0e-10 F/m。\n'
         '時間応答：開始時は1 Vの定常状態、25 nsで切替、外部抵抗・追加容量なし。曲線は連続時間の理論解。', fontsize=11)
fig.savefig(out / 'theory_overview.png', dpi=170)
fig.savefig(out / 'theory_overview.svg')
np.savetxt(out / 'radial_theory.csv', np.column_stack((r, alpha(r))), delimiter=',', header='radius_mesh_unit,potential_V_at_electrode_1V', comments='')
np.savetxt(out / 'time_theory.csv', np.column_stack((t, clamp, opened, alpha(500)*clamp, alpha(500)*opened)), delimiter=',', header='time_ns,electrode_clamp_V,electrode_open_V,r500_clamp_V,r500_open_V', comments='')
