#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(14, 18))
ax.set_xlim(0, 10)
ax.set_ylim(0, 24)
ax.axis('off')

# フォント設定
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS']

# 色定義
color_input = '#E8F4F8'
color_process = '#B8E0D2'
color_decision = '#FFE5CC'
color_output = '#D4A5A5'
color_fem = '#C9B1FF'

def draw_box(ax, x, y, width, height, text, color, fontsize=10):
    """ボックスを描画"""
    box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                         boxstyle="round,pad=0.1", 
                         edgecolor='black', facecolor=color, linewidth=2)
    ax.add_patch(box)
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize, 
            wrap=True, weight='bold')

def draw_arrow(ax, x1, y1, x2, y2, label=''):
    """矢印を描画"""
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
                           arrowstyle='->', mutation_scale=30, 
                           linewidth=2, color='black')
    ax.add_patch(arrow)
    if label:
        mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mid_x + 0.3, mid_y, label, fontsize=8, style='italic')

# フロー図の描画
y_pos = 23

# スタート
draw_box(ax, 5, y_pos, 2, 0.6, 'START', '#90EE90', 10)
draw_arrow(ax, 5, y_pos - 0.3, 5, y_pos - 1)
y_pos -= 1.5

# main() 開始
draw_box(ax, 5, y_pos, 3, 0.8, 'main()', '#FFD700', 11)
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# ファイル入力
draw_box(ax, 5, y_pos, 3.5, 0.8, 'ElectrostaticAnalyzer\nanalyzer', color_input, 9)
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# read_mesh
draw_box(ax, 5, y_pos, 3, 0.8, 'read_mesh("in.dat")', color_input, 9)
ax.text(7.5, y_pos, '読込: ノード座標\n要素情報', fontsize=7, va='center')
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# read_material_and_bc
draw_box(ax, 5, y_pos, 3.5, 0.8, 'read_material_and_bc\n("sin.dat")', color_input, 9)
ax.text(7.8, y_pos, '読込: 材料特性\n境界条件', fontsize=7, va='center')
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# read_initial_potential
draw_box(ax, 5, y_pos, 3.5, 0.8, 'read_initial_potential\n("tmate.dat")', color_input, 9)
ax.text(7.8, y_pos, '読込: 初期電位\n(材料1に1.0V)', fontsize=7, va='center')
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# read_analysis_config
draw_box(ax, 5, y_pos, 3.5, 0.8, 'read_analysis_config\n("sina.dat")', color_input, 9)
ax.text(7.8, y_pos, '読込: 解析パラメータ\n反復回数など', fontsize=7, va='center')
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# solve()
draw_box(ax, 5, y_pos, 2.5, 0.8, 'solve()', color_process, 10)
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# potentials初期化
draw_box(ax, 5, y_pos, 3, 0.8, 'potentials[]初期化\n全ノード = 0.0', color_process, 9)
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# 境界条件適用
draw_box(ax, 5, y_pos, 3.5, 0.8, 'apply_boundary\n_conditions()', color_process, 9)
ax.text(7.9, y_pos, '初期電位(材料1)と\n境界条件を適用', fontsize=7, va='center')
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# 大域行列組立
draw_box(ax, 5, y_pos, 3.5, 0.8, 'assemble_global\n_matrix()', color_fem, 9)
ax.text(7.9, y_pos, 'FEM: 局所剛性行列\nGauss積分で組立', fontsize=7, va='center')
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# FEM処理詳細
draw_box(ax, 2, y_pos, 3, 1.2, '各要素(e=0~749):\n・Jacobian計算\n・形状関数導関数\n・局所剛性行列Ke\n・大域行列に組立', '#E6CCFF', 8)
ax.text(2, y_pos - 0.8, 'σ∇φ·∇φ', fontsize=7, style='italic', ha='center')
y_pos -= 2

# 線形方程式求解
draw_box(ax, 5, y_pos, 3.5, 0.8, 'solve_linear\n_system()', color_fem, 9)
ax.text(7.9, y_pos, 'Gauss-Seidel反復法\nで連立方程式を求解', fontsize=7, va='center')
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# 反復ループ
draw_box(ax, 5, y_pos, 4, 1, '各反復(iter=0~nstep):\n固定ノードはスキップ\nφ[i] を更新\n収束判定(max_diff < 1e-6)', color_decision, 8)
draw_arrow(ax, 5, y_pos - 0.5, 5, y_pos - 1.3)
y_pos -= 2

# 出力
draw_box(ax, 5, y_pos, 3.5, 0.8, 'write_potential\n_distribution()', color_output, 9)
ax.text(7.9, y_pos, '出力: tempa.dat001\n6列の指数表記', fontsize=7, va='center')
draw_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
y_pos -= 1.6

# 終了
draw_box(ax, 5, y_pos, 2, 0.6, 'END', '#FFB6C6', 10)

# タイトル
ax.text(5, 23.8, 'Electrostatic FEM Solver - Flow Chart', 
        ha='center', fontsize=14, weight='bold')

# 凡例
legend_elements = [
    mpatches.Patch(facecolor=color_input, edgecolor='black', label='ファイル入力'),
    mpatches.Patch(facecolor=color_process, edgecolor='black', label='データ処理'),
    mpatches.Patch(facecolor=color_fem, edgecolor='black', label='FEM計算'),
    mpatches.Patch(facecolor=color_output, edgecolor='black', label='出力'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

plt.tight_layout()
plt.savefig('/home/ohkura/test_c+/fem/Program_Flow.pdf', format='pdf', dpi=300, bbox_inches='tight')
print("Generated: Program_Flow.pdf")
plt.close()

# 詳細フロー図も作成
fig, ax = plt.subplots(1, 1, figsize=(14, 16))
ax.set_xlim(0, 10)
ax.set_ylim(0, 20)
ax.axis('off')

y_pos = 19.5
ax.text(5, y_pos, 'Detailed FEM Calculation Flow', 
        ha='center', fontsize=14, weight='bold')
y_pos -= 0.8

# 大域行列組立のフロー
sections = [
    ("assemble_global_matrix()", [
        "n = nodes.size() (216)",
        "global_K[n×n] = 0",
        "global_F[n] = 0"
    ]),
    ("For each element e (0~749)", [
        "nodes[4]: n0, n1, n2, n3を取得",
        "導電率 σ を材料から取得",
        "Jacobian J[3×3]を計算",
        "detJ = Jacobian行列式",
        "volume = detJ / 6"
    ]),
    ("For each Gauss点(4点積分)", [
        "形状関数導関数を計算",
        "物理座標での勾配 B = J^-1 * dN/dξ",
        "局所剛性行列 Ke += σ * B^T * B * detJ * weight",
    ]),
    ("大域行列に組立", [
        "global_K[gi][gj] += Ke[i][j]",
        "global_F[gi] += 電流源項"
    ]),
]

colors_section = ['#C9B1FF', '#D4C5FF', '#E0D5FF', '#ECE5FF']
for idx, (title, steps) in enumerate(sections):
    draw_box(ax, 5, y_pos, 4, 0.7, title, colors_section[idx], 9)
    y_pos -= 0.9
    for step in steps:
        ax.text(1, y_pos, '• ' + step, fontsize=8, va='center')
        y_pos -= 0.4
    y_pos -= 0.2

# 連立方程式求解のフロー
ax.text(5, y_pos, 'solve_linear_system() - Gauss-Seidel Method', 
        ha='center', fontsize=11, weight='bold')
y_pos -= 0.6

sections2 = [
    ("初期条件", [
        "potentials[] に初期電位と境界条件を適用済み",
        "max_iter = nstep = 100",
        "収束判定: max_diff < 1e-6"
    ]),
    ("反復処理(iter loop)", [
        "max_diff = 0",
        "For each node i (0~215):",
        "  - fixed_node[i]ならスキップ(固定電位)",
        "  - global_K[i][i] == 0ならスキップ",
        "  - sum = global_F[i]",
        "  - sum -= Σ(global_K[i][j] * potentials[j]) (i≠j)",
        "  - new_phi = sum / global_K[i][i]",
        "  - diff = |new_phi - potentials[i]|",
        "  - potentials[i] = new_phi"
    ]),
    ("収束判定", [
        "max_diff < tol なら収束",
        "各反復後に max_diff を監視",
        "反復が進むにつれ電位分布が安定化"
    ]),
]

colors_section2 = ['#FFE5CC', '#FFE5CC', '#FFE5CC']
for idx, (title, steps) in enumerate(sections2):
    draw_box(ax, 5, y_pos, 4, 0.6, title, colors_section2[idx], 9)
    y_pos -= 0.8
    for step in steps:
        ax.text(1, y_pos, '• ' + step, fontsize=7.5, va='center')
        y_pos -= 0.35
    y_pos -= 0.15

plt.tight_layout()
plt.savefig('/home/ohkura/test_c+/fem/Detailed_Calculation_Flow.pdf', format='pdf', dpi=300, bbox_inches='tight')
print("Generated: Detailed_Calculation_Flow.pdf")
plt.close()

print("Done!")
