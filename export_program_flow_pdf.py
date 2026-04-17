#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib import font_manager

# 日本語フォントを優先で探す
font_names = [
    'IPAexGothic', 'IPAPGothic', 'TakaoPGothic', 'Noto Sans CJK JP', 'DejaVu Sans'
]
font_path = None
for name in font_names:
    try:
        font_path = font_manager.findfont(name, fallback_to_default=False)
        break
    except Exception:
        continue

if font_path:
    plt.rcParams['font.family'] = font_manager.FontProperties(fname=font_path).get_name()
else:
    plt.rcParams['font.family'] = 'sans-serif'

text = '''プログラムの動作フロー

1. 初期化フェーズ（main.cpp）

main()
  ├─ ElectrostaticAnalyzer analyzer; // クラスインスタンス生成
  ├─ analyzer.read_mesh("in.dat")             // ノード座標・要素情報を読込
  ├─ analyzer.read_material_and_bc("sin.dat") // 材料・境界条件を読込
  ├─ analyzer.read_initial_potential("tmate.dat") // 初期電位を読込
  ├─ analyzer.read_analysis_config("sina.dat") // 解析パラメータを読込
  ├─ analyzer.solve()                         // 解析実行
  └─ analyzer.write_potential_distribution("tempa.dat001") // 結果出力

2. ファイル読込フェーズ

read_mesh("in.dat")
  - ノード数: 216
  - 要素数: 750
  - 各要素の4節点番号を読込
  - 各節点の x, y, z 座標を読込

read_material_and_bc("sin.dat")
  - 20個の境界条件を読込
  - 10個の材料定義ブロックを読込
  - 2つの材料の導電率を読込
  - 電流源数: 0

read_initial_potential("tmate.dat")
  - materials_id 1 に属する要素の全ノードに 1.0V を設定
  - material 1 に属する要素: 145〜150
  - これらに含まれるノード番号: 37, 38, 39, 40, 43, 46, 49, 50, 54, 73, 74, 76, 79, 82, 85, 86, 90

read_analysis_config("sina.dat")
  - 反復回数 nstep = 100

3. 解析実行フェーズ（solve()）

solve()
  - potentials[] を 0 で初期化
  - apply_boundary_conditions() で固定電位ノードを設定
  - assemble_global_matrix() で大域剛性行列を組立
  - solve_linear_system() で Gauss-Seidel 反復解法を実行

apply_boundary_conditions()
  - 初期電位 (材料1) を fixed_node にマーク
  - fixed_value を 1.0 V に設定
  - type=4 の境界条件は両端を固定 (0.0V)
  - 初期電位があるノードは変更されない

4. FEM計算フェーズ（assemble_global_matrix()）

各要素について:
  - 要素の4節点座標を取得
  - Jacobian 行列 J を計算
  - detJ = det(J) を算出
  - 逆行列 J_inv を計算
  - 体積 volume = |detJ| / 6

Gauss 4点積分で:
  - 形状関数導関数 dN を定義
  - 物理座標勾配 B = J_inv * dN
  - 局所剛性 Ke[i][j] += σ * Σ_k(B[i][k] * B[j][k]) * detJ * weight

大域行列 global_K に組立:
  - global_K[gi][gj] += Ke[i][j]
  - global_F は電流源項を保持

5. 線形方程式求解フェーズ（solve_linear_system()）

Gauss-Seidel 反復:
  - each iteration:
    - fixed_node ならスキップ
    - sum = global_F[i] - Σ_{j≠i} global_K[i][j] * potentials[j]
    - new_phi = sum / global_K[i][i]
    - potentials[i] = new_phi
    - max_diff を更新
  - 収束判定: max_diff < 1e-6

6. 出力フェーズ

write_potential_distribution():
  - potetials[] を 6列の指数表記で tempa.dat001 に書き出し
'''

fig = plt.figure(figsize=(11.7, 16.5))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111)
ax.axis('off')

ax.text(0.01, 0.99, 'Electrostatic FEM Solver - 動作フロー', fontsize=18, weight='bold', va='top')
ax.text(0.01, 0.95, text, fontsize=10, va='top', family='sans-serif')

plt.subplots_adjust(left=0.03, right=0.98, top=0.97, bottom=0.03)
output_path = '/home/ohkura/test_c+/fem/Current_State_Flow_Explanation.pdf'
fig.savefig(output_path, format='pdf', dpi=300)
print(output_path)
