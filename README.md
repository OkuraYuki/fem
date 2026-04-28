# fem

3次元静電界の有限要素法（FEM）ソルバー。ポアソン方程式（∇·(σ∇φ) = 0）を4節点四面体要素で離散化し、節点電位分布を求める。

参考実装：<a href="https://dora.bk.tsukuba.ac.jp/~takeuchi/?サイエンス/２次元電位分布の数値計算#x0404645">筑波大学 竹内研・2次元電位分布の数値計算</a>

---

## ファイル構成

| ファイル | 役割 |
|---|---|
| `main.cpp` | エントリーポイント |
| `electrostatic_analyzer.h/.cpp` | FEM解析クラス本体 |
| `pdesolver2d.h/.cpp` | 2次元PDE反復ソルバー（参考実装） |
| `electrostatic` | コンパイル済み実行バイナリ |

---

## 入力ファイル

### `in.dat` — メッシュデータ

```
<節点数> <要素数> <ダミー1> <ダミー2> <単位スケール>
<要素定義> ... （4節点の番号を列挙、1-based）
<節点座標> ... （x y z を列挙）
```

### `sin.dat` — 材料・境界条件

```
<ntblk> <ntmblk> <ntmate> <npt2> <nlv>
[境界条件ブロック: 節点1 節点2 タイプ 値]  ×ntblk
[材料ブロック: 開始要素 終了要素 材料番号]  ×ntmblk
[材料物性: σx σy σz ρ]                      ×ntmate
<電流源数>
[電流源: 要素番号 Ix Iy Iz]                 ×電流源数
```

- 境界条件タイプ：`1` or `4` = Dirichlet（固定電位）、`5` = Neumann（自然境界）

### `tmate.dat` — 初期電位

```
<エントリ数>
<材料番号> <電位値>  ×エントリ数
```

指定した材料番号に属する全節点を固定電位として登録する。

### `sina.dat` / `sina.dat.t` — 解析制御パラメータ

フラグ群・反復ステップ数（`nstep`）・時間刻み・収束判定値などを指定する。`sina.dat` が存在する場合はそちらを優先。

---

## 出力ファイル

### `tempa.dat001` — 節点電位分布

```
<電位値> ...  （1行に6列の固定小数点形式、全節点分）
```

---

## 解析フロー

```
main()
 ├─ read_mesh("in.dat")                  # 節点・要素読み込み
 ├─ read_material_and_bc("sin.dat")      # 材料・境界条件読み込み
 ├─ read_initial_potential("tmate.dat")  # 初期電位（固定節点）設定
 ├─ read_analysis_config("sina.dat")     # 解析設定読み込み
 ├─ solve()
 │   ├─ apply_boundary_conditions()      # Dirichlet BC を剛性行列に適用
 │   ├─ assemble_global_matrix()         # 全体剛性行列 K の組み立て
 │   └─ solve_linear_system()            # Gauss-Seidel法で反復求解
 └─ write_potential_distribution("tempa.dat001")
```

---

## クラス仕様

### `ElectrostaticAnalyzer`

| メンバ | 型 | 説明 |
|---|---|---|
| `nodes` | `vector<Node>` | 節点座標 (x, y, z) |
| `elements` | `vector<Element>` | 4節点四面体要素（節点番号4つ + 材料ID） |
| `materials` | `vector<Material>` | 材料物性（σx, σy, σz, ρ, Cv） |
| `boundaries` | `vector<BoundaryCondition>` | 境界条件 |
| `potentials` | `vector<double>` | 各節点の電位 φ |
| `fixed_node` | `vector<bool>` | 固定節点フラグ |
| `fixed_value` | `vector<double>` | 固定節点の電位値 |
| `global_K` | `vector<vector<double>>` | 全体剛性行列（密行列） |
| `global_F` | `vector<double>` | 荷重ベクトル（電流源） |
| `nstep` | `int` | 最大反復回数（デフォルト200） |

### データ構造

```cpp
struct Node     { double x, y, z; };
struct Element  { int nodes[4]; int material; };
struct Material { double sx, sy, sz; double ro, cv; };
struct BoundaryCondition { int node1, node2; int type; double value; };
```

---

## 数値解法

| 項目 | 内容 |
|---|---|
| 要素タイプ | 4節点四面体（一次要素） |
| 積分 | 解析的（体積 = \|detJ\|/6） |
| 行列 | 密行列（全体剛性行列） |
| ソルバー | Gauss-Seidel反復法 |
| 収束判定 | 更新量の最大値 &lt; 1e-8 |
| σの等方化 | σ = (σx + σy + σz) / 3 |

---

## ビルド方法

```bash
g++ -O2 -o electrostatic main.cpp electrostatic_analyzer.cpp pdesolver2d.cpp -lm
```

---

- 技術解説資料: `electrostatic_fem_solver_technical_guide.md`
