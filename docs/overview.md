# FEM 静電界解析ソルバー — システム概要

## 概要

このプログラムは、有限要素法（FEM: Finite Element Method）を用いて  
**3次元静電界（電位分布）** を数値計算するソルバーです。

Poisson 方程式

```
∇·(σ ∇φ) = -I
```

を4節点四面体要素で離散化し、Gauss-Seidel 法で解きます。

---

## ファイル構成

```
fem/
├── main.cpp                    メインエントリポイント
├── electrostatic_analyzer.h    ElectrostaticAnalyzer クラス宣言
├── electrostatic_analyzer.cpp  ElectrostaticAnalyzer クラス実装
├── pdesolver2d.h               PDESolver2D 抽象クラス宣言
├── pdesolver2d.cpp             PDESolver2D クラス実装
├── mat.py                      結果可視化スクリプト（Python）
├── in.dat                      メッシュデータ
├── sin.dat                     導電率・境界条件データ
├── tmate.dat                   初期電位データ
├── sina.dat                    解析設定データ
└── tempa.dat001                計算結果（電位分布）出力先
```

---

## 実行の流れ

```
main()
 │
 ├─ read_mesh("in.dat")              節点座標・要素情報を読み込む
 ├─ read_material_and_bc("sin.dat")  導電率・境界条件・電流源を読み込む
 ├─ read_initial_potential("tmate.dat") 初期電位を読み込む
 ├─ read_analysis_config("sina.dat") 反復回数・時間刻幅を読み込む
 │
 ├─ solve()
 │   ├─ apply_boundary_conditions()  Dirichlet 境界条件を適用
 │   ├─ assemble_global_matrix()     大域剛性行列 K・荷重ベクトル F を組立
 │   └─ solve_linear_system()        Gauss-Seidel 法で K φ = F を解く
 │
 └─ write_potential_distribution("tempa.dat001")  結果を出力
```

---

## 入力ファイル形式

### in.dat（メッシュデータ）

```
<節点数> <要素数> <dummy> <dummy> <単位>
<要素0のノードID×4>
<要素1のノードID×4>
...
<節点0の x y z>
<節点1の x y z>
...
```

### sin.dat（導電率・境界条件）

```
<ntblk> <ntmblk> <ntmate> <npt2> <nlv>
<node1> <node2> <type> <電位値φ>  (境界条件 ntblk 行)
<start> <end> <材料番号>          (要素材料割り当て ntmblk 行)
<σx σy σz ro cv>                  (材料データ ntmate 行)
<nqt>
<要素番号> <Ix Iy Iz>             (電流源 nqt 行)
```

境界条件タイプ:

| type | 意味 |
|------|------|
| 0    | Dirichlet（片端固定電位） |
| 1    | 両端固定 |
| 4    | 標準固定値 |

### tmate.dat（初期電位）

```
<エントリ数>
<材料番号> <電位値>
...
```

### sina.dat（解析設定）

```
<flags × 7>
<params × 7>      params[0] = 最大反復回数
<dt> <time0>
<周波数>
```

---

## 出力ファイル形式

### tempa.dat001（電位分布）

各節点の電位 φ [V] を1行6値の固定小数点形式で出力します。

```
  1.2345E+02  0.0000E+00  ...
```

---

## 結果の可視化

`mat.py` を使ってバイナリ形式の結果ファイルを可視化できます。

```bash
python mat.py
```

numpy と matplotlib が必要です。ファイル名は `mat.py` 内の  
`"solver5after.data"` を実際のファイル名に変更してください。

---

## ビルド方法

```bash
g++ -O2 -o fem main.cpp electrostatic_analyzer.cpp pdesolver2d.cpp -lm
```

## 実行方法

```bash
./fem
```

入力ファイル（`in.dat`, `sin.dat`, `tmate.dat`, `sina.dat`）が  
カレントディレクトリに存在する必要があります。
