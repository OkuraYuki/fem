# fem

3次元静電界の有限要素法（FEM）ソルバー。ポアソン方程式（∇·(σ∇φ) = 0）を4節点四面体要素で離散化し、節点電位分布を求める。

参考実装：<a href="https://dora.bk.tsukuba.ac.jp/~takeuchi/?サイエンス/２次元電位分布の数値計算#x0404645">筑波大学 竹内研・2次元電位分布の数値計算</a>

---

## クイックスタート

### ビルド
```bash
g++ -std=c++11 main.cpp electrostatic_analyzer.cpp -o fem.exe
```

### 実行

**すべてのモデルを実行**
```powershell
.\fem.exe models.dat
```

**特定のモデルのみ実行**
```powershell
.\fem.exe models.dat model_1
```

---

## ディレクトリ構造

```
fem/
  models/
    model_1/
      in.dat              ← メッシュデータ
      sin.dat             ← 材料・境界条件
      tmate.dat           ← 初期電位
      sina.dat.t          ← 解析設定（時間ステップ）
    model_2/
      in.dat
      sin.dat
      tmate.dat
      sina.dat.t
    model_3/ ...
  output/
    model_1/
      tempa.dat001~010    ← ステップごとの結果
      analysis_summary.txt ← 解析レポート
    model_2/
      ...
  models.dat              ← モデル指定リスト
  main.cpp
  electrostatic_analyzer.cpp
  electrostatic_analyzer.h
  fem.exe
  batch_summary.csv       ← 全モデルの集約結果（バッチ実行時）
```

### models.dat — モデル指定ファイル

```
# Simple model name list
# Each model should have its data in models/<model_name>/
model_1
model_2
model_3
```

- 各行に**モデル名**を記入
- コメント行は `#` で始まる
- 空行は無視される
- 入力ファイル（in.dat, sin.dat, tmate.dat, sina.dat.t）は自動的に `models/<model_name>/` から読み込まれる
- 出力は自動的に `output/<model_name>/` に書き込まれる

---

## ファイル構成

| ファイル | 役割 |
|---|---|
| `main.cpp` | エントリーポイント・batch runner |
| `electrostatic_analyzer.h/.cpp` | FEM解析クラス本体 |
| `models.dat` | モデル指定リスト（実行するモデルのリスト） |

---

## 入力ファイル

各モデルのデータは `models/<model_name>/` に配置します。

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
[材料物性: σx σy σz εx εy εz ρ Cv]        ×ntmate  （εx εy εz は時間依存用）
<電流源数>
[電流源: 要素番号 Ix Iy Iz]                 ×電流源数
```

- 境界条件タイプ：`1` or `4` = Dirichlet（固定電位）、`5` = Neumann（自然境界）
- 異方性導電率（σx, σy, σz）・誘電率（εx, εy, εz）に対応

### `tmate.dat` — 初期電位

```
<エントリ数>
<材料番号> <電位値>  ×エントリ数
```

指定した材料番号に属する全節点を固定電位として登録する。

### `sina.dat.t` — 解析制御パラメータ

```
<dt>        # 時間ステップ（> 0: 時間依存、≤ 0: 定常解析）
<nstep>     # ステップ数
<t0>        # 初期時間
```

時間依存解析時の Backward Euler 陰解法パラメータ。

---

## 出力ファイル

### `output/<model_name>/tempa.dat001～` — 節点電位分布

```
<電位値 1> <電位値 2> ... <電位値 N>  （1行に6列の固定小数点形式）
```

ステップごとに `tempa.dat001`, `tempa.dat002`, ... が生成されます。

### `output/<model_name>/analysis_summary.txt` — 解析レポート

```
=== Analysis Report ===
Mesh:           216 nodes, 750 elements
Boundary Conditions: 0
Materials:      2
[...]
Convergence:    Step 10, t=0.01 s, max_diff=9.91e-04, converged_early=no
```

### `batch_summary.csv` — バッチ実行結果集約（複数モデル実行時）

```
model_name,mesh_file,material_file,initial_file,config_file,output_dir,...,final_step,final_time,final_max_diff,converged_early
model_1,models/model_1/in.dat,...,output/model_1,...,10,0.01,0.000990879,no
model_2,models/model_2/in.dat,...,output/model_2,...,10,0.01,0.000990879,no
```

すべてのモデル実行結果を1行ずつ記録。解析後の比較・検証に使用。

---

## 解析フロー

```
main()
 ├─ load_model_specs(models.dat)         # モデル指定リスト読み込み
 ├─ for each model:
 │   ├─ ElectrostaticAnalyzer analyzer
 │   ├─ read_mesh(models/<name>/in.dat)
 │   ├─ read_material_and_bc(models/<name>/sin.dat)
 │   ├─ read_initial_potential(models/<name>/tmate.dat)
 │   ├─ read_analysis_config(models/<name>/sina.dat.t)
 │   ├─ solve()
 │   │   ├─ apply_boundary_conditions()
 │   │   ├─ assemble_global_matrix()    # K_sigma + K_epsilon/dt
 │   │   ├─ [time loop] for i=1 to nstep:
 │   │   │   ├─ build_time_step_rhs()  # F + (K_epsilon/dt)·φ^n
 │   │   │   ├─ solve_linear_system()  # CG法で φ^(n+1) を求解
 │   │   │   ├─ write_potential_distribution()
 │   │   │   └─ convergence check
 │   │   └─ write_analysis_report()
 │   └─ append to batch_summary.csv
 └─ output all results
```

---

## クラス仕様

### `ElectrostaticAnalyzer`

| メンバ | 型 | 説明 |
|---|---|---|
| `nodes` | `vector<Node>` | 節点座標 (x, y, z) |
| `elements` | `vector<Element>` | 4節点四面体要素（節点番号4つ + 材料ID） |
| `materials` | `vector<Material>` | 材料物性（σx, σy, σz, εx, εy, εz） |
| `boundaries` | `vector<BoundaryCondition>` | 境界条件 |
| `potentials` | `vector<double>` | 各節点の電位 φ |
| `fixed_node` | `vector<bool>` | 固定節点フラグ |
| `fixed_value` | `vector<double>` | 固定節点の電位値 |
| `global_K` | CSR形式 | 全体剛性行列（異方性対応） |
| `global_E` | CSR形式 | 全体誘電率行列（時間項用） |
| `global_F` | `vector<double>` | 荷重ベクトル（電流源） |

### マトリックス記号

- **K_sigma**：導電率行列（σx, σy, σz で異方性）
- **K_epsilon**：誘電率行列（εx, εy, εz で異方性）
- **時間依存定式化**：$(K_\sigma + \frac{K_\epsilon}{\Delta t})\phi^{n+1} = F + \frac{K_\epsilon}{\Delta t}\phi^n$
- **時間積分**：Backward Euler（陰解法）

---

## 数値解法

| 項目 | 内容 |
|---|---|
| 要素タイプ | 4節点四面体（一次要素） |
| 積分 | 解析的（体積 = \|detJ\|/6） |
| 行列形式 | CSR（圧縮スパース行列） |
| ソルバー | 共役勾配法（CG） |
| 収束判定 | max差分 < 1e-8 V または nstep に到達 |
| σ異方性 | σx, σy, σz を個別に使用 |
| ε異方性 | εx, εy, εz を個別に使用（時間項 K_epsilon に反映） |

---

## ビルド・実行

### ビルド
```bash
g++ -std=c++11 main.cpp electrostatic_analyzer.cpp -o fem.exe
```

### 実行例

```bash
# すべてのモデル実行
./fem.exe models.dat

# model_1 のみ実行
./fem.exe models.dat model_1

# model_2 のみ実行
./fem.exe models.dat model_2
```

### 回帰チェック（model_1）
```bash
python3 check_model1_regression.py
```

---

## 新しいモデル追加手順

1. **モデルフォルダ作成**
   ```bash
   mkdir models/model_3
   ```

2. **入力ファイルを配置**
   ```bash
   cp models/model_1/in.dat models/model_3/
   # sin.dat, tmate.dat, sina.dat.t も同様に配置
   ```

3. **models.dat に追加**
   ```
   model_1
   model_2
   model_3  ← 新規追加
   ```

4. **実行**
   ```bash
   ./fem.exe models.dat model_3
   ```

---
- 実装背景: `agents.md`
