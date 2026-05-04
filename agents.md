# Electrostatic FEM Solver - Processing Agents & Workflow

## Project Overview
異方性導電率・時間依存電気静力学場をFEM（有限要素法）で解く3Dシミュレータ。

## Governing Equations

### 支配方程式（弱形式）
$$(K_\sigma + \frac{K_\epsilon}{\Delta t})\phi^{n+1} = F + \frac{K_\epsilon}{\Delta t}\phi^n$$

- **$K_\sigma$**：導電率行列（異方性: σx, σy, σz）
- **$K_\epsilon$**：誘電率行列（異方性: εx, εy, εz）
- **$\phi^{n+1}, \phi^n$**：次ステップ・現在ステップの電位
- **$F$**：電流ソース項
- **$\Delta t$**：時間ステップ

時間積分：Backward Euler（陰解法）

---

## Processing Pipeline

### 1. **Input Agent** - データ入力・解析
```
read_mesh()              → メッシュ（座標・要素接続性）
read_material_and_bc()   → 材料特性（σx,σy,σz,εx,εy,εz）+ 境界条件
read_initial_potential() → 初期・固定電位条件
read_analysis_config()   → 解析設定（dt, nstep, t0）
```

**入力ファイル形式**：
- `in.dat`：ノード座標 + 四面体要素
- `sin.dat`：材料ブロック + 材料特性 + 電流ソース
- `tmate.dat`：初期電位（材料単位での固定値）
- `sina.dat`/`sina.dat.t`：時間ステップ・解析期間

### 2. **Assembly Agent** - 行列組立
```
assemble_global_matrix()
├─ 要素レベル計算
│  ├─ Jacobian行列（座標変換）
│  ├─ 勾配ベクトル（基準四面体 → 物理座標）
│  └─ 局所剛性行列ke = (σx·gx + σy·gy + σz·gz)·vol
├─ 局所誘電率行列ke_eps = (εx·gx + εy·gy + εz·gz)·vol
├─ グローバル行列集約
│  ├─ K_total = K_sigma + K_epsilon/dt（全体剛性行列）
│  └─ E（誘電率行列、時間項用）
└─ RHS組立（電流ソース項）
```

**出力**：
- `global_K_*`：CSR形式グローバル剛性行列
- `global_E_*`：CSR形式誘電率行列
- `global_F`：RHS初期値（電流ソース）

### 3. **Boundary Condition Agent** - 境界条件適用
```
apply_boundary_conditions()  → Dirichlet固定値の設定
apply_dirichlet_constraints()→ 行列・RHSへの埋め込み
```

**処理**：
- 固定ノード行：単位行化
- 固定ノード列：RHS調整

### 4. **Solver Agent** - 線形方程式求解
```
solve_linear_system()
├─ 初期値設定（前ステップ解 or 0）
├─ 共役勾配法（CG）
│  ├─ max_iter = max(50, nstep)
│  ├─ tol = 1e-10
│  └─収束判定
└─ 解更新
```

**線形ソルバー**：CG（Conjugate Gradient）
- 対称行列用
- 反復解法
- スパース行列対応

### 5. **Time Integration Agent** - 時間ステッピング
```
solve()  （dt > 0の場合）
└─ for step = 1 to nstep:
   ├─ global_rhs ← global_F + (K_epsilon/dt)·φ^n
   ├─ solve_linear_system()  → φ^(n+1)を求解
   ├─ φ_prev ← φ^(n+1)
   ├─ current_time += dt
   └─ ステップ完了
```

**時間ステップ戦略**：
- 陰解法→安定
- 大きなΔtでも発散なし
- 精度⇔時間の取捨選択

### 6. **Output Agent** - 結果出力
```
write_potential_distribution()
├─ ノード電位を出力形式へ
└─ tempa.dat001（6列フォーマット）
```

---

## Data Structures

### Mesh & Material
```cpp
struct Node {
    double x, y, z;
};

struct Element {
    int nodes[4];      // 四面体要素
    int material;      // 材料ID
};

struct Material {
    double sx, sy, sz; // 導電率（異方性）
    double ex, ey, ez; // 誘電率（異方性）
    double ro, cv;     // レガシー（熱容量等）
};

struct BoundaryCondition {
    int node1, node2;
    int type;
    double value;
};
```

### Linear System
```cpp
// CSR (Compressed Sparse Row) Format
std::vector<int>    global_K_row_ptr;  // 行ポインタ
std::vector<int>    global_K_col_idx;  // 列インデックス
std::vector<double> global_K_values;   // 値

// 誘電率行列
std::vector<int>    global_E_row_ptr;
std::vector<int>    global_E_col_idx;
std::vector<double> global_E_values;

// RHS・解
std::vector<double> global_F;           // F（ソース）
std::vector<double> global_rhs;         // RHS（時間ステップ時）
std::vector<double> potentials;         // φ^(n+1)
std::vector<double> potentials_prev;    // φ^n
```

---

## Call Flow

```
main()
  ├─ read_mesh()
  ├─ read_material_and_bc()
  ├─ read_initial_potential()
  ├─ read_analysis_config()
  └─ solve()
      ├─ apply_boundary_conditions()
      ├─ assemble_global_matrix()
      └─ [time loop]
          ├─ 時間項計算：rhs += (E/dt)·φ_prev
          ├─ solve_linear_system()
          │   ├─ apply_dirichlet_constraints()
          │   └─ CG反復
          ├─ φ_prev ← φ
          └─ write_potential_distribution()（最終ステップ後）
```

---

## Current Status

### ✅ Implemented
- ✅ 異方性導電率（σx, σy, σz個別使用）
- ✅ 異方性誘電率（εx, εy, εz）
- ✅ 時間依存定式化（K_sigma + K_epsilon/dt）
- ✅ 陰解法時間積分（Backward Euler）
- ✅ CSR疎行列フォーマット
- ✅ 境界条件（Dirichlet固定値）
- ✅ 電流ソース項
- ✅ コメント処理（入力ファイルのコメント対応）

### ⚠️ In Progress / Known Issues
- ⚠️ テスト検証：多時間ステップでの精度確認
- ⚠️ パフォーマンス：大規模メッシュ（>100k nodes）での効率化

### ⬜ Future Enhancements
- ⬜ 並列化（OpenMP/MPI）
- ⬜ 可視化モジュール
- ⬜ Neumann/Robin境界条件
- ⬜ 非線形導電率対応
- ⬜ 適応メッシュ細分化
- ⬜ 誤差推定・収束解析

---

## Input File Reference

### `in.dat` - Mesh
```
[num_nodes] [num_elements] [...]
x0 y0 z0
x1 y1 z1
...
[elem node_ids...]
...
```

### `sin.dat` - Material & BC
```
[ntblk] [ntmblk] [ntmate] [npt2] [nlv]    # Header

# Boundary conditions (if ntblk > 0)
[node1] [node2] [type] [value]
...

# Material block assignment (if ntmblk > 0)
[start] [end] [material_id]
...

# Material properties (ntmate entries)
[σx] [σy] [σz] [εx] [εy] [εz]
...

# Current sources (nqt)
[nqt]
[elem] [ix] [iy] [iz]
...
```

### `tmate.dat` - Initial Potential
```
[ninit]
[material_id] [fixed_value]
...
```

### `sina.dat.t` - Analysis Config
```
[dt]      # Time step (≤0 → steady-state)
[nstep]   # Number of steps
[time0]   # Initial time
```

---

## Testing & Validation

**Current Test Case**：
- Mesh：216 nodes, 750 elements
- Material：2 types, anisotropic
- BC：Material 1 fixed @ 1.0V
- Time：dt=1.0e-6, nstep=5

**Expected Results**：
- Steady-state + transient responses
- φ varies smoothly from boundaries
- No divergence/NaN for reasonable inputs

---

## Performance Notes

- **Matrix Assembly**：O(n_elem × 16) for 4-node tet elements
- **CG Solver**：O(nnz × max_iter) per time step
- **Memory**：O(n_nodes + nnz) for sparse matrices
- **Typical Case**：132k nodes → ~2M nonzeros, <1s solve per step

---

## Model Management System (v2.0)

### Architecture

実装は**モデルベースの管理システム**へ移行しました。複数の異なるメッシュ・材料条件を同時に管理・実行できます。

#### ディレクトリ構造

```
fem/
  models/                    ← モデルデータ（入力）
    model_1/
      in.dat               ← メッシュ
      sin.dat              ← 材料・BC
      tmate.dat            ← 初期電位
      sina.dat.t           ← 解析設定（時間ステップ）
    model_2/
      (同じ構造)
    model_3/ ...
  
  output/                    ← 結果出力
    model_1/
      tempa.dat001~010     ← ステップごとの電位分布
      analysis_summary.txt ← 解析レポート
    model_2/
      (同じ構造)
  
  models.dat                 ← モデル指定リスト
  batch_summary.csv          ← バッチ実行結果集約（全モデル統計）
```

### models.dat — モデル指定リスト

```
# Simple model name list
# Each model should have its data in models/<model_name>/
model_1
model_2
model_3
```

- **シンプル形式**：1行に1つのモデル名
- **自動パス解決**：
  - 入力：`models/<model_name>/` から in.dat, sin.dat, tmate.dat, sina.dat.t を自動読み込み
  - 出力：`output/<model_name>/` に自動書き込み

### 実行方法

#### すべてのモデル実行（バッチモード）
```bash
./fem.exe models.dat
```

#### 特定のモデルのみ実行
```bash
./fem.exe models.dat model_1      # model_1 のみ
./fem.exe models.dat model_2      # model_2 のみ
```

#### 単一モデル実行（レガシーモード、models.dat なし）
```bash
./fem.exe
# → models/ の存在チェック
# → ない場合は root の in.dat, sin.dat 等から読み込み
```

### Output Files

各モデルの結果は `output/<model_name>/` に保存：

1. **tempa.dat001～**
   - ステップごとの節点電位分布
   - 6列形式（小数点以下7桁）

2. **analysis_summary.txt**
   - メッシュ情報（ノード数・要素数）
   - 境界条件・材料数
   - 要素範囲（材料ブロック）
   - 収束情報（最終ステップ・時間・最大差分・早期収束フラグ）

3. **batch_summary.csv** （複数モデル実行時のみ）
   ```
   model_name,mesh_file,material_file,output_dir,nodes,elements,...,final_step,final_time,final_max_diff,converged_early
   model_1,models/model_1/in.dat,...,output/model_1,...,10,0.01,9.91e-04,no
   model_2,models/model_2/in.dat,...,output/model_2,...,10,0.01,9.91e-04,no
   ```

### 新しいモデル追加手順

1. **フォルダ作成**
   ```bash
   mkdir models/model_3
   ```

2. **入力ファイル配置**
   ```bash
   cp models/model_1/in.dat models/model_3/
   # sin.dat, tmate.dat, sina.dat.t も同様
   # または異なるメッシュ・材料条件を用意
   ```

3. **models.dat に追加**
   ```
   model_1
   model_2
   model_3  ← 新規
   ```

4. **実行**
   ```bash
   ./fem.exe models.dat model_3
   ```

### Implementation Details

**main.cpp** の主要機能：
- `load_model_specs()`: models.dat を解析、モデル名リストを抽出
- `run_one_model()`: 各モデルをシーケンシャル実行
- `create_directory_recursive()`: output フォルダを自動作成（Windows/Unix対応）
- `join_path()`: フォルダパスの連結

**ElectrostaticAnalyzer** の拡張：
- `set_output_paths()`: 出力フォルダをモデルごとに配置
- `get_*()`: バッチCSV用のメタデータ取得メソッド

### Advantages

✅ **モデル分離**  
- 入力・出力がモデルごとに整理される  
- 上書きの心配なし

✅ **バッチ実行**  
- `models.dat` に追加するだけで複数モデル同時実行  
- 結果は `batch_summary.csv` に集約

✅ **スクリプト化対応**  
- シェルスクリプトや Python で自動実行可能  
- パラメータスイープが容易

✅ **トレーサビリティ**  
- analysis_summary.txt で各モデルの設定記録  
- 入力ファイルパスが batch_summary.csv に記録

### Migration from v1.0

v1.0（単一モデル）から v2.0（マルチモデル）への移行：

```bash
# v1.0: root の in.dat, sin.dat を移動
mkdir -p models/default
mv in.dat sin.dat tmate.dat sina.dat models/default/

# models.dat 作成
echo "default" > models.dat

# 実行
./fem.exe models.dat default
```

カレントディレクトリの in.dat などが存在しなる場合は自動的にレガシーモード（エラー）になります。

---

## Current Status (v2.0)

### ✅ Implemented
- ✅ 異方性導電率（σx, σy, σz個別使用）
- ✅ 異方性誘電率（εx, εy, εz）
- ✅ 時間依存定式化（Backward Euler）
- ✅ CSR疎行列フォーマット
- ✅ 境界条件（Dirichlet固定値）
- ✅ 電流ソース項
- ✅ 収束判定（max_diff < 1e-8 V）
- ✅ **モデル管理システム（v2.0新機能）**
- ✅ **バッチ実行・モデル指定実行**
- ✅ **分析レポート出力**
- ✅ **バッチCSV集約**

### ⚠️ Known Limitations
- ⚠️ シーケンシャル実行のみ（並列化未対応）
- ⚠️ 大規模メッシュ（>500k nodes）での性能未検証

### 🔜 Future Enhancements
- 🔜 モデル間の依存関係・データ共有
- 🔜 並列実行（OpenMP/MPI）
- 🔜 パラメータスイープ自動化
- 🔜 可視化モジュール（ParaView連携）
- 🔜 Neumann/Robin 境界条件
- 🔜 適応メッシュ細分化

