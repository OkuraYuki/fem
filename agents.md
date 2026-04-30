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

