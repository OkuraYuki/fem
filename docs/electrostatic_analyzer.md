# ElectrostaticAnalyzer クラス解説

## 概要

`ElectrostaticAnalyzer` は、**3次元静電界の電位分布を FEM で解くクラス**です。

支配方程式（Poisson 方程式）:

```
∇·(σ ∇φ) = -I
```

- σ: 導電率テンソル [S/m]
- φ: 電位 [V]
- I: 電流源密度 [A/m³]

4節点四面体要素で空間を離散化し、  
大域剛性行列を組み立てて **Gauss-Seidel 法** で解きます。

---

## データ構造

### Node（節点）

```cpp
struct Node {
    double x, y, z;   // 3次元座標
};
```

### Element（要素）

```cpp
struct Element {
    int nodes[4];   // 4節点四面体の節点ID
    int material;   // 材料番号
};
```

### Material（材料）

```cpp
struct Material {
    double sx, sy, sz;  // 導電率 σx, σy, σz [S/m]
    double ro, cv;      // ダミー（後方互換性）
};
```

### BoundaryCondition（境界条件）

```cpp
struct BoundaryCondition {
    int node1, node2;   // 対象節点ID
    int type;           // 0: 片端固定, 1: 両端固定, 4: 標準固定
    double value;       // 固定電位値 φ [V]
};
```

---

## メンバー一覧

| 名前 | 型 | 説明 |
|------|----|------|
| `nodes` | `vector<Node>` | 全節点の座標 |
| `elements` | `vector<Element>` | 全要素の節点接続 |
| `materials` | `vector<Material>` | 材料特性（導電率） |
| `boundaries` | `vector<BoundaryCondition>` | 境界条件リスト |
| `potentials` | `vector<double>` | 各節点の電位 φ |
| `fixed_node` | `vector<bool>` | 固定電位ノードフラグ |
| `fixed_value` | `vector<double>` | 固定電位値 |
| `current_sources` | `vector<double>` | 各節点の電流源 |
| `global_K` | `vector<vector<double>>` | 大域剛性行列 K |
| `global_F` | `vector<double>` | 大域荷重ベクトル F |
| `dt` | `double` | 時間刻幅（定常解析では参照のみ） |
| `nstep` | `int` | 最大反復回数 |
| `unit` | `double` | 長さの単位 |

---

## メソッド詳細

### `read_mesh(const char *filename)`

**in.dat** からメッシュ情報を読み込みます。

読み込み内容:
1. ヘッダー: 節点数・要素数・単位
2. 要素ごとの4節点ID
3. 節点ごとの3次元座標 (x, y, z)

---

### `read_material_and_bc(const char *filename)`

**sin.dat** から材料特性と境界条件を読み込みます。

読み込み内容:
1. 境界条件リスト（固定電位ノードとその値）
2. 要素ごとの材料番号
3. 材料ごとの導電率 (σx, σy, σz)
4. 電流源データ（要素番号と電流ベクトル）

電流源は各要素の4節点に等分配されます。

---

### `read_initial_potential(const char *filename)`

**tmate.dat** から初期電位を読み込みます。

材料番号と電位値のペアを読み込み、  
該当材料の全要素・全節点に初期電位を設定します。

---

### `read_analysis_config(const char *filename)`

**sina.dat** から解析設定を読み込みます。

主な読み込み内容:
- `dt`: 時間刻幅（定常解析では実質不使用）
- `nstep`: 最大反復回数（Gauss-Seidel の上限）

---

### `solve()`

電位分布の計算を実行します。処理順序:

```
1. apply_boundary_conditions()  境界条件の適用
2. assemble_global_matrix()     大域剛性行列の組立
3. solve_linear_system()        線形方程式を解く
```

---

### `apply_boundary_conditions()`

境界条件リストに基づいて固定電位ノードを設定します。

| type | 処理 |
|------|------|
| 1, 4 | `node1` と `node2` 両方を固定 |
| 0, その他 | `node1` のみを固定 |

固定されたノードの `potentials[i]` に `fixed_value[i]` をセットします。

---

### `assemble_global_matrix()`

**4点 Gauss 積分**を用いて各要素の局所剛性行列 Ke を計算し、  
大域剛性行列 K と荷重ベクトル F に組み立てます。

#### 積分点（四面体 4 点 Gauss 積分）

```
a = (5 - √5) / 20 ≈ 0.1382
b = (5 + 3√5) / 20 ≈ 0.5854
w = 1/24
```

| 積分点 | (ξ, η, ζ) |
|--------|----------|
| 0 | (a, a, a) |
| 1 | (b, a, a) |
| 2 | (a, b, a) |
| 3 | (a, a, b) |

#### 局所剛性行列の計算

```
Ke[i][j] += σ · (B[i] · B[j]) · detJ · weight
```

- B[i]: i番節点の形状関数勾配（物理座標系）
- σ: 材料の導電率 σx（等方性を仮定）
- detJ: Jacobi 行列式（= 6 × 四面体の体積）

#### 形状関数（一次四面体）

| 節点 i | N_i(ξ, η, ζ) | ∂N_i/∂ξ | ∂N_i/∂η | ∂N_i/∂ζ |
|-------|-------------|---------|---------|---------|
| 0 | 1 - ξ - η - ζ | -1 | -1 | -1 |
| 1 | ξ | 1 | 0 | 0 |
| 2 | η | 0 | 1 | 0 |
| 3 | ζ | 0 | 0 | 1 |

---

### `solve_linear_system()`

**Gauss-Seidel 反復法**で連立方程式 K φ = F を解きます。

```
φ_i^(k+1) = (F_i - Σ_{j≠i} K_ij φ_j) / K_ii
```

- 固定電位ノードは更新をスキップ
- 収束判定: 最大更新差分 < 1e-6
- 最大反復回数: `nstep`（`sina.dat` から）

> **注意**: 現状の実装は密行列（full matrix）を使用しているため、  
> 大規模メッシュではメモリと計算時間が節点数の2乗に比例して増大します。

---

### `write_potential_distribution(const char *filename)`

**tempa.dat001** に電位分布を出力します。

- 各行に最大6値
- 書式: `%12.4E`（指数表記、小数4桁）

---

## 計算フローの概念図

```
節点数 n

    [K] {φ} = {F}

    n×n の剛性行列    n×1 の電位ベクトル    n×1 の荷重ベクトル

          ↓ Gauss-Seidel 反復

    各節点の電位 φ_0, φ_1, ..., φ_{n-1}
```
