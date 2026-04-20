# PDESolver2D クラス解説

## 概要

`PDESolver2D` は、**2次元正方領域で偏微分方程式（PDE）を反復法で解く抽象基底クラス**です。

- 領域サイズ: n × n ピクセル
- 境界条件: 固定端（Dirichlet）または自由端（Neumann）
- 解法: チェッカーボード Gauss-Seidel 反復法

継承クラスは `calc_new_value()` と `create_instance()` を実装することで  
具体的な PDE（Laplace 方程式、Poisson 方程式など）に対応します。

---

## クラス宣言（pdesolver2d.h）

```
PDESolver2D(int n)
```
n×n の正方領域を確保します。  
`pixels` 配列（電位値など）と `masks` 配列（固定端フラグ）が初期化されます。

---

## メンバー一覧

### フィールド

| 名前 | 型 | 説明 |
|------|----|------|
| `n` | `int` | 領域の一辺ピクセル数 |
| `pixels` | `double*` | 各ピクセルの値（NaN = 自由端） |
| `masks` | `int*` | 固定端フラグ（非ゼロ = 固定端） |
| `speed` | `double` | 収束加速定数（デフォルト 1.0） |
| `interrupted` | `int` | 強制終了フラグ |

### ピクセルアクセス

```cpp
double pixel(int p, int q)          // 値を取得（p = x列, q = y行）
void   pixel(int p, int q, double v) // 値を設定
int    mask(int p, int q)            // 固定端フラグを取得
void   mask(int p, int q, int v)     // 固定端フラグを設定
void   clear()                       // 全ピクセルをゼロクリア
```

配列の並びは **行優先（row-major）**: `pixels[n*q + p]`

---

## 主要メソッド

### `solve(double delta, FILE *log = NULL)`

収束するまで反復計算を行います。

- **終了条件**: 全ピクセルの更新差分の相対値が `delta` 未満
- 更新差分の相対値: `|new - old| / (|new| + |old|)`
- 途中で `interrupt()` が呼ばれると計算を中断
- `log` が非 NULL なら進行状況をファイルにも書き出す
- **戻り値**: 正常終了で非ゼロ、中断でゼロ

### `set_speed(double s)`

過緩和（SOR）定数を設定します。

- `s > 1`: 加速（収束が速くなる場合がある）
- `s < 1`: 減速（安定性向上）
- デフォルト: `1.0`（通常の Gauss-Seidel）

### `dump(const char *file_name)`

`pixels` 配列をバイナリファイルに出力します。  
`mat.py` による可視化で使用します。

### `create_double(int delete_this = 0)`

現在の計算値を引き継ぎ、**2倍の解像度**の領域を新たに生成します。

- 低解像度で粗い収束を得てから高解像度に引き継ぐ用途
- `delete_this = 1` にすると元インスタンスを破棄

---

## 反復アルゴリズムの詳細

### チェッカーボード更新（`half_step`）

領域を市松模様（チェッカーボード）に分割し、  
「偶数ピクセル」と「奇数ピクセル」を交互に更新することで  
並列性と収束安定性を両立します。

```
step()
 ├─ half_step(0)  // 偶数ピクセルを更新
 └─ half_step(1)  // 奇数ピクセルを更新
```

### 更新式

```
new_value = calc_new_value(c, t, r, b, l, p, q)
pixel(p, q) += speed * (new_value - old_value)
```

- `c`: 現在値
- `t, r, b, l`: 上・右・下・左の隣接ピクセル値
- 自由端（NaN）の隣接ピクセルは現在値で置換（Neumann 条件）

### 自由端の扱い

`pixels` に `NaN` を書き込むとそのピクセルは計算対象外（自由端）になります。  
自由端の隣接点として使われるときは `replace_nan()` によって  
自身の現在値に置き換えられます（勾配ゼロ = Neumann 条件）。

---

## 継承クラスで実装すべきメソッド

### `calc_new_value(c, t, r, b, l, p, q)`

周囲4点の値から中心ピクセルの新しい値を計算して返します。

例（Laplace 方程式の場合）:
```cpp
double calc_new_value(double c, double t, double r, double b, double l, int p, int q) {
    return (t + r + b + l) / 4.0;
}
```

### `create_instance(int n)`

サイズ `n` の自分と同じ型のインスタンスを生成して返します。  
`create_double()` 内部で使用されます。

---

## 境界条件の設定例

```cpp
PDESolver2D *solver = new MyDerivedSolver(64); // 64×64

// 上辺を固定端（電位 1.0）に設定
for (int p = 0; p < 64; p++) {
    solver->mask(p, 0, 1);
    solver->pixel(p, 0, 1.0);
}

// 下辺を自由端に設定
for (int p = 0; p < 64; p++) {
    solver->pixel(p, 63, NAN);
}

// 収束計算（相対差分 1e-6 以下で終了）
solver->solve(1e-6, stderr);
```
