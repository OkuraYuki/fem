# 円筒電極パルス解析の理論式とFEM比較

## 目的

`models/cylinder_electrode` は、中心円柱電極をパルス電圧で駆動し、外側円筒電極を 0 V 接地した同軸円筒モデルである。

このメモでは、緩和時間

```text
tau = epsilon / sigma
```

を使いながら、まず比較用の理論式を整理し、FEM出力との比較条件を固定する。

## 今回使うモデル条件

幾何:

```text
inner electrode radius a = 0.5 mm
outer grounded electrode inner radius b = 4.5 mm
height = 20 mm
comparison plane z = 10 mm
```

境界条件:

```text
phi(a,t) = g(t)
phi(b,t) = 0
```

入力パルス:

```text
period = 50 ns
high = 1 V for 25 ns
low  = 0 V for 25 ns
dt = 1 ns
nstep = 150
```

土壌領域の代表物性:

```text
sigma = 1.0e-2 S/m
epsilon_r = 20
epsilon = epsilon_r * epsilon0 = 1.770838e-10 F/m
tau = epsilon / sigma = 1.770838e-8 s = 17.708 ns
```

土壌物性は含水率・塩分・周波数で大きく変化する。今回は、検証用に「湿り気のある非塩性土壌」程度の仮値として `sigma=0.01 S/m`, `epsilon_r=20` を採用した。

参考として、乾燥土壌の比誘電率はおおむね数程度で、水は約80とされるため、含水土壌ではその中間的な値になる。また、土壌のバルク導電率は含水率や塩分に強く依存し、低含水・非塩性では `0.1 S/m` 未満の値が広く見られる。

## 支配方程式

電場と電位の関係を

```text
E = -grad(phi)
```

とする。導電電流と変位電流を含めた電流密度は

```text
J = sigma E + epsilon dE/dt
```

である。電荷保存から、内部に電流源がない場合は

```text
div(J) = 0
```

したがって

```text
div(sigma E + epsilon dE/dt) = 0
```

である。`E = -grad(phi)` を代入すると

```text
div(sigma grad(phi) + epsilon d/dt grad(phi)) = 0
```

一様媒質で `sigma`, `epsilon` が空間的に一定なら

```text
sigma Laplacian(phi) + epsilon d/dt Laplacian(phi) = 0
```

となる。

この式は、教授との議論で出てきた形

```text
du/dt + (sigma/epsilon) u = f
```

または

```text
du/dt + (1/tau) u = f
```

と対応しており、

```text
tau = epsilon / sigma
```

が電荷緩和時間になる。

## 同軸円筒の静的理論解

有限長端部の影響を無視し、中心高さ付近を無限長同軸円筒として近似する。

円筒対称なので、電位は半径 `r` のみの関数とみなせる。静的なラプラス方程式は

```text
(1/r) d/dr (r dphi/dr) = 0
```

これを積分すると

```text
phi(r) = A ln(r) + B
```

境界条件

```text
phi(a) = V
phi(b) = 0
```

を入れると

```text
phi(r) = V * ln(b/r) / ln(b/a)
```

したがって、時間依存境界電圧 `g(t)` をそのまま入れた瞬時応答近似は

```text
phi_instant(r,t) = alpha(r) g(t)

alpha(r) = ln(b/r) / ln(b/a)
```

である。

電場は

```text
E_r(r,t) = -dphi/dr = g(t) / (r ln(b/a))
```

となる。

## 一次遅れ近似

写真の議論に合わせて、各点の応答を一次遅れとして近似する場合は

```text
tau dy/dt + y = g(t)
```

を使う。周波数領域では

```text
Y(omega) / G(omega) = 1 / (1 + j omega tau)
```

である。したがって、円筒の空間分布を掛けた比較式は

```text
phi_RC(r,t) = alpha(r) y(t)
```

となる。

今回の比較スクリプトでは、FEMと同じ `dt=1 ns` で Backward Euler 離散化した。

```text
y[n+1] = (y[n] + (dt/tau) g[n+1]) / (1 + dt/tau)
```

## 比較点

境界電極上の点は固定値そのものになり、妥当性比較としては自明になりやすい。
そのため、土壌領域の中心高さ `z=10 mm` で、半径方向に3点を選んだ。

```text
near inner soil : r = 0.921053 mm, alpha = 0.721963
middle soil     : r = 2.605263 mm, alpha = 0.248743
near outer soil : r = 4.078947 mm, alpha = 0.044710
```

対応する節点:

```text
near inner soil : node 53182
middle soil     : node 54142
near outer soil : node 54982
```

## 生成ファイル

FEMと理論値の比較は次のスクリプトで作成した。

```powershell
python scripts\compare_cylinder_theory.py
```

出力:

```text
output/cylinder_electrode/theory_comparison.csv
output/cylinder_electrode/theory_comparison.svg
```

`theory_comparison.csv` には、各比較点について以下を出している。

```text
FEM result
instant coaxial theory
RC-filtered coaxial theory
absolute error against instant theory
absolute error against RC theory
```

## 現時点の比較結果

平均絶対誤差:

```text
near inner soil
  instant theory: 0.14685 V
  RC theory     : 0.37515 V

middle soil
  instant theory: 0.16937 V
  RC theory     : 0.26674 V

near outer soil
  instant theory: 0.05024 V
  RC theory     : 0.06851 V
```

最大絶対誤差:

```text
near inner soil
  instant theory: 0.67011 V
  RC theory     : 0.79329 V

middle soil
  instant theory: 0.40335 V
  RC theory     : 0.59123 V

near outer soil
  instant theory: 0.12004 V
  RC theory     : 0.15381 V
```

現状では、単純な一次遅れ近似より、瞬時同軸円筒解の方がFEM結果に近い。
ただし、FEM結果と同軸円筒の理論値にはまだ無視できない差がある。

## 解釈と次の切り分け

理論上、無限長・一様媒質・完全同軸・境界電圧指定の条件では、空間分布は

```text
ln(b/r) / ln(b/a)
```

になるはずである。

そのため、FEM結果がこの分布から大きく外れる場合、主な候補は次のどれかである。

1. 有限長円筒の端部効果
2. 円筒メッシュの四面体分割による異方的な数値誤差
3. 中心電極・外側電極を「体積材料」として固定していることによる境界位置のずれ
4. 要素行列の勾配変換または異方性係数の扱い
5. 時間項 `K_epsilon/dt` と Dirichlet 埋め込みの相互作用

次の検証としては、同じ半径・高さで以下を行うとよい。

```text
1. 定常解析 dt=0 で同軸円筒理論との半径分布を比較
2. 高さを長くして端部効果が減るか確認
3. 半径方向分割と周方向分割を増やして収束性を見る
4. 2D軸対称または1D半径方向の簡易モデルと比較する
```

特に、まず `dt=0` の定常解析で

```text
phi(r) = ln(b/r) / ln(b/a)
```

に近づくかを見るのが重要である。ここで合わない場合、時間応答以前に空間離散化または行列組立側を疑うべきである。

## 参考値の出典

- NASA Technical Reports Server, "Dielectric properties of soils as a function of moisture content": dry soilの比誘電率は小さく、水分増加で大きくなることが示されている。https://ntrs.nasa.gov/api/citations/19750018483/downloads/19750018483.pdf
- NOAA repository report on bulk electrical conductivity and soil water content: 多様な土壌で低いバルク導電率が `0.1 S/m` 未満として観測されることに触れている。https://repository.library.noaa.gov/view/noaa/66597/noaa_66597_DS1.pdf
- SoilSensor.com, "Dielectric Permittivity": 乾燥土壌の比誘電率は典型的に約4、水は25度で約78.54と説明している。https://soilsensor.com/dielectric-permittivity/
