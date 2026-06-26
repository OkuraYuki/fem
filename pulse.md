# Pulse Voltage Support

## 概要

`pulse.dat` があるモデルでは、材料IDで指定した固定電位を時間ステップごとに台形パルスとして更新します。

従来の `tmate.dat` は固定電位を与える材料と節点集合を作るために引き続き使います。`pulse.dat` は、その固定材料の電位値を時刻に応じて上書きします。

`pulse.dat` がないモデルは従来通り、`tmate.dat` の定電圧固定として動きます。

## 入力ファイル

モデルディレクトリに任意で `pulse.dat` を置きます。

```text
models/<model_name>/pulse.dat
```

書式は次の通りです。

```text
material_id v_low v_high t_rise t_high t_fall t_low [phase]
```

各列の意味:

```text
material_id : sin.dat の材料ID（1始まり）
v_low       : low 側電圧 [V]
v_high      : high 側電圧 [V]
t_rise      : v_low から v_high への立ち上がり時間 [s]
t_high      : v_high を維持する時間 [s]
t_fall      : v_high から v_low への立ち下がり時間 [s]
t_low       : v_low を維持する時間 [s]
phase       : 位相シフト [s]。省略時は 0
```

周期は以下で決まります。

```text
period = t_rise + t_high + t_fall + t_low
```

コメント行や行末コメントには `#` が使えます。

## 波形

台形波は次の順に進みます。

```text
v_low
  -> t_rise で v_high へ線形上昇
  -> t_high だけ v_high を維持
  -> t_fall で v_low へ線形下降
  -> t_low だけ v_low を維持
  -> 繰り返し
```

`t_rise = 0` かつ `t_fall = 0` にすると、矩形波として使えます。

## model_14 の設定

今回の `models/model_14/pulse.dat` は以下です。

```text
# material_id v_low(V) v_high(V) t_rise(s) t_high(s) t_fall(s) t_low(s) phase(s)
# material 1 remains at 0 V from tmate.dat; material 2 is the pulsed electrode.
2 0.0 1.0 3.0e-9 20.0e-9 3.0e-9 30.0e-9 0.0
```

これは材料2を次の波形で動かします。

```text
0-3 ns     : 0 V から 1 V へ線形立ち上がり
3-23 ns    : 1 V を維持
23-26 ns   : 1 V から 0 V へ線形立ち下がり
26-56 ns   : 0 V を維持
以後繰り返し
```

材料1は `tmate.dat` の 0 V 固定のままなので、材料2が 0 V 基準に対して 0-1 V のパルス電極になります。

## 解析時の挙動

`pulse.dat` が存在する場合、過渡解析では各ステップの時刻 `t = step * dt` に対してパルス電圧を評価し、該当材料に属する固定節点の `fixed_value` を更新します。

この処理は固定節点の値を書き換えるだけなので、行列構築やソルバーに比べて計算コストは小さいです。

## 外側収束判定

パルス電圧が有効な場合、時間ステップ間の早期収束判定は無効になります。

理由は、入力電圧そのものが時間で変化するため、`max_diff < 1e-8` のような「前ステップとの差が小さい」という判定が物理的な終了条件にならないためです。

そのため、`pulse.dat` がある解析は `sina.dat.t` の `nstep` で指定した回数まで必ず実行します。

ログには次のように出ます。

```text
Pulse voltage enabled: early convergence disabled; running all <nstep> steps.
```

## 出力サマリ

`analysis_summary.txt` には `Pulse Voltages (pulse.dat)` セクションが追加されます。

そこには、読み込まれた各パルス条件と、対象になった固定節点数が出力されます。

## 確認済み内容

小さいテストモデルで以下を確認しました。

```text
dt = 1 ns
t_rise = 3 ns
```

この条件で、固定電位が以下のように更新されました。

```text
step 1: 0.333333 V
step 2: 0.666667 V
step 3: 1.000000 V
```

また、`max_diff` が 0 になっても早期停止せず、指定した `nstep` まで実行されることを確認しています。
