# 高さ方向の一様性と電源OFFの扱い（2026-10-29）

通常の`fem` / `fem.exe`から、電圧を下げる動作と、電源を切り離して放電する動作の両方を扱えます。負電位の切り上げや結果の平滑化は行っていません。

## すぐ実行できる円柱モデル

```bash
./bin/fem model_entyu_open
./bin/fem model_entyu_fall_sym
./bin/fem model_entyu_clamp_sym
```

| モデル | 電源OFFの条件 | 結果の保存先 |
|---|---|---|
| model_entyu_open | 25 nsで理想電源を開放。媒質の導電性を通って外側接地電極へ放電 | output/model_entyu_open/ |
| model_entyu_fall_sym | 25〜45 nsで電源電圧を1 Vから0 Vへ線形に低下。既存model_entyu_3の波形を使用 | output/model_entyu_fall_sym/ |
| model_entyu_clamp_sym | 25 nsで電極を瞬時に0 Vへ固定。比較用 | output/model_entyu_clamp_sym/ |

開放の参照ケースでは外部放電抵抗や寄生容量を仮定していません。媒質に既に設定されているσとεだけで放電します。高さ方向の一様性を保つこの条件では、時定数はε/σ_xy=10 nsです。実機に外部抵抗・容量がある場合は、その値を入力する必要があります。

元の入力・結果を保持し、新しい結果は上記の別フォルダに保存しました。同じメッシュと物性を持つ既存のmodel_entyu、model_entyu_1、model_entyu_2、model_entyu_3、model_entyu_teijoにはsymmetry.datを追加しています。既存の出力は再計算していません。

## symmetry.dat：高さ方向の同一電位

```text
axis z
```

同じXY座標を持つ節点を一つの電位自由度にまとめ、縮約した行列で解きます。出力は元の全節点数・節点順に展開されます。材料物性を等方に置き換える処理ではありません。

各高さでXY座標が完全に一致する押出し三角形メッシュが対象です。材料配置・固定電位条件・電流ソースが高さ方向で一致することを検査し、不一致なら停止します。XY座標は入力値の完全一致で照合します。一般の三次元形状に自動適用する機能ではありません。

## 電圧を下げる場合

従来のpulse.datでt_fallを指定します。理想電圧源で電極を駆動する場合はcircuit.datは不要です。

```text
# material v_low v_high t_rise t_high t_fall t_low phase
1 0 1 0 25e-9 20e-9 1000e-9 0
```

この波形では25〜45 nsで電圧が下がります。model_entyu_fall_symでは既存モデルの同じ立下り時間を使用しています。実測した電源波形が別の時間ならt_fallを変更してください。

## circuit.dat：有限の電源抵抗と電源の切り離し

model_entyu_openの設定は以下です。

```text
electrode_material 1
series_resistance_ohm ideal
disconnect_time_s 25e-9
off_resistance_ohm open
external_capacitance_F 0
coordinate_scale_to_m 0.001
initial_state steady
```

| キー | 意味 |
|---|---|
| electrode_material | 電極材料ID。単一の電極回路を指定 |
| series_resistance_ohm | 切り離し前の電源直列抵抗。正の数、ideal（0 Ω）、open（開放） |
| disconnect_time_s | 電源を切り離す絶対時刻。切り離さず波形駆動する場合はnever |
| off_resistance_ohm | 切り離し後に電極から0 V基準へ接続する放電抵抗。正の数、ideal、open。切り離す場合は必須 |
| external_capacitance_F | 電極と0 V基準間に追加する外部容量。省略時0。FEMの場の容量は既に含まれる |
| coordinate_scale_to_m | メッシュ座標からmへの換算係数。回路との結合では必須 |
| initial_state | steady：開始時の導電定常状態。charge_free：初期自由電荷のない誘電場 |

対象材料に対するpulse.datの波形は、電極を直接固定する値ではなく、直列抵抗の電源側電圧になります。電極を構成する全節点を一つの電位にまとめ、FEMの導電電流・変位電流と外部回路のKCLを同時に解きます。切り離し後は元の電源波形から切り離され、off_resistance_ohmで指定した経路へ接続します。

電源抵抗を正の数にしてdisconnect_time_sをneverにすると、有限抵抗で電極を駆動しながら、pulse.datの立下りを与えられます。これにより電圧低下と容量による遅れを組み合わせられます。現在は一度だけ切り離す単一電極回路が対象です。

初期化のsteadyは、解析開始前に現在の電源条件で十分通電していた状態です。charge_freeは、有限電源抵抗で0 Vから充電する検証などに使います。OFF時にはどちらも直前の電位・電荷を保持し、初期化し直しません。

### 切替時刻と時間刻み

切り離し時刻はtime0より後で、時間刻みの格子上に置く必要があります。格子からずれる指定はエラーになります。切替瞬間の電極電位は保持し、その次の時間区間から放電を積分します。例えば25 nsで開放する場合、25 nsでは1 Vを保ち、26 nsから減衰します。

時定数10 ns、dt=1 nsの後退Euler法では、26 nsの電極電位は1/(1+1/10)=0.9090909 Vです。dtを小さくすると連続時間の指数関数へ近づきます。

### 単位と対応範囲

外部抵抗はΩ、外部容量はFです。coordinate_scale_to_mでFEMの導電率行列と誘電率行列の長さスケールを換算してから結合します。参照モデルの0.001はFortran版がヘッダーfact=1000で座標を割る規約に合わせた値です。実物の長さ単位が別なら変更してください。外部抵抗を加えない一様媒質のε/σによる時定数は、一様な長さ換算では変わりません。

回路結合は過渡解析のみ対応し、電極外に0 V固定の基準が必要です。対象電極への明示的な固定境界や別の固定材料との節点共有はエラーになります。非零の分布電流ソース、複数電極回路、任意の初期電荷ファイル、途中時刻からの再開には未対応です。開始時刻は切り離し前に指定します。

## 保存する情報

- tempa.dat001〜：元のメッシュの全節点電位。
- analysis_summary.txt：回路条件、高さ方向の拘束、実際の自由度数、反復情報、指定節点の履歴。
- circuit_history.csv：初期状態を含む時刻、電極電位、接続先基準電位、電極電荷、抵抗枝から電極へ流れ込む電流。source_V欄は切り離し前が電源側、切り離し後が放電先の0 V基準です。

理想電圧源の電流欄は空欄です。0 Ωに対するV/Rを表示していません。openの枝電流は0です。切替時刻の回路履歴には切替直後の枝条件を記録します。

## 検証

```bash
mkdir -p /tmp/fem-electrical-tests
g++ -O2 -std=c++17 tests/physics_regression.cpp src/electrostatic_analyzer.cpp -o /tmp/fem-physics-test
/tmp/fem-physics-test /tmp/fem-electrical-tests
g++ -O2 -std=c++17 tests/circuit_regression.cpp src/electrostatic_analyzer.cpp -o /tmp/fem-circuit-test
/tmp/fem-circuit-test /tmp/fem-electrical-tests
python3 tests/check_electrical_options.py ./bin/fem output/poweroff_implementation_20261029
```

抵抗を通した充電と放電、電源抵抗のある立下り、初期導電定常場、追加容量、単位換算、切替電位の連続性、電荷保存、再計算時の節点順序、全節点出力を理論値と比較しました。円柱モデルも通常のmain.cpp経由で3条件×150ステップを計算しました。
