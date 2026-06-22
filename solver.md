# Solver Improvement Log

作成日: 2026-05-28
目的: ソルバーの処理時間短縮と大規模モデル対応のための調査・実装履歴

## 高レベル計画
- A: プロファイリング計測を追加してボトルネックを定量化する
- B: ICCG（Incomplete Cholesky + CG）プロトタイプ実装
- C: 時間ステップ行列のCSRマージ最適化（毎ステップの再構築回避）
- D: 並列化（OpenMP）と外部ライブラリ検討（Eigen, SuiteSparse, PETSc）

## 今回の変更履歴
- 2026-05-28: 初期調査と軽微修正（AI実施）
  - `electrostatic_analyzer.cpp` にて以下を実装:
    - コンストラクタで `mesh_node_count` と `mesh_element_count` を初期化
    - `read_mesh` で読み込んだノード/要素数を保存
    - `apply_dirichlet_constraints` を O(nnz) の単一走査に改善
  - プロファイリング計測を追加（実行時に各主要処理の所要時間を printf 出力）:
    - `assemble_global_matrix`
    - `build_sparse_matrix`
    - `build_time_step_matrix`
    - `solve_linear_system`
  - 目的: どの処理が時間とメモリのボトルネックか定量化するため

## 次のアクション（短期）
1. ビルドと小さなモデル実行で計測ログを取得（私がやります）
2. 計測結果を解析して、ICCG を先に実装するか CSR マージ最適化を先に行うか決定

---

注: 私（AI）は今後このリポジトリに対する実装を行うたびにこの `solver.md` を更新します。

## 実行ログ（model_4）
実行日: 2026-05-28

概要:
- Mesh: 216 nodes, 750 elements
- Materials: 2
- Time stepping: dt=1.500e-20, nstep=500

主要タイミング:
- `build_sparse_matrix`: 0.000741 s
- `assemble_global_matrix` (全体): 0.001995 s
- `build_time_step_matrix` (scale=6.666667e+19): 0.000901 s
- `solve_linear_system` (CG): 0.057554 s

ソルバ挙動:
- CG 収束: 96 反復で収束
- 1ステップで収束し、解析は早期収束（max_diff=1.408988e-12）

メモ:
- 現時点では行列組立は非常に速く、ソルバ（CG）の時間が支配的。ICCG や前処理導入の効果が見込めます。

次の推奨アクション:
- ICCG プロトタイプ実装（IC(0) を前処理にして CG を実行）を優先して比較ベンチを取る。
- 併せて `build_time_step_matrix` の CSR マージ最適化も並行で評価可能。

## ICCG 実装結果
実装日: 2026-05-28

概要:
- IC(0) を構築し ICCG（IC(0) 前処理付き CG）を実装・有効化しました。

ベンチ（model_4）:
- ICCG の解法時間: 0.011806 s (CG の 0.057554 s から改善)
- ICCG 反復: 20 反復で収束（従来 CG は 96 反復）

所見:
- IC(0) により反復数が大幅に減り、ソルバ時間が約5倍改善されました（小モデルでの結果）。
- 次は中〜大規模モデル（数千〜数万ノード）での挙動確認とメモリ負荷・factorization 時間の計測が必要です。

## 解析結果への反復表示
実装日: 2026-05-28

変更内容:
- `analysis_summary.txt` にソルバ情報を追記しました。
- 追加項目:
  - `solver_method`
  - `last_solver_iterations`
  - `last_solver_time_sec`
  - `Solver Iterations Per Step` の表

検証:
- `model_4` で `analysis_summary.txt` に反復回数が出力されることを確認済み。
- 出力例: `step,iterations,solve_time_sec` の CSV 形式テーブル。

## model_8 実行確認
実行日: 2026-05-28

結果:
- `./fem_test models.dat model_8` は解析本体に入る前に終了。
- 理由: `models/model_8/in.dat` が存在せず、mesh file not found で弾かれた。

所見:
- 今回の失敗は要素数過多ではなく、入力ファイル欠落によるもの。
- `model_8` を本当に大きいモデルとして試すには、`models/model_8/in.dat` を配置してから再実行が必要。

## model_7 実行確認
実行日: 2026-05-28

概要:
- `model_7` は正常に起動し、解析本体まで到達しました。
- Mesh: 109500 nodes, 616464 elements

主要タイミング:
- `build_sparse_matrix`: 0.346186 s
- `assemble_global_matrix`: 0.731389 s
- `build_time_step_matrix`: 0.518723 〜 0.606903 s/step
- `solve_linear_system` (ICCG): 0.313391 〜 0.424852 s/step

収束状況:
- 各ステップで ICCG は 50 反復まで到達
- `model_7` では `max_iter = max(50, nstep)` の上限に当たっているため、現状の前処理では十分に収束し切っていない可能性が高い

所見:
- 大規模モデルでも「解析開始前に弾かれる」状態ではなく、少なくとも `model_7` は実行できた。
- ボトルネックは引き続き `build_time_step_matrix` とソルバ反復数。
- 次の改善候補は、`build_time_step_matrix` の CSR マージ最適化と、ICCG 前処理の強化（もしくは別前処理）です。

## CSR マージ最適化
実装日: 2026-05-28

変更内容:
- `build_time_step_matrix` を、行ベクタ再構築 + ソート方式から CSR 直マージ方式へ変更。
- 毎ステップの `rows` 生成をやめ、`global_K_row_ptr/global_E_row_ptr` をそのままマージするようにした。

検証（model_4）:
- 旧実装: `build_time_step_matrix` が約 0.000901 s
- 新実装: `build_time_step_matrix` が約 0.000033 s
- 体感で約 27 倍程度短縮

所見:
- 小規模モデルでも明確に効く。
- `model_7` のような大規模ケースでは、毎ステップの行列再構築コスト削減により、全体時間の短縮が期待できる。

## model_7 再ベンチ（CSR 最適化後）
実行日: 2026-05-28

結果:
- `model_7` は引き続き正常に実行できた。
- Mesh: 109500 nodes, 616464 elements

主要タイミング:
- `build_sparse_matrix`: 0.397305 s
- `assemble_global_matrix`: 0.760364 s
- `build_time_step_matrix`: 0.013125 〜 0.023905 s/step
- `solve_linear_system` (ICCG): 0.405772 〜 0.463521 s/step

比較メモ:
- CSR 直マージ前の `build_time_step_matrix` はおよそ 0.55 s/step だったため、約 20〜40 倍程度の改善が見られる。
- ソルバ本体は依然として 50 反復に張り付いているため、残る主要ボトルネックは ICCG の収束性と前処理の質。

出力確認:
- `analysis_summary.txt` に `solver_method`, `last_solver_iterations`, `last_solver_time_sec`, 反復表が出ていることを確認済み。

## 前処理の追加調整（対角スケーリング付き IC）
実装日: 2026-05-28

変更内容:
- IC(0) の前に対角スケーリングを入れて、係数のスケール差を抑えるようにした。
- `model_7` で安定性と反復数の変化を確認した。

検証結果:
- `model_4`: 収束・出力ともに正常。
- `model_7`: `ICCG` の反復数は依然として 50 反復まで到達し、収束性の改善は限定的だった。

所見:
- 前処理の数値安定化は入ったが、`model_7` の主課題はまだ残っている。
- 次は SSOR / symmetric Gauss-Seidel の前処理か、ノード順序付け（RCM など）の導入が有力。

## RCM 順序付け
実装日: 2026-05-28

変更内容:
- 連立方程式を解く前に Reverse Cuthill-McKee (RCM) 順序付けを適用するようにした。
- `global_K`, `global_E`, `global_F`, `fixed_node`, `fixed_value`, `potentials` を同じ順序で並べ替えて解く。
- 出力時は元のノード順へ戻す。

検証結果（model_7）:
- 解析は正常終了した。
- 反復数は引き続き 50 反復に到達し、今回の RCM だけでは収束性改善は限定的だった。
- `final_max_diff` は `1.459713e-04` で、前回よりやや改善。

所見:
- RCM は壊れていないが、ICCG の収束性を劇的に改善するほどではない。
- 次の候補は SSOR / symmetric Gauss-Seidel、あるいはより強い前処理の導入。

## 定常解析ソルバの改善（model_10）
実装・検証日: 2026-06-22

### 背景

- `model_10` は `model_7` と同じメッシュ・材料条件を使い、`dt=0` とした定常解析モデル。
- Mesh: 109500 nodes, 616464 elements
- 導電率は、金属部の `5.8e7` に対して低導電率部が `8.0e-3 〜 1.0e-2` で、最大約 72.5 億倍の差がある。
- 従来の大規模用 SSOR-CG では、相対残差が約 `1.32e-1` で停滞し、5000 反復後も収束しなかった。
- 定常解析において `nstep=500` が線形ソルバの上限 `nstep * 10 = 5000` に流用されていた。

### 変更内容

1. **定常解析のソルバ切り替え**
   - `dt <= 0` の場合は、節点数が 50000 以上でも SSOR-CG ではなく IC(0) 前処理付き CG（ICCG）を使用する。
   - 時間依存解析の大規模モデルは、従来どおり SSOR-CG を使用する。

2. **対称対角スケーリング**
   - 定常線形システムに対し、以下の変数変換を適用する。

     ```text
     S_ii = 1 / sqrt(A_ii)
     A_scaled = S A S
     b_scaled = S b
     x = S y
     ```

   - 変換後の `A_scaled y = b_scaled` を ICCG で解き、解析後に `x = S y` で元の電位へ戻す。
   - 対称性と正定値性を保つ変換のため、CG の適用条件を維持できる。
   - 対称対角スケーリング後の行列に対して IC(0) 前処理を構築する。ソルバ表示は `ICCG (scaled)` とする。

3. **定常解析の反復条件を分離**
   - 反復上限: 2000
   - 相対残差許容値: `1.0e-8`
   - 定常解析では `sina.dat.t` の `nstep` を線形ソルバの反復上限に使用しない。

4. **収束判定とログの明確化**
   - 収束時はソルバ名、最終相対残差、解法時間を表示する。
   - 未収束時は `WARNING: linear solver did not converge` を表示する。
   - ターミナルの反復ログは 1〜10 反復を毎回、以降を 100 反復ごとに表示する。表示を省略した反復も計算自体は実行される。
   - `analysis_summary.txt` に以下を追加する。
     - `last_solver_converged`
     - `last_solver_relative_residual`

### model_10 検証結果

| 項目 | 変更前 | 変更後 |
|---|---:|---:|
| ソルバ | SSOR-CG | ICCG (scaled) |
| 反復数 | 5000 | 318 |
| 最終相対残差 | 約 `1.32e-1` | `9.919235e-9` |
| 解法時間 | 約 28.8〜39.6 s | `1.551057 s` |
| 収束 | no | yes |

最新の `output/model_10/analysis_summary.txt` で、以下を確認した。

```text
solver_method: ICCG (scaled)
last_solver_converged: yes
last_solver_relative_residual: 9.919235e-09
last_solver_iterations: 318
last_solver_time_sec: 1.551057e+00
```

### 注意点

- `final_max_diff` は時間ステップ間の電位差であり、1 回だけ解く定常解析の線形ソルバ精度を示す値ではない。
- 定常解析の精度確認には `last_solver_converged` と `last_solver_relative_residual` を使用する。
- 対角スケーリングは導電率差そのものを物理的に変更する処理ではない。線形方程式を数値的に解きやすい同値な形へ変換している。
- 変更は `experiment/steady-solver-scaling` ブランチで切り分け、実装前の基準状態をコミット `c6f1d20` に保存した。
