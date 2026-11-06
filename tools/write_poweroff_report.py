from pathlib import Path
import csv,json,shutil,math
root=Path(__file__).resolve().parents[1];out=root/'output/poweroff_20261026'
def read(name):
 with (out/name).open() as f:return list(csv.DictReader(f))
def sci(x):return '{:.6e}'.format(float(x))
status=read('case_status.csv')
if len(status)!=9 or any(s['completed']!='True' or s['warnings']!='0' for s in status):raise SystemExit('All nine completed cases with zero warnings are required.')
metrics=read('metrics.csv');glob=read('global_metrics.csv');before=read('before_after.csv');conv=read('time_convergence.csv')
g={r['case']:r for r in glob};base=[r for r in metrics if r['case']=='base'];r100=base[0]
lines=['# 電源OFF応答の修正・再検証報告書（2026-10-26）','',
'## 結論','',
'円柱モデルの大きな負電位は、修正前の初期化と四面体の勾配変換に問題があった結果であり、旧報告書の「物理的に妥当」とする結論は撤回する。修正後は中央断面の代表5点で大きな負電位が解消した。一方、端面には小さな正負の応答が残る。', '',
'追加の理論確認：全151層でXY座標と材料の節点集合が一致し、明示境界は外側接地電極のみだった（extrusion_theory_check.json）。媒質がz方向に一様、端面が零流束、sigma_x=sigma_y、誘電率が等方という連続体条件では、2次元調和関数h(x,y)を用いたphi(x,y,z,t)=V(t)h(x,y)が支配方程式を満たす。0<=h<=1かつV(t)>=0なので、今回の理想化条件では負電位は生じないはずである。これは今回の形状・物性からの理論的推論である。端面に残る負電位を物理的な端部効果とする説明は根拠が不十分であり、三次元空間離散化または残る実装問題を優先して疑う。時間刻み収束だけではこの疑いは解消しない。','',
'25 nsのOFF直後、r=100の電位は旧結果−0.068750 Vに対して修正後 '+sci(r100['V_first_off'])+' V。全領域最小値は旧結果−0.079881 Vから修正後 '+sci(g['base']['post_off_global_min_V'])+' Vになった。基準ケースの150 ns時点では全領域最小値 '+sci(g['base']['last_global_min_V'])+' V、最大値 '+sci(g['base']['last_global_max_V'])+' V。','',
'## 実施した1〜7の作業','',
'| 項目 | 実施結果 |','|---|---|',
'| 1. 円柱モデル再解析 | 元の363,759節点・2,095,200要素で基準ケース150 nsを完走 |',
'| 2. 修正前後比較 | 同じ5節点・同じ時刻で電位を比較し、CSVと図を保存 |',
'| 3. 時間刻み検証 | 1、0.5、0.25 nsで共通の0〜75 nsを比較 |',
'| 4. 初期条件比較 | 即時ON、5 ns遅延ON、5 ns立上りを比較 |',
'| 5. 電源波形比較 | 立下り0/5/20 ns、High時間10/25/50 nsを比較 |',
'| 6. 追加点検 | 境界・材料割当・単位・収束・節点順序・出力パスを点検し、不具合を修正 |',
'| 7. 報告書更新 | 本報告書、元の単発パルス報告書への訂正、理論文書の補足を更新 |','',
'## モデルと解析条件','',
'材料1を内側パルス電極、材料3を接地電極とする。材料2の導電率は(0.01, 0.01, 0.008) S/m、誘電率は各方向1e-10 F/m。電流ソースは0。電源OFFは電極を0 Vへ固定する条件であり、浮遊電極・外部回路の切断は今回の対象に含めない。','',
'開始時には、自由節点に対して `(K_epsilon phi(time0))_free = 0` を解く無自由電荷の誘電場を設定した。開始前から長時間通電していた導電定常状態とは異なる。時間積分はBackward Euler、線形求解は対称スケーリング付きICCG。全ケースで線形収束警告は0件である。','',
'| ケース | dt [ns] | 終了 [ns] | 立上り [ns] | High保持 [ns] | 立下り [ns] | 位相 [ns] |','|---|---:|---:|---:|---:|---:|---:|']
spec=json.loads((out/'cases.json').read_text())
for c in spec['cases']:lines.append('| '+' | '.join(map(str,c))+' |')
lines+=['','Low保持時間は各ケース1000 nsとし、解析中の再印加を避けた。基準ケースの最初の150 nsの入力は元モデルと同じ。遅延ONでは位相5 ns、立上りケースでは立上り後に25 nsのHighを保持するため、OFFは30 nsとなる。矩形波はステップ終端時刻の境界値を用いる。','',
'## 修正した不具合と検証','',
'1. 四面体の勾配を`J^(-1)`ではなく`J^(-T)`で変換する。せん断した要素の全24節点順序、3方向、異方性導電率・誘電率でアフィン場の厳密エネルギーと一致した。',
'2. 電極のみHigh・内部0 Vの人工的な初期電荷を避け、境界に整合する初期誘電場を解く。均質な小規模モデルでは内部電位が入力の1/3になる厳密解と一致した。パルスなしの固定電圧過渡解析にも適用した。',
'3. type 5・value 0の自然境界は節点を固定しない。非零流束や未対応typeはエラーにする。',
'4. 開始時刻time0を波形と出力に反映する。設定が欠けた場合の配列範囲外参照を防止する。',
'5. 長い出力パスの切り詰めをなくす。同じ解析オブジェクトで再計算できるよう、RCM後の境界節点番号も入力順へ戻す。単一モデル実行で電位差モニタを読み込む。',
'6. CG法で正の微小分母を一律に破綻扱いしない。真に未収束の場合は以後の結果を出さず解析を停止する。','',
'更新した通常実行ファイルで円柱モデルの最初のステップも実行し、363,759節点の表示精度内の値が比較用プログラムと完全一致した。回帰試験は合格。通常実行経路の5試験（正常、非零開始時刻、不正設定、未対応流束、未対応境界type）も合格した。`-Wall -Wextra -Wpedantic`付きビルドで警告はなかった。AddressSanitizerは実行環境のランタイムライブラリが不足しており実施できなかった。','',
'## 修正前後の電位','',
'旧CSVは元の節点出力の1/24/25/50/75/150 nsと照合し、対象5点で差0だった。新しい比較は内部の倍精度値から12桁でCSVへ出力した。従来のtempa.datは有効数字5桁なので、微小な刻み差の比較にはCSVを用いた。','',
'| 25 nsの評価点 | 旧結果 [V] | 修正後 [V] |','|---|---:|---:|']
for r in before:
 if float(r['time_ns'])==25:lines.append('| '+r['probe']+' | '+sci(r['old_V'])+' | '+sci(r['fixed_V'])+' |')
lines+=['','![修正前後の比較](before_after.png)','','## 時間刻みの影響','','| dt [ns] | OFF後の全領域最小値 [V] | 全領域最小値の減衰近似 [ns] |','|---:|---:|---:|']
for c,dt in [('base',1),('dt05',.5),('dt025',.25)]:lines.append('| '+str(dt)+' | '+sci(g[c]['post_off_global_min_V'])+' | '+'{:.4f}'.format(float(g[c]['global_min_tau_fit_ns']))+' |')
fine=abs(float(g['dt025']['post_off_global_min_V']));err1=abs(float(g['base']['post_off_global_min_V'])-float(g['dt025']['post_off_global_min_V']));err05=abs(float(g['dt05']['post_off_global_min_V'])-float(g['dt025']['post_off_global_min_V']))
lines+=['','0.25 ns結果を基準にしたピークの相対差は1 nsで{:.3f}%、0.5 nsで{:.3f}%。'.format(100*err1/fine,100*err05/fine),
'0.25 nsの解そのものを厳密解とは扱わない。3水準の細分化による変化を確認したものであり、ピークの刻み依存性は残る。減衰近似はOFFから50 nsの全領域最小値を指数関数に当てはめた値で、位置が固定された固有モードの厳密時定数ではない。','',
'| 比較対象 | 中央断面5点の最大絶対差（0〜75 ns）[V] | OFF後の最大絶対差 [V] |','|---|---:|---:|']
for case in ['base','dt05']:
 rows=[r for r in conv if r['case']==case]
 lines.append('| '+case+' 対 dt025 | '+sci(max(float(r['max_abs_error_V']) for r in rows))+' | '+sci(max(float(r['post_off_max_abs_error_V']) for r in rows))+' |')
lines+=['','共通の時刻にある値だけを比較し、矩形波の切替点をまたぐ補間は行っていない。','',
'![時間刻み比較](time_convergence.png)','','## 初期条件・立下り時間・通電時間の比較','','| ケース | 最初の完全OFF [ns] | OFF後の全領域最小値 [V] | 最小時刻 [ns] |','|---|---:|---:|---:|']
for case in ['base','delay','rise5','fall5','fall20','high10','high50']:
 r=g[case];lines.append('| '+case+' | '+'{:.3f}'.format(float(r['off_ns']))+' | '+sci(r['post_off_global_min_V'])+' | '+'{:.3f}'.format(float(r['minimum_time_ns']))+' |')
lines+=['','立下りを0→5→20 nsと延ばすと、全領域の負電位ピークの絶対値は約0.168→0.149→0.0951 mVと小さくなった。High保持を10→25→50 nsと延ばすと、約0.102→0.168→0.192 mVと大きくなった。この比較はdt=1 nsであり、ピークの定量値には上記の刻み依存性を含む。','','遅延ONでは全節点0 Vから最初のONを時間積分内で扱う。即時ONでは時刻0の無自由電荷の誘電場を初期化する。Backward Eulerがステップ終端値を使うため、切替時刻の1ステップ分の導電緩和の扱いには離散化差があり、有限dtで完全一致する条件ではない。','',
'![電源波形比較](waveforms.png)','','## 端面の残留応答と解釈の範囲','',
'基準ケース25 nsの最小節点は884、座標は(-50, -86.603, 0)。最大節点は828、座標は(50, 86.603, 0)。単位換算前のメッシュ座標であり、いずれも端面にある。中央断面の代表点のみで全領域の応答を代表させることはできない。','',
'均質かつ導電率テンソルと誘電率テンソルが比例する系では、適切な初期条件のもとで電位が境界波形に即応する例がある。単純にepsilon/sigmaの一次遅れを仮定して、ON時に必ずゆっくり上昇するとは言えない。今回の材料は異方性比が異なり、端部形状と離散化も含む。端面の小さな応答をすべて実験上の物理現象と断定するには、空間メッシュ収束と実回路条件との照合が必要である。','',
'## 追加点検の結果と残る仕様','',
'- 材料割当は全要素でちょうど1回。節点番号は範囲内で、縮退要素は0。電極と接地の共有節点・明示接地との衝突は0。',
'- 導電率・誘電率はいずれも正の有限値。対象5点はすべて材料2内部にある。',
'- メッシュヘッダーのunit=1000は座標換算に適用されていない。今回の零ソース条件では一様な座標スケーリングがK_sigmaとK_epsilonの双方に同じ係数を掛けるため、電位と時定数の比較は変わらない。電場・電流のSI値には元の長さ単位の確認が必要。',
'- 任意の初期自由電荷や事前通電履歴の入力、浮遊電極、外部回路、非零Neumann流束は今回実装していない。電流ソース入力は対象モデルで0のため、その非零入力の物理的な定義・妥当性は本検証の対象外。',
'- パルスありでは早期収束停止を無効にして全ステップを計算。パルスなしの早期停止は従来の最大ステップ差1e-8 Vによる。','',
'## 証拠と再実行','',
'- [全ケースの完了・警告確認](case_status.csv)',
'- [代表点の比較表](before_after.csv)、[代表点の指標](metrics.csv)、[全領域の指標](global_metrics.csv)',
'- [時間刻み比較](time_convergence.csv)、[極値の位置](extrema_locations.csv)、[メッシュ点検](mesh_audit.json)',
'- [通常実行ファイルの全節点照合](installed_cli_full_mesh_check.json)、[回帰試験ログ](regression.log)、[通常実行経路の試験](cli_checks.json)、[入力・ビルド情報](build_manifest.json)',
'- 各ケースのサブフォルダにtimeseries.csv、analysis_summary.txt、run.log、tempa.dat各ステップを保存。',
'- 途中で見つかったCG停止判定の不具合を修正する前の試行はpilot_pre_cg_fixに保存し、最終比較から除外した。元のmodels/model_entyu_1とoutput/model_entyu_1は変更していない。','',
'再計算は、プロジェクト直下で `python3 tools/run_poweroff_suite.py /absolute/new/output/path` を実行する。出力先は新規ディレクトリを指定する。解析は同じモデルから9ケースを作り、最大4プロセスで計算する。','']
(out/'report.md').write_text('\n'.join(lines))
oldreport=root/'docs/model_entyu_1_single_pulse_report.md';backup=out/'single_pulse_report_before_update.md'
if not backup.exists():shutil.copy2(oldreport,backup)
oldreport.write_text('''# 円柱電極モデルの単発パルス応答（修正版・2026-10-26）

初期条件と四面体の勾配変換を修正し、円柱モデルを再解析しました。
旧報告書の「OFF後の大きな負電位が物理的に妥当」という結論は撤回します。

r=100の25 ns時点の電位は、旧結果−0.068750 Vに対して修正後 '''+sci(r100['V_first_off'])+''' Vです。
端面に残る微小な負電位も、今回の理論条件では空間離散化または残る実装問題を優先して疑うべき結果です。

[修正・9ケースの再解析・比較図を含む報告書](output/poweroff_20261026/report.md)を参照してください。

[更新前の報告書](output/poweroff_20261026/single_pulse_report_before_update.md)は検証履歴として保存しています。
''')
p=root/'docs/theory_validation.md';s=p.read_text();marker='## 2026-10-26 電源OFF応答の再検証による訂正'
if marker not in s:p.write_text(s+'\n\n'+marker+'''\n
四面体の勾配は、本コードのJacobian定義では`J^{-T}`で変換する必要がある。
また、電極だけHigh・内部0 Vを初期状態とすると、人工的な初期自由電荷が生じ得る。
修正版では開始時の境界値に整合する無自由電荷の誘電場を初期化する。

旧円柱解析の大きなOFF後負電位とepsilon/sigmaに近い減衰時間だけを根拠に、
物理的妥当性を結論することはできない。修正後の9ケースの計算結果と限界は
[再検証報告書](output/poweroff_20261026/report.md)を参照。
''')
p=root/'docs/power_off_investigation.md';s=p.read_text();marker='## 1〜7の自律検証完了（2026-10-26）'
if marker not in s:p.write_text(s+'\n\n'+marker+'\n\n円柱モデル9ケースの計算、修正前後・時間刻み・初期条件・波形の比較、追加点検と修正、報告書更新を完了。最終結果は[再検証報告書](output/poweroff_20261026/report.md)を参照。\n')
print('Wrote report.md and updated the original report and theory notes.')
