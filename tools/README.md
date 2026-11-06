# 補助ツール

通常の解析はルートの[README](../README.md)に従ってください。
ここには過去の電源OFF調査・回路解析の補助ツールがあります。
PythonツールはPython 3を使用し、数値集計にはNumPy、グラフ生成にはMatplotlibが必要です。
個々のスクリプト先頭のimportと出力先も確認してください。

- `run_poweroff_suite.py`：保存済みの9条件を新しい出力先で再計算。
  `python3 tools/run_poweroff_suite.py /absolute/path/to/new-results`
  （g++と `output/poweroff_20261026/cases.json` が必要）
- `check_poweroff_cli.py`：`bin/fem`のCLIを検証し、過去の調査フォルダの検証JSONを更新。
- `analyze_poweroff.py`、`audit_poweroff.py`、`write_poweroff_report.py`：
  過去の調査結果の集計・監査・レポート再生成。
- `analyze_circuit_examples.py`：保存済み回路解析結果の比較・集計。
- `diagnose_extrusion.py`、`locate_poweroff_extrema.py`、
  `summarize_operator_probe.py`：円柱メッシュ・極値・演算子の調査。
- C++ / Fortranファイル：上記調査の専用プローブ。
  Fortran比較の一部は外部の参照実装が必要です。

スクリプトはプロジェクトルートを自分の配置から解決しますが、
保存済みモデル・過去の出力や専用バイナリを前提とするものがあります。
既存の調査結果を更新するツールを通常解析の手順として実行する必要はありません。

## 報告書の日本語グラフを再生成

```bash
python3 tools/localize_report_graphs.py
```

NumPy、Matplotlib、Noto Sans CJK JPフォントが必要です。保存済みCSV・電位データから、報告書の図と同じ解析フォルダの補足図をPNG・SVG形式で生成します。FEM解析、数値集計、報告書本文は更新しません。過去の解析スクリプトで図を作り直した場合も、最後にこのコマンドで日本語の図へ更新してください。
