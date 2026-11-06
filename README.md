# fem

異方性の導電率・誘電率を扱う3次元電場の有限要素法（FEM）ソルバーです。
4節点四面体要素を使い、定常解析、時間依存解析、パルス電圧、外部回路との結合に対応します。

報告書は[報告書一覧](reports/README.md)に、説明付きでまとめています。新しい報告書もここに追加します。

## フォルダ構成

```text
fem/
├── README.md           # 最初に読む案内・実行方法
├── CMakeLists.txt      # ビルドと回帰テスト
├── src/                # C++ソース・ヘッダー
├── config/             # バッチ実行用のモデル一覧
│   ├── models.dat
│   └── pulse_models.dat
├── models/             # 入力データ（モデルごとのフォルダ）
├── output/             # 計算結果・過去の検証結果
│   ├── <モデル名>/
│   └── batch_summary.csv
├── reports/            # 報告書へのリンク・説明の一覧
├── docs/               # 入力仕様・理論・検証報告
├── tools/              # 解析・集計・検証の補助ツール
├── tests/              # 物理・回路・CLIの回帰テスト
├── bin/                # ビルドした実行ファイル
│   └── legacy/         # 整理前の実行ファイルを保存
└── build/              # CMakeの中間生成物・テスト結果（自動生成）
```

入力モデル名と既存の計算結果のフォルダ名は維持しています。
大容量の入力メッシュ（`in.dat`など）や解析の生データは、この更新のコミット対象から除外しています。
モデルによっては必要な入力がリポジトリにそろっていないため、実行・再検証時は別途保管したデータを元のパスに配置してください。
報告書・図・集計データは `output/` に含めています。
詳しい資料は[ドキュメント一覧](docs/README.md)から参照できます。

## ビルド

必要なものはC++17対応コンパイラ（g++ / clang++など）です。
以下のコマンドはすべて、このREADMEがある**プロジェクトルート**で実行してください。
CMakeを使う場合はCMake 3.10以上も必要です。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j 2
```

Linux / macOSでは `bin/fem`、Windowsでは `bin/fem.exe` が生成されます。
CMakeを使わない場合も直接ビルドできます。

```bash
# Linux / macOS
mkdir -p bin
g++ -O2 -std=c++17 src/main.cpp src/electrostatic_analyzer.cpp -o bin/fem
```

```powershell
# Windows / PowerShell（g++を使用）
New-Item -ItemType Directory -Force bin
g++ -O2 -std=c++17 src/main.cpp src/electrostatic_analyzer.cpp -o bin/fem.exe
```

`bin/legacy/`は整理前のバイナリの保存先です。現在のソースからビルドしたものを使用してください。

## 実行

```bash
# モデルを1つ実行（一覧への登録不要）
./bin/fem model_1

# config/models.dat に列挙されたモデルを実行
./bin/fem config/models.dat

# 一覧から1つに絞って実行
./bin/fem config/models.dat model_10

# パルス比較用のモデル一覧を実行
./bin/fem config/pulse_models.dat

# 回路結合の例
./bin/fem model_entyu_open
./bin/fem model_entyu_fall_sym
```

Windows PowerShellでは `./bin/fem` を `.\bin\fem.exe` に置き換えます。
モデル一覧の内容は1行1モデル名で、空行と `#` から始まる行は無視されます。
一覧の場所に関係なく、入力・出力パスは**カレントディレクトリを基準**に解決されます。
`bin/`や`docs/`に移動して実行せず、ルートから起動してください。

引数なしの場合は、従来の単一モデル方式でカレントディレクトリの
`in.dat`などを読み、同じ場所に結果を書きます。通常はモデル名を指定してください。

## 入力と出力

入力は `models/<モデル名>/` にまとめます。

| ファイル | 用途 |
|---|---|
| `in.dat` | 必須：メッシュ・要素接続・節点座標 |
| `sin.dat` | 材料物性・境界条件・電流源 |
| `tmate.dat` | 材料ごとの固定電位 |
| `sina.dat.t` | 時間刻み・ステップ数などの解析制御 |
| `pulse.dat` | 任意：パルス電圧 |
| `voltage_difference.dat` | 任意：節点間電位差の記録 |
| `symmetry.dat` | 任意：同一電位の拘束 |
| `circuit.dat` | 任意：外部回路との結合 |

解析制御ファイルには複数の数値形式があるため、
[入力値リファレンス](docs/input_values_reference.md)の書式を参照してください。

通常の出力先は `output/<モデル名>/` です。
`tempa.dat001`以降に節点電位、`analysis_summary.txt`に条件と収束情報を保存します。
電位差の履歴は解析レポートに記録され、回路結合時は `circuit_history.csv`も生成されます。
一覧ファイルを指定した実行では、1モデルに絞った場合も
`output/batch_summary.csv`を更新します。モデル名だけの実行では更新しません。

同じモデルを再実行すると同名ファイルは上書きされます。
以前の長い実行で生成された余分なステップファイルは自動削除されないため、
比較に必要な結果は実行前に別のフォルダへ保存してください。

新しいモデルは `models/<新しい名前>/`を作って入力を配置し、
`./bin/fem <新しい名前>`で実行できます。バッチに含める場合は
`config/models.dat`などへ名前を追加してください。

## 検証と開発

```bash
ctest --test-dir build -C Release --output-on-failure
# CLIの設定検証（Python 3が必要）
python3 tests/check_electrical_options.py ./bin/fem build/test-output
```

CMakeの回帰テストは `build/test-output/`に結果を保存し、既存のモデル結果は更新しません。
補助ツールの用途と依存関係は[tools/README.md](tools/README.md)を参照してください。
VS Codeでは「Build FEM」「Test FEM」「Run model」タスクが利用できます。

## 旧構成からの変更

- ルートのC++ファイル → `src/`
- 解説・報告のMarkdown → `docs/`
- モデル一覧 → `config/`
- 整理前の `fem` / `fem.exe` → `bin/legacy/`
- バッチ集約CSV → `output/batch_summary.csv`
- 実行コマンドは `./fem models.dat` → `./bin/fem config/models.dat`

過去の `output/`内のレポート・実行記録は当時の記録として保持しています。
そこに書かれた旧パスのコマンドを再実行する場合は、上の対応に従ってください。
