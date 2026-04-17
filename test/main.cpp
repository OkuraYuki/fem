#include "electrostatic_analyzer.h"
#include <cstdio>

int main()
{
    ElectrostaticAnalyzer analyzer;

    // モデル入力ファイル
    analyzer.read_mesh("in.dat");
    analyzer.read_material_and_bc("sin.dat");
    analyzer.read_initial_potential("tmate.dat");

    // 設定ファイルは sina.dat を優先し、無ければ sina.dat.t を使う
    FILE *cfg = fopen("sina.dat", "r");
    if (cfg) {
        fclose(cfg);
        analyzer.read_analysis_config("sina.dat");
    } else {
        analyzer.read_analysis_config("sina.dat.t");
    }

    // 解析実行
    analyzer.solve();

    // 結果出力
    analyzer.write_potential_distribution("tempa.dat001");

    return 0;
}