#include "electrostatic_analyzer.h"

int main()
{
    ElectrostaticAnalyzer analyzer;

    // ファイル読み込み
    analyzer.read_mesh("in.dat");
    analyzer.read_material_and_bc("sin.dat");  // 導電率と境界条件
    analyzer.read_initial_potential("tmate.dat"); // 初期電位
    analyzer.read_analysis_config("sina.dat");

    // 電位分布を計算
    analyzer.solve();

    // 結果出力
    analyzer.write_potential_distribution("tempa.dat001");

    return 0;
}