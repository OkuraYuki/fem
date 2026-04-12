#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <float.h>
#include <algorithm>
#include <time.h>
#include <signal.h>

#ifndef __PDESOLVER2D
#define __PDESOLVER2D

// ２Ｄ正方領域で偏微分方程式 (PDE: Partial Differential Equation) を
// 解くためのクラス
//
// 領域サイズは n ピクセル x n ピクセル
// pixel に NaN を書き込むと、そこは自由端になる
// mask に非ゼロを書き込むと、そこは固定端になる
// 領域の周囲は必ず自由端もしくは固定端にしておかなければならない
//
// PDESolver2D 自体は抽象クラスなので、子クラスにて calc_new_value 
// を実装して、周囲４点の値から新しい値を求める漸化式を与えないと
// 実体化できない
//
class PDESolver2D
{
protected:
  int n;
  double *pixels;             // NaN を代入すると自由端になる
  int *masks;                 // 非ゼロを代入すると固定端になる
  int interrupted;            // 強制終了する
  double speed;               // 加速、減速定数 (speed > 1 で加速)
public:
  PDESolver2D(int n);         // n x n の正方領域を確保する
  virtual ~PDESolver2D();
  int size()      { return n; }
  double pixel(int p, int q)
                  { return pixels[n*q+p]; }
  void pixel(int p, int q, double v)
                  { pixels[n*q+p] = v; }
  int mask(int p, int q)
                  { return masks[n*q+p]; }
  void mask(int p, int q, int v)
                  { masks[n*q+p] = v; }
  void clear()    { memset(pixels, 0, n*n*sizeof(double));
                    memset(masks, 0, n*n*sizeof(int)); }
  void interrupt(){ interrupted = 1; }
  void set_speed(double s) { speed = s; }
  
  // 全ピクセルで更新差分の相対値が delta 未満になるまで繰り返し計算する
  // 更新差分の相対値とは更新差分をピクセル値で割った値のこと
  // 計算が終了すれば非ゼロを返す
  // 途中でキャンセルされればゼロを返す
  // log がゼロでなければ stderr に表示する進行状況を log にも出力する
  int solve(double delta, FILE *log = NULL);

  // 指定されたファイルにデータを書き込む
  void dump(const char *file_name);

  // 現在の計算値を引き継いで２倍の解像度の領域を確保する
  virtual PDESolver2D *create_double(int delete_this = 0);

protected:  // 子クラスで実装すべき virtual 関数

  // 自身の値 (c) と上下左右 (t, b, l, r) と座標値 (p,q) より、新しい値を計算して返す
  virtual double calc_new_value(double c, double t, double r, double b, double l, int p, int q) = 0;

  // サイズ n のインスタンスを作成して返す
  virtual PDESolver2D *create_instance(int n) = 0;

protected:  // 内部で用いるヘルパー関数

  // create_double 内から呼ばれる。２倍領域のピクセル座標を
  // 元画像に与えると、元画像のピクセル値を補完した値を返す
  double interpolate(int p, int q);

  // 自由端処理のためのヘルパー関数
  double replace_nan(double value, double replace)
                  { return isnan(value) ? replace : value; }
  double replace_nan(int p, int q, double replace)
                  { return replace_nan(pixel(p, q), replace); }

  // 進行状況を表示する
  void log_step(int i, double diff, time_t started_at, FILE *log);

  // 奇数ピクセルあるいは偶数ピクセルだけを計算する
  // ピクセル値の更新差分の相対値の最大値を返す
  double half_step(int even_or_odd);

  // 計算を１回分進める
  // 内部では half_step を２回呼ぶ
  // ピクセル値の更新差分の相対値の最大値を返す
  double step();
};
#endif