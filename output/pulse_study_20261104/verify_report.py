from pathlib import Path
import csv,json
import numpy as np
D=Path(__file__).resolve().parent
rows=list(csv.DictReader((D/'summary.csv').open()))
checks=[]
for row in rows:
    name=row['case']
    for step in [1,100,200]:
        v=np.fromstring((D/'output'/name/f'tempa.dat{step:03d}').read_text(),sep=' ')
        assert len(v)==363759 and np.isfinite(v).all()
        z_error=float(np.ptp(v.reshape(151,2409),axis=0).max())
        assert z_error==0
        checks.append({'case':name,'step':step,'nodes':len(v),'height_max_difference_V':z_error})
log=(D/'physics_regression.log').read_text();assert 'PASS:' in log
(D/'additional_checks.json').write_text(json.dumps(checks,indent=2))
p=D/'report.md';s=p.read_text()
intro='''## 主な結果

9条件・合計1,800ステップを完走した。節点181500の最大電位は基準条件で0.74768 V、振幅0.5 Vで0.37384 V、2 Vで1.4954 Vとなり、出力丸めの範囲で振幅に比例した。High時間を10 nsから40 nsへ延ばすと、200 ns平均電位は0.149536 Vから0.598144 Vへ増加した。今回の理想電圧駆動では、内部電位は入力波形にほぼ直接追従する。

'''
s=s.replace('## 解析対象と目的',intro+'## 解析対象と目的')
s+='''
## 追加検証と空間精度の注意点

各条件の1・100・200ステップ（計27個の全節点ファイル）を読み、363,759個の電位が有限値であること、および同じXY座標の高さ方向151層で電位差が0 Vであることを確認した。既存の物理回帰試験も通過した（physics_regression.log）。検査記録はadditional_checks.jsonに保存した。

既存の[同じ円柱モデルの理論比較](../cylinder_theory_validation_20261102/report.md)では、代表点の電位に約6.77〜7.32%の過大評価が報告されている。今回も節点181500の1 V入力時の最大電位0.74768 Vに対し、半径a=25、b=2500、観測半径r=100の連続同軸円筒理論値はln(2500/100)/ln(2500/25)=0.698970 Vで、約6.97%高い。今回の結果は条件間の相対比較に使用できるが、絶対電位の高精度評価にはメッシュ改善と収束確認が必要である。

保存電位は指数表記の有限桁出力であり、2 V条件と1 V条件の単純な倍率比較には最大4e-5 V程度の丸め差が含まれる。微小負電位は全条件で絶対値1e-10 V未満であり、有意な負電位の発生として解釈しない。
'''
p.write_text(s)
print('Additional checks passed; report updated')
