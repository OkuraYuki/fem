"""保存データから報告書の日本語グラフのみを再生成する。
実行: python3 tools/localize_report_graphs.py
解析の再実行、CSV・数値結果・過去の検証記録の更新は行わない。
"""
from pathlib import Path
import json
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FONT = "Noto Sans CJK JP"
font_manager.findfont(font_manager.FontProperties(family=FONT), fallback_to_default=False)
plt.rcParams.update({"font.family": FONT, "font.size": 10, "axes.unicode_minus": False,
                     "svg.fonttype": "path", "axes.prop_cycle": matplotlib.cycler(
                         color=["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#444444"])})
warnings.filterwarnings("error", message="Glyph .* missing from current font")
GENERATED = []
MODEL_LABELS = {
    "model_entyu": "周期矩形波",
    "model_entyu_1": "単発矩形波",
    "model_entyu_2": "立下り時間 5 ns",
    "model_entyu_3": "立下り時間 20 ns",
    "model_entyu_teijo": "定常（電極電位 1 V）",
}
CASE_LABELS = {
    "base": "基準：開始時から印加（時間刻み 1 ns）",
    "dt05": "時間刻み 0.5 ns", "dt025": "時間刻み 0.25 ns",
    "delay": "印加開始を 5 ns 遅延", "rise5": "立上り時間 5 ns",
    "fall5": "立下り時間 5 ns", "fall20": "立下り時間 20 ns",
    "high10": "High保持時間 10 ns", "high50": "High保持時間 50 ns",
}
STYLES = ["-", "--", "-.", ":"]
PROBES = ["r100_V", "r500_V", "r1200_V", "r2000_V", "r2400_V"]


def table(path):
    return np.atleast_1d(np.genfromtxt(str(path), delimiter=",", names=True, dtype=None, encoding="utf8"))


def line(ax, x, y, label, index=0, **kwargs):
    return ax.plot(x, y, label=label, linestyle=STYLES[index % 4], **kwargs)


def finish(fig, path, note):
    for ax in fig.axes:
        if ax.axison:
            ax.grid(alpha=.25)
            handles, labels = ax.get_legend_handles_labels()
            if labels:
                ax.legend(fontsize=9, loc="best")
    fig.text(.5, .015, note, ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=(0, .08, 1, .95))
    fig.savefig(str(path.with_suffix(".png")), dpi=160)
    fig.savefig(str(path.with_suffix(".svg")))
    GENERATED.append(str(path.relative_to(ROOT)))
    plt.close(fig)


def pulse():
    folder = OUT / "pulse_study_20261104"
    labels = {"baseline": "基準（1 V・High 25 ns・周期 50 ns）",
              "amplitude_half": "振幅 0.5 V", "amplitude_double": "振幅 2 V",
              "width_short": "High保持時間 10 ns", "width_long": "High保持時間 40 ns",
              "period_short": "周期 35 ns", "period_long": "周期 100 ns",
              "edge_5ns": "立上り・立下り 各5 ns", "edge_10ns": "立上り・立下り 各10 ns"}
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    groups = [
        (["baseline", "amplitude_half", "amplitude_double"], "振幅の比較"),
        (["baseline", "width_short", "width_long"], "High保持時間の比較（周期 50 ns）"),
        (["baseline", "period_short", "period_long"], "周期の比較（High保持時間 25 ns）"),
        (["baseline", "edge_5ns", "edge_10ns"], "立上り・立下りの比較（周期 50 ns）")]
    for ax, (names, title) in zip(axes.flat, groups):
        for i, name in enumerate(names):
            d = table(folder / "output" / name / "probe_timeseries.csv")
            line(ax, d["time_ns"], d["node_181500_V"], labels[name], i)
        ax.set(title=title, xlabel="時間 [ns]", ylabel="観測点の電位 [V]")
    fig.suptitle("円柱モデル：パルス条件による電位応答の違い", fontsize=15)
    finish(fig, folder / "comparison",
           "観測点：節点181500（半径100・高さ750、メッシュ座標）／時間刻み 1 ns・200ステップ\n"
           "傾斜時間の比較ではHigh保持時間25 nsを固定し、Low時間を短縮。外部回路なし。")


def theory():
    folder = OUT / "cylinder_theory_validation_20261102"
    d = table(folder / "steady_all_free_nodes.csv")
    center = np.isclose(d["z"], 750)
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    axes[0].scatter(d["radius"][center], d["fem_V"][center], s=5, label="有限要素法（中央断面）")
    radius = np.geomspace(25, 2500, 400)
    line(axes[0], radius, np.log(2500 / radius) / np.log(100), "同軸円筒の理論式", 1, color="black")
    axes[0].set(title="定常電位：中央断面の半径分布と理論式", xlabel="半径 [メッシュ座標単位]", ylabel="電位 [V]")
    ts = table(folder / "timeseries.csv")
    for i, name in enumerate(list(MODEL_LABELS)[:4]):
        rows = ts[ts["model"] == name]
        line(axes[1], rows["time_ns"], rows["node181254_fem_V"], MODEL_LABELS[name], i)
    axes[1].set(title="過渡応答：入力波形による観測点電位の違い", xlabel="時間 [ns]", ylabel="観測点の電位 [V]")
    finish(fig, folder / "comparison",
           "内側電極1 V・外側電極0 V、半径25〜2500（メッシュ座標）／外部回路なし\n"
           "上：高さ750の自由節点。下：節点181254（半径500・高さ750）、時間刻み1 ns。")


def old_models():
    folder = OUT / "old_models_comparison_20261031"
    fig, axes = plt.subplots(3, 2, figsize=(14, 13))
    series = {}
    for name, ax in zip(MODEL_LABELS, axes.flat):
        d = series[name] = table(folder / "output" / name / "comparison_timeseries.csv")
        if name.endswith("teijo"):
            x = [100, 500, 1200, 2000, 2400]
            for i, prefix in enumerate(["old", "new"]):
                y = [d[prefix + "_node%d_V" % n][0] for n in [181500, 181254, 180954, 180733, 180684]]
                line(ax, x, y, ["修正前", "再計算（10月31日）"][i], i, marker=["o", "x"][i])
            ax.set_xlabel("半径 [メッシュ座標単位]")
        else:
            for i, prefix in enumerate(["old", "new"]):
                line(ax, d["time_ns"], d[prefix + "_node181500_V"], ["修正前", "再計算（10月31日）"][i], i)
            ax.axvline(25, color="gray", ls=":", lw=1)
            ax.set_xlabel("時間 [ns]")
        ax.set(title=MODEL_LABELS[name] + "\n" + name, ylabel="電位 [V]")
    axes.flat[-1].axis("off")
    fig.suptitle("旧円柱モデル5種類：修正前と再計算後の電位比較", fontsize=15)
    finish(fig, folder / "before_after",
           "過渡：節点181500（半径100・高さ750、メッシュ座標）／定常：高さ750の代表5点\n"
           "過渡の時間刻み1 ns。点線は25 nsの切替位置。外部回路なし。")
    for name, case, stem, title in [
        ("model_entyu_1", "base", "node181500_single_pulse_comparison", "単発パルス：修正段階による電位応答の比較"),
        ("model_entyu_3", "fall20", "node181500_fall20_before_after", "立下り時間20 ns：修正段階による電位応答の比較")]:
        d = series[name]
        middle = table(OUT / "poweroff_20261026" / case / "timeseries.csv")
        fig, ax = plt.subplots(figsize=(12, 6))
        line(ax, d["time_ns"], d["old_node181500_V"], "修正前", 0)
        line(ax, middle["time_ns"], middle["r100_V"], "中間修正（10月26日）", 1)
        line(ax, d["time_ns"], d["new_node181500_V"], "再計算（10月31日）", 2)
        if case == "base":
            ax.axvline(25, color="gray", ls=":", label="電極を0 Vに固定（25 ns）")
        else:
            ax.set_xlim(0, 100)
        ax.set(title=title, xlabel="時間 [ns]", ylabel="観測点の電位 [V]")
        finish(fig, folder / stem, "節点181500（半径100・高さ750、メッシュ座標）／時間刻み1 ns・外部回路なし")
    fig, ax = plt.subplots(figsize=(12, 6))
    for i, (name, label) in enumerate([
        ("model_entyu_clamp_sym", "25 nsで電極を0 Vに固定"),
        ("model_entyu_fall_sym", "25〜45 nsで電源電圧を低下"),
        ("model_entyu_open", "25 nsで電源切り離し・媒質を通して放電")]):
        d = table(OUT / name / "probe_timeseries.csv")
        line(ax, d["time_ns"], d["r100_V"], label, i)
    ax.set(title="電源OFF方法の比較：電位を瞬時に下げる場合と緩やかに下げる場合",
           xlabel="時間 [ns]", ylabel="観測点の電位 [V]", xlim=(20, 100))
    finish(fig, folder / "node181500_gradual_decay",
           "節点181500（半径100・高さ750、メッシュ座標）／時間刻み1 ns・高さ方向の同一電位拘束あり")


def poweroff():
    folder = OUT / "poweroff_20261026"
    data = {name: table(folder / name / "timeseries.csv") for name in CASE_LABELS}
    base = data["base"]
    old = np.loadtxt(str(OUT / "model_entyu_1/single_pulse_nodes_timeseries.csv"), delimiter=",", skiprows=1)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for i, col in enumerate(PROBES[:3]):
        radius = col[1:-2]
        color = "C%d" % i
        line(axes[0], old[:, 1], old[:, 3+i], "修正前：半径" + radius, 1, color=color)
        line(axes[0], base["time_ns"], base[col], "修正後：半径" + radius, 0, color=color)
    axes[0].axvline(25, color="gray", ls=":")
    axes[0].set(title="勾配変換・初期条件の修正前後", xlabel="時間 [ns]", ylabel="電位 [V]")
    for i, col in enumerate(PROBES):
        line(axes[1], base["time_ns"], base[col]*1e9, "半径" + col[1:-2], i)
    axes[1].set(title="修正後：電極を0 Vにした後の中央断面", xlabel="時間 [ns]", ylabel="電位 [nV]",
                xlim=(25,100), ylim=(-30,5))
    finish(fig, folder / "before_after", "円柱モデルの中央断面（高さ750、半径はメッシュ座標）／時間刻み1 ns・25 nsで電極を0 Vに固定")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for i, name in enumerate(["base", "dt05", "dt025"]):
        ext = table(folder / name / "extrema_timeseries.csv")
        line(axes[0], ext["time_ns"], ext["global_min_V"]*1e3, CASE_LABELS[name], i)
        line(axes[1], data[name]["time_ns"], data[name]["r1200_V"]*1e9, CASE_LABELS[name], i)
    axes[0].set(title="時間刻みによる全領域最小電位の違い", xlabel="時間 [ns]", ylabel="全領域の最小電位 [mV]", xlim=(25,75))
    axes[1].set(title="中央断面・半径1200の電位応答", xlabel="時間 [ns]", ylabel="電位 [nV]", xlim=(25,75), ylim=(-30,5))
    finish(fig, folder / "time_convergence", "25 nsで電極を0 Vに固定／右：節点180954（半径1200・高さ750、メッシュ座標）")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for i, name in enumerate(["base", "delay", "rise5", "fall5", "fall20", "high10", "high50"]):
        d = data[name]
        ext = table(folder / name / "extrema_timeseries.csv")
        line(axes[0], d["time_ns"], d["r100_V"], CASE_LABELS[name], i)
        line(axes[1], ext["time_ns"], ext["global_min_V"]*1e3, CASE_LABELS[name], i)
    axes[0].set(title="初期条件・電源波形による観測点電位の違い", xlabel="時間 [ns]", ylabel="観測点の電位 [V]")
    axes[1].set(title="初期条件・電源波形による残留負電位の違い", xlabel="時間 [ns]", ylabel="全領域の最小電位 [mV]")
    finish(fig, folder / "waveforms", "左：節点181500（半径100・高さ750、メッシュ座標）／時間刻み1 ns・外部回路なし\n"
           "基準：開始時から1 Vを25 ns保持。凡例の項目を変更し、他の基準条件は共通。")


def spatial():
    folder = OUT / "poweroff_20261029"
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    ref = np.fromfile(str(folder / "collapsed/phi0.bin"), dtype=np.float64)
    for i, (name, dz) in enumerate([("coarse",20), ("full",10), ("fine",5)]):
        field = np.fromfile(str(folder / name / "phi0.bin"), dtype=np.float64).reshape(-1,2409)
        line(axes[0], np.linspace(0,1500,len(field)), (field[:,883]-ref[883])*1e3,
             "高さ方向の刻み " + str(dz), i)
    summary = table(folder / "cause_summary.csv")
    values = [-float(summary[summary["case"] == name]["off25_min_V"][0])*1e3 for name in ["coarse","full","fine"]]
    line(axes[1], [20,10,5], values, "3次元四面体の有限要素法", 0, marker="o")
    line(axes[1], [0,20], [0,values[0]], "原点と粗い分割の結果を結ぶ直線（目安）", 1, color="gray")
    axes[0].set(title="端面付近の不要な高さ方向電位差", xlabel="高さ [メッシュ座標単位]", ylabel="初期電位 − 高さ方向同一電位の参照値 [mV]")
    axes[1].set(title="高さ方向の細分化による負電位の減少", xlabel="高さ方向の刻み [メッシュ座標単位]",
                ylabel="電源OFF直後の負電位ピークの絶対値 [mV]", xlim=(0,21), ylim=(0,.34))
    finish(fig, folder / "spatial_cause", "時間刻み1 ns・1 Vを25 ns保持／左：底面節点884と同じXY位置の初期電位\n"
           "高さ方向同一電位の参照値との差。右：25 ns時点の全領域最小電位の絶対値。")


def circuit():
    folder = OUT / "poweroff_implementation_20261029"
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for i, (name, label) in enumerate([
        ("model_entyu_open", "電源切り離し・媒質を通して放電"),
        ("model_entyu_fall_sym", "電源電圧を20 nsかけて低下"),
        ("model_entyu_clamp_sym", "電極を瞬時に0 Vへ固定")]):
        d = table(OUT / name / "probe_timeseries.csv")
        line(axes[0], d["time_ns"], d["r100_V"], label, i)
    c = table(OUT / "model_entyu_open/circuit_history.csv")
    t = c["time_s"]*1e9
    line(axes[1], t, c["electrode_V"], "有限要素法＋外部回路（時間刻み1 ns）")
    line(axes[1], t, np.exp(-np.maximum(t-25,0)/10), "連続時間の抵抗・容量回路（時定数10 ns）", 1)
    axes[0].set(title="電源OFF時の境界条件による観測点電位の違い", xlabel="時間 [ns]", ylabel="観測点の電位 [V]", xlim=(15,100))
    axes[1].set(title="電源切り離し後の電極電位と理論式", xlabel="時間 [ns]", ylabel="電極の電位 [V]", xlim=(20,100))
    finish(fig, folder / "implemented_discharge", "左：節点181500（半径100・高さ750、メッシュ座標）／25 nsからOFF動作・時間刻み1 ns\n"
           "円柱モデル・高さ方向の同一電位拘束あり。右：切り離し直後は電位を保持してから放電。")
    fig, ax = plt.subplots(figsize=(12,6))
    d = table(OUT / "model_entyu_open/probe_timeseries.csv")
    line(ax, d["time_ns"], d["r100_V"], "観測点の電位")
    ax.axvline(25, color="gray", ls="--", label="電源切り離し（25 ns）")
    ax.set(title="円柱モデル：電源切り離し前後の観測点電位", xlabel="時間 [ns]", ylabel="観測点の電位 [V]")
    finish(fig, OUT / "model_entyu_open/node181500_potential_time",
           "節点181500（半径100・高さ750、メッシュ座標）／時間刻み1 ns・媒質を通して放電")


def main():
    pulse()
    theory()
    old_models()
    poweroff()
    spatial()
    circuit()
    log = ROOT / "build" / "japanese_graphs"
    log.mkdir(parents=True, exist_ok=True)
    (log / "generated.json").write_text(json.dumps(GENERATED, ensure_ascii=False, indent=2)+"\n")
    print("日本語グラフ %d 点をPNG・SVGで生成しました。" % len(GENERATED))


if __name__ == "__main__":
    main()
