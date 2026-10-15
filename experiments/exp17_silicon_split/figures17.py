"""Exp17 figure (results/exp17/fig_exp17.png) and summary table (results/exp17/analysis.json)."""
import glob
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import law  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp17"


def load(pattern):
    return {Path(f).stem: json.loads(Path(f).read_text(encoding="utf-8")) for f in sorted(glob.glob(str(OUT / "runs" / pattern)))}


def main():
    c0 = law.constants()
    runs = load("*.json")
    mech = json.loads((OUT / "mechanism.json").read_text(encoding="utf-8"))
    cache, table = {}, []
    for name, r in runs.items():
        p = r["job"].get("p")
        key = json.dumps(p, sort_keys=True)
        if key not in cache:
            cache[key] = law.constants(p) if p else c0
        c = cache[key]
        h = r["hist"]
        birth = next((x for x in h if x["n_fil"] > 0), None)
        s = r["splits"][0] if r["splits"] else None
        n_before = None
        if s:
            k = next(i for i, x in enumerate(h) if x["vdd"] >= s["vdd"])
            n_before = h[k - 1]["n_fil"]
        table.append(dict(name=name, set=("U1" if name.startswith("U1") else "U2" if name.startswith("U2") else
                                          "N" if name.startswith("N") else "base"),
                          lam=r["job"]["lam"], r_D=r["job"]["r_D"], R_L=r["job"]["R_L"], mm=r["job"]["sigma_scale"],
                          lam_c=law.lambda_c(r["job"]["r_D"], c) if c["regime"] == "Maxwell" else None,
                          pred=r["prediction"]["w_split"], w_obs=s["parent_w"] if s else None, n_before=n_before,
                          birth_w=max(birth["widths_um"]) if birth else None, n_splits=len(r["splits"]), n_max=r["n_max"],
                          wtot_ratio=float(np.median([sum(x["widths_um"]) / (x["I_uA"] / (c["j_on"] * 1e6))
                                                      for x in h if x["n_fil"] > 0 and x["I_uA"] > 50]))))
    (OUT / "analysis.json").write_text(json.dumps(dict(constants=c0, table=table), indent=1, default=float), encoding="utf-8")

    fig, axs = plt.subplots(2, 3, figsize=(18, 10.5))
    # (a) universal curve
    ax = axs[0, 0]
    x = np.linspace(1.02, 5, 400)
    ax.plot(x, 2 * x * np.log(x / (x - 1)), "k-", lw=2, label="kuram: w_s/λ_c = 2x ln[x/(x−1)], x = λ/λ_c")
    mk = {"base": ("o", "C0"), "U1": ("s", "C1"), "U2": ("D", "C2")}
    for t in table:
        if t["w_obs"] is None or t["lam_c"] is None or t["set"] == "N":
            continue
        m, col = mk[t["set"]]
        clean = t["R_L"] >= 4e3 and t["mm"] == 1.0 and t["n_before"] is not None and t["n_before"] <= 2
        if clean:
            ax.plot(t["lam"] / t["lam_c"], t["w_obs"] / t["lam_c"], m, c=col, ms=8, mec="k",
                    mfc=col if t["n_before"] == 1 else "white")
        elif t["R_L"] < 4e3:
            ax.plot(t["lam"] / t["lam_c"], t["w_obs"] / t["lam_c"], m, c=col, ms=5, alpha=0.3)
    for s_, (m, col) in mk.items():
        ax.plot([], [], m, c=col, mec="k", label={"base": "temel set (r_D 5/10/20)", "U1": "U1: SS 30, F_ii 3, g 1e-12",
                                                    "U2": "U2: SS 15, F_ii 30"}[s_])
    ax.plot([], [], "o", c="0.5", mfc="white", mec="k", label="boş işaret: bölünürken ortamda 2. filaman var")
    ax.plot([], [], "o", c="0.5", alpha=0.3, label="sert yük, R_L = 2 kΩ (soluk)")
    ax.axhline(2, c="0.6", ls=":", lw=1)
    ax.set_xlim(1, 4.5)
    ax.set_ylim(0, 7)
    ax.set_xlabel("λ / λ_c")
    ax.set_ylabel("w_s / λ_c")
    ax.set_title("(a) olgun denklem, evrensel eğri: üç cihaz seti, üç r_D aynı eğride\n"
                 "(λ_c = √(2(V_M−V_down)/(r_D j_on)), sabitler cihaz modelinden)", fontsize=10)
    ax.legend(fontsize=7)
    # (b) mechanism histogram
    ax = axs[0, 1]
    groups = {"temiz (R_L 5k)": ("C0", "temiz (R_L 5 kΩ)"), "evrensellik (U)": ("C1", "evrensellik (U1, U2)"),
              "hızlı (R_L 2k)": ("0.6", "sert yük (R_L 2 kΩ)")}   # normal mismatch only; low-leakage case: panel (f)
    bins = np.linspace(0.6, 1.4, 33)
    for g, (col, lab) in groups.items():
        v = [e["ratio"] for e in mech["events"] if e["protocol"] == g and e["mm"] == 1.0]
        if v:
            ax.hist(v, bins=bins, alpha=0.6, color=col, label=f"{lab}: n = {len(v)}, ort. {np.mean(v):.3f} ± {np.std(v):.3f}")
    ax.axvline(1, c="k", lw=1.5)
    ax.set_xlabel("bölünme anında merkez–kenar drenaj farkı / (V_M − V_down)")
    ax.set_ylabel("olay sayısı")
    ax.set_title("(b) mekanizma: bölünme anında merkez–kenar farkı / (V_M − V_down)\n(komşular ve cihaz uçları dahil gerçek düzen; normal uyumsuzluk)", fontsize=10)
    ax.legend(fontsize=7)
    # (c) lineage
    ax = axs[0, 2]
    r = runs["H_lam60_s1"]
    for xh in r["hist"]:
        for cc, w in zip(xh["centres_um"], xh["widths_um"]):
            ax.plot([xh["I_uA"], xh["I_uA"]], [cc - w / 2, cc + w / 2], c="C0", lw=1.2, alpha=0.8)
    for e in r["splits"]:
        ax.plot(e["I_uA"], e["parent_c"], "r*", ms=12)
    ax.set_xlabel("toplam akım I (µA)")
    ax.set_ylabel("konum y (µm)")
    ax.set_title("(c) soy ağacı (λ = 60 µm, W = 600 µm, temiz protokol):\ntek filaman → 2 → 4; kızlar yeniden bölünüyor (★)", fontsize=10)
    # (d) data collapse
    ax = axs[1, 0]
    for name, r in runs.items():
        t = next(t for t in table if t["name"] == name)
        if not (t["R_L"] >= 4e3 and t["mm"] == 1.0 and t["set"] != "N" and np.isfinite(t["pred"])):
            continue
        p = r["job"].get("p")
        c = cache[json.dumps(p, sort_keys=True)]
        xx = [h["I_uA"] * 1e-6 / (c["j_on"] * t["pred"]) for h in r["hist"]]
        ax.step(xx, [h["n_fil"] for h in r["hist"]], where="post", lw=1.2, alpha=0.8)
    xx = np.linspace(0, 6, 200)
    ax.plot(xx, np.maximum(1, xx / 0.75), "k--", lw=1.5, label="kılavuz çizgi N = maks(1, x / 0,75), uydurma değil")
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 6)
    ax.set_xlabel("boyutsuz akım I / (j_on w_s)")
    ax.set_ylabel("filaman sayısı N")
    ax.set_title("(d) filaman sayısı, boyutsuz akıma karşı: bütün temiz koşular (3 cihaz seti)", fontsize=10)
    ax.legend(fontsize=8)
    # (e) Turing-like triplet + array formula
    ax = axs[1, 1]
    r = runs.get("F_lowmm_lam80")
    if r:
        for xh in r["hist"]:
            for cc, w in zip(xh["centres_um"], xh["widths_um"]):
                ax.plot([xh["vdd"], xh["vdd"]], [cc - w / 2, cc + w / 2], c="C1", lw=1.2, alpha=0.8)
        for e in r["splits"]:
            ax.plot(e["vdd"], e["parent_c"], "r*", ms=12)
        ev = next((e for e in mech["events"] if e["run"] == "F_lowmm_lam80"), None)
        if ev:
            ax.text(0.02, 0.02, f"bölünme anında gerçek düzende merkez–kenar farkı {ev['du_mV']:.2f} mV\n"
                                f"(V_M − V_down = {mech['dV_c_mV']:.2f} mV); tek filaman formülü {ev['du_isolated_mV']:.1f} mV verirdi",
                    transform=ax.transAxes, fontsize=8, va="bottom")
    ax.set_xlabel("V_DD (V)")
    ax.set_ylabel("y (µm)")
    ax.set_title("(e) neredeyse düzgün cihaz (uyumsuzluk × 0,1, λ = 80 µm): V_up'ta 3 filaman birden,\neşit aralıklı (Turing benzeri); bölünmede gerçek düzenle merkez–kenar farkı", fontsize=10)
    # (f) validity limit: very low leakage (g_sub = 1e-15) -> body discharge time ~ sweep step
    ax = axs[1, 2]
    for name, col, lab in (("N_ballast_lam60", "C3", "1 s/adım: bölünme yok"), ("N_ballast_lam60_slow", "C2", "30 s/adım: geç bölünme")):
        rn = runs.get(name)
        if not rn:
            continue
        xs = [h["I_uA"] for h in rn["hist"] if h["n_fil"] > 0]
        ws = [max(h["widths_um"]) for h in rn["hist"] if h["n_fil"] > 0]
        ax.plot(xs, ws, c=col, lw=1.5, label=f"{lab} (en geniş filaman)")
    pr = runs.get("N_ballast_lam60")
    if pr:
        ax.axhline(pr["prediction"]["w_split"], c="k", ls="--", lw=1, label=f"yarı statik öngörü {pr['prediction']['w_split']:.0f} µm")
    ax.set_xlabel("I (µA)")
    ax.set_ylabel("filaman genişliği (µm)")
    ax.set_title("(f) geçerlilik sınırı: çok düşük sızıntı (g_sub = 1e-15, c_b/g_sub ≈ 1 s)\ngövde boşalması taramaya yetişemiyor → bölünme gecikiyor", fontsize=10)
    ax.legend(fontsize=7)
    fig.suptitle("Exp17: silikonda (SKY130 nFET, 4 K model) kendiliğinden ve tekrar tekrar bölünen akım filamanları — olgun denklem",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "fig_exp17.png", dpi=95)
    clean = [t for t in table if t["R_L"] >= 4e3 and t["mm"] == 1.0 and t["n_before"] is not None and t["n_before"] <= 2
             and t["w_obs"] and t["lam_c"] and t["set"] != "N"]
    err = [(t["w_obs"] - t["pred"]) / t["pred"] for t in clean]
    print("clean isolated first splits:", len(clean), "mean rel err", round(float(np.mean(err)), 3), "rmse",
          round(float(np.sqrt(np.mean(np.square(err)))), 3))
    for t in sorted(clean, key=lambda t: (t["set"], t["r_D"], t["lam"])):
        print(f"  {t['name']:14s} set={t['set']:4s} lam/lam_c={t['lam']/t['lam_c']:.2f} pred={t['pred']:.1f} obs={t['w_obs']:.1f}")


if __name__ == "__main__":
    main()
