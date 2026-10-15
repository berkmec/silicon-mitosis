"""Exp15 sekilleri: dogrulama, 1D birlesme, 2D haritalar, boyut karsilastirmasi, 3D kesitler, faz uzayi."""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pic
import structures as S

RES = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                   "results", "exp15_phase_space"))
RUNS = os.path.join(RES, "runs")
plt.rcParams.update({"font.size": 8, "figure.dpi": 110})


def ld(name):
    p = os.path.join(RUNS, name + ".npz")
    return np.load(p) if os.path.exists(p) else None


def fig_validation():
    v = json.load(open(os.path.join(RES, "validation.json"), encoding="utf-8"))
    s = np.load(os.path.join(RES, "validation_series.npz"))
    fig, ax = plt.subplots(2, 2, figsize=(9, 6.5))
    # (a) Langmuir
    a = ax[0, 0]
    kdx = np.linspace(0.02, 1.6, 200)
    h = 4 * np.pi / 64
    for dt, c in ((0.1, "C0"), (0.05, "C1")):
        a.plot(kdx, [pic.langmuir_freq_pred(k / h, h, dt)[1] for k in kdx], c, lw=1,
               label=f"ongoru: CIC + leapfrog, dt={dt}")
    for r in v["langmuir"]:
        a.plot(r["k_dx"], r["omega_fit"], "o", color="C0" if r["dt"] == 0.1 else "C1", ms=6, mfc="none")
    a.axhline(1.0, color="k", lw=0.6, ls=":")
    a.set_xlabel("k dx"); a.set_ylabel("omega / omega_p")
    a.set_title("(a) soguk Langmuir: olculen (daire) vs ongoru")
    a.legend(fontsize=7)
    # (b) iki demet
    a = ax[0, 1]
    a.plot(s["disp_k"], s["disp_cold"], "k-", lw=1, label="soguk dispersiyon koku")
    kw = np.linspace(0.05, 1.2, 231)
    with np.errstate(all="ignore"):
        gw = pic.warm_two_stream_curve(kw, 1.0, 0.2)
    a.plot(kw, gw, "k--", lw=1, label="sicak (v_th=0.2), deney parametresi")
    ts = v["two_stream_cold"]
    a.plot([r["k"] for r in ts["single_mode"]], [r["gamma_meas"] for r in ts["single_mode"]], "o", color="C3",
           ms=6, mfc="none", label="PIC tek-mod tohum")
    mm = [r for r in ts["modes"] if r["gamma_meas"] == r["gamma_meas"]]
    a.plot([r["k"] for r in mm], [r["gamma_meas"] for r in mm], "x", color="C2", label="PIC cok-mod tohum")
    a.set_xlabel("k v0 / omega_p"); a.set_ylabel("gamma / omega_p")
    a.set_title("(b) simetrik soguk iki-demet buyume orani")
    a.legend(fontsize=7)
    # (c) Landau
    a = ax[1, 0]
    t, sig, a1 = s["landau"]
    L = v["landau"]
    a.semilogy(t, np.abs(sig), lw=0.7, label="|E_1| (izdusum)")
    tt = np.linspace(0, 45, 10)
    A0 = np.abs(sig[np.argmin(np.abs(t - 6))]) * 1.3
    a.semilogy(tt, A0 * np.exp(L["theory_gamma"] * (tt - 6)), "k--", lw=1,
               label=f"teori gamma={L['theory_gamma']:.4f}")
    a.set_xlabel("t omega_p"); a.set_ylabel("|E_1|")
    a.set_title(f"(c) Landau k lambda_D=0.5: uyum gamma={L['fit_gamma']:.4f}, omega={L['fit_omega']:.4f}")
    a.legend(fontsize=7)
    # (d) enerji kaymasi
    a = ax[1, 1]
    an = os.path.join(RES, "analysis.json")
    names, vals = [], []
    for r in v["langmuir"]:
        names.append(f"Langmuir m{r['m']} dt{r['dt']}"); vals.append(r["energy_drift_max"])
    names.append("iki-demet (soguk)"); vals.append(ts["energy_drift_max"])
    names.append("Landau"); vals.append(L["energy_drift_max"])
    if os.path.exists(an):
        A = json.load(open(an, encoding="utf-8"))["runs"]
        for nm, r in A.items():
            names.append(nm); vals.append(r["energy_drift_max"])
    a.barh(range(len(vals)), vals, color="C0")
    a.set_xscale("log")
    a.set_yticks(range(len(vals))); a.set_yticklabels(names, fontsize=6)
    a.set_xlabel("max |E_top(t) - E_top(0)| / E_top(0)")
    a.set_title("(d) toplam enerji korunumu")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_validation.png"))
    plt.close(fig)


def fig_1d():
    names = [n for n in ("1d_L100_s1", "1d_L100_s2", "1d_L100_s3", "1d_L40_seed") if ld(n) is not None]
    if not names:
        return
    A = json.load(open(os.path.join(RES, "analysis.json"), encoding="utf-8"))
    ser = np.load(os.path.join(RES, "analysis_series.npz"))
    z = ld(names[0])
    fig = plt.figure(figsize=(10, 6.5))
    a = fig.add_subplot(2, 2, 1)
    phi = z["phi"]
    L = float(z["L"][0])
    a.imshow(phi, aspect="auto", origin="lower", cmap="RdBu_r", extent=[0, L, z["t"][0], z["t"][-1]],
             vmin=-2, vmax=2)
    a.set_xlabel("x / (v0/omega_p)"); a.set_ylabel("t omega_p")
    a.set_title(f"(a) {names[0]}: phi(x,t) (kirmizi = delik)")
    a = fig.add_subplot(2, 2, 2)
    for nm in names:
        t = ser[f"{nm}__t"]
        a.plot(t, ser[f"{nm}__N_th30"], lw=1, label=f"{nm} (esik 0.3)")
    a.set_xscale("log"); a.set_xlabel("t omega_p"); a.set_ylabel("delik sayisi N")
    a.set_title("(b) 1D: delik sayisi N(t) (phi 7/omega_p zaman ortalamali)")
    a.legend(fontsize=6)
    # (c,d) faz uzayi
    for j, tsel in enumerate((30, float(z["ps_t"][-1]))):
        a = fig.add_subplot(2, 2, 3 + j)
        i = int(np.argmin(np.abs(z["ps_t"] - tsel)))
        a.imshow(z["ps_h"][i].T, origin="lower", aspect="auto", extent=[0, L, -2.5, 2.5], cmap="viridis")
        a.set_xlabel("x"); a.set_ylabel("v_x / v0")
        a.set_title(f"({'cd'[j]}) {names[0]} faz uzayi t={z['ps_t'][i]:.0f}")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_1d_merging.png"))
    plt.close(fig)


def fig_2d_maps():
    rows = [n for n in ("2d_L40_wide", "2d_L40_wide_B2", "2d_L40_wide_quietkick", "2d_L40_wide_quietkick_B2",
                        "2d_L40_wide_noseed", "2d_L20_seed") if ld(n) is not None]
    if not rows:
        return
    times = (15, 30, 50, 75, 100, 150)
    fig, ax = plt.subplots(len(rows), len(times), figsize=(1.75 * len(times), 1.85 * len(rows) + 0.3),
                           squeeze=False)
    for i, nm in enumerate(rows):
        z = ld(nm)
        L = z["L"]
        for j, tt in enumerate(times):
            a = ax[i, j]
            k = int(np.argmin(np.abs(z["t"] - tt)))
            if abs(z["t"][k] - tt) > 1:
                a.axis("off"); continue
            a.imshow(z["phi"][k].T, origin="lower", cmap="RdBu_r", vmin=-2, vmax=2,
                     extent=[0, L[0], 0, L[1]])
            a.set_xticks([]); a.set_yticks([])
            a.set_title(f"t={tt}, max={z['phi'][k].max():.2f}", fontsize=6)
            if j == 0:
                a.set_ylabel(nm, fontsize=7)
    fig.suptitle("phi(x,y): demetler x yonunde; kirmizi = elektron deligi (pozitif potansiyel)", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_2d_maps.png"))
    plt.close(fig)


def fig_dimension():
    A = json.load(open(os.path.join(RES, "analysis.json"), encoding="utf-8"))["runs"]
    ser = np.load(os.path.join(RES, "analysis_series.npz"))
    groups = [("1d_L40_seed", "k-"), ("2d_L40_narrow", "C7-"), ("2d_L40_wide", "C3-"), ("2d_L40_wide_B2", "C0-"),
              ("2d_L40_wide_quietkick", "C1-"), ("2d_L40_wide_quietkick_B2", "C4-"),
              ("2d_L40_wide_noseed", "C5-"), ("1d_L20_seed", "k--"), ("2d_L20_seed", "C3--"),
              ("2d_L20_quietkick", "C1--"), ("3d_L20_seed", "C2--"), ("3d_L20_quietkick", "C8--"),
              ("3d_L20_seed_B2", "C9--")]
    groups = [(n, st) for n, st in groups if n in A]
    fig, ax = plt.subplots(2, 3, figsize=(13, 7))
    for nm, st in groups:
        t = ser[f"{nm}__t"]
        WE = ser[f"{nm}__WE"]
        area = np.prod(A[nm]["L"])
        ax[0, 0].semilogy(t, WE / area, st, lw=1, label=nm)
        ax[0, 1].plot(t, ser[f"{nm}__N_th30"], st, lw=1, label=nm)
        if f"{nm}__C" in ser:
            ax[0, 2].plot(t, ser[f"{nm}__C"], st, lw=1, label=nm)
            Kp = ser[f"{nm}__Kperp"]
            ax[1, 0].plot(t, (Kp - Kp[0]) / area, st, lw=1, label=nm)
        if f"{nm}__kink" in ser:
            ax[1, 1].plot(t, ser[f"{nm}__kink"], st, lw=1, label=nm)
    ax[0, 0].set_ylabel("W_E / hacim"); ax[0, 0].set_title("(a) alan enerjisi yogunlugu (delik gucu)")
    ax[0, 1].set_ylabel("N (phi > 0.3 phi_ref)"); ax[0, 1].set_title("(b) yapi sayisi")
    ax[0, 2].set_ylabel("C = W(k_perp=0)/W"); ax[0, 2].set_title("(c) enine uyum (1 = duzlemsel)")
    ax[1, 0].set_ylabel("(K_perp - K_perp(0)) / hacim"); ax[1, 0].set_title("(d) dik isinma yogunlugu")
    ax[1, 1].set_ylabel("std_y x_c(y)"); ax[1, 1].set_title("(e) 2D: delik merkez hattinin kink genligi")
    for a in ax.ravel()[:5]:
        a.set_xlabel("t omega_p")
        a.legend(fontsize=6)
    # (f) olay sayilari (doygunluk + 5 sonrasi, esik 0.3)
    a = ax[1, 2]
    nms = [n for n, _ in groups]
    xx = np.arange(len(nms))
    sp = [A[n]["th30"]["splits_persistent"] for n in nms]
    sr = [A[n]["th30"]["splits"] for n in nms]
    mg = [A[n]["th30"]["merges"] for n in nms]
    dt_ = [A[n]["th30"]["deaths"] for n in nms]
    a.bar(xx - 0.3, mg, 0.2, label="birlesme")
    a.bar(xx - 0.1, sr, 0.2, label="bolunme (ham)")
    a.bar(xx + 0.1, sp, 0.2, label="bolunme (kalici, 5/omega_p)")
    a.bar(xx + 0.3, dt_, 0.2, label="kaybolma")
    a.set_xticks(xx); a.set_xticklabels(nms, rotation=60, fontsize=6, ha="right")
    a.set_title("(f) olaylar (esik 0.3 phi_ref)")
    a.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_dimension.png"))
    plt.close(fig)


def fig_3d():
    rows = [n for n in ("3d_L20_seed", "3d_L20_quietkick", "3d_L20_seed_B2") if ld(n) is not None]
    if not rows:
        return
    times = (15, 30, 50, 80, 120)
    fig, ax = plt.subplots(2 * len(rows), len(times), figsize=(1.8 * len(times), 3.6 * len(rows)), squeeze=False)
    for i, nm in enumerate(rows):
        z = ld(nm)
        L = z["L"]
        for j, tt in enumerate(times):
            k = int(np.argmin(np.abs(z["t"] - tt)))
            ph = z["phi"][k]
            if abs(z["t"][k] - tt) > 1:
                ax[2 * i, j].axis("off"); ax[2 * i + 1, j].axis("off"); continue
            a = ax[2 * i, j]
            a.imshow(ph[:, :, 0].T, origin="lower", cmap="RdBu_r", vmin=-2, vmax=2, extent=[0, L[0], 0, L[1]])
            a.set_xticks([]); a.set_yticks([])
            if j == 0:
                a.set_ylabel(f"{nm}\nx-y kesit z=0", fontsize=6)
            a.set_title(f"t={tt}, max phi={ph.max():.2f}", fontsize=7)
            a = ax[2 * i + 1, j]
            ix = int(np.argmax(ph.max(axis=(1, 2))))
            a.imshow(ph[ix].T, origin="lower", cmap="RdBu_r", vmin=-2, vmax=2, extent=[0, L[1], 0, L[2]])
            a.set_xticks([]); a.set_yticks([])
            if j == 0:
                a.set_ylabel("y-z kesit (en derin x)", fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_3d.png"))
    plt.close(fig)


def fig_phase_space():
    names = [n for n in ("1d_L40_seed", "2d_L40_narrow", "2d_L40_wide", "2d_L40_wide_B2", "2d_L40_wide_quietkick",
                         "2d_L40_wide_quietkick_B2", "1d_L20_seed", "2d_L20_seed", "3d_L20_seed",
                         "3d_L20_quietkick") if ld(n) is not None]
    if not names:
        return
    fig, ax = plt.subplots(2, len(names), figsize=(2.0 * len(names), 4.2), squeeze=False)
    for j, nm in enumerate(names):
        z = ld(nm)
        L = float(z["L"][0])
        for i, tsel in enumerate((20, 150)):
            k = int(np.argmin(np.abs(z["ps_t"] - tsel)))
            a = ax[i, j]
            a.imshow(z["ps_h"][k].T, origin="lower", aspect="auto", extent=[0, L, -2.5, 2.5], cmap="viridis")
            a.set_title(f"{nm}\nt={z['ps_t'][k]:.0f}", fontsize=6)
            a.set_xticks([]);
            if j > 0:
                a.set_yticks([])
    ax[0, 0].set_ylabel("v_x"); ax[1, 0].set_ylabel("v_x")
    fig.suptitle("x-v_x faz uzayi (dik koordinatlar uzerinden toplanmis)", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_phase_space.png"))
    plt.close(fig)


if __name__ == "__main__":
    which = sys.argv[1:] or ["validation", "1d", "2d", "dimension", "3d", "phase"]
    for w in which:
        dict(validation=fig_validation, **{"1d": fig_1d, "2d": fig_2d_maps, "3d": fig_3d},
             dimension=fig_dimension, phase=fig_phase_space)[w]()
        print("ok", w, flush=True)
