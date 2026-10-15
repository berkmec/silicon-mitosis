"""Exp11 comparison: which drives cell division harder, the supercapacitor idea (A) or dissipative structures (B)?
Reads results/exp11/runs/*.json (+ npz), model_a_phys.json, phys.json; writes classification.json, b_rho.json,
comparison.csv, fig_exp11_sim.png, fig_exp11_phys.png."""
import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "exp11"
RUNS = RES / "runs"
SPINODAL = (3 - np.sqrt(3)) / 6
EPS = (0.08, 0.11, 0.14, 0.17, 0.20)
DNU = (-0.006, -0.004, -0.002, 0.0, 0.002, 0.004)
A_COLS = [(0.0, "inf")] + [(A, d) for A in (0.001, 0.003, 0.01, 0.02, 0.08) for d in ("inf", "3", "1")]


def a_name(drive, A, d):
    """Run name of a supercapacitor (A) run; weak-Coulomb runs carry the A2_ prefix."""
    if drive == "leaky" and A == 0.0:
        return "B_e0.17_d+0.000"  # leaky charging without Coulomb = the dnu = 0 point of the B scan
    pre = "A2_" if A in (0.001, 0.003, 0.01) else "A_"
    return f"{pre}{drive}_A{A:g}_d{d}"


def load(name):
    p = RUNS / f"{name}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def doubling_time(ts, ns):
    """Fit log2 N = a + t / T_d over the growth phase 2 <= N <= max(N)/2 (needs >= 3 points and a factor >= 2 of
    growth inside the window, otherwise None: flat counts gave nonsense values of 1e18-1e19 before the review)."""
    ts, ns = np.asarray(ts, float), np.asarray(ns, float)
    nmax = ns.max()
    sel = (ns >= 2) & (ns <= max(2.0, nmax / 2))
    if sel.sum() < 3 or ns[sel].max() < 2 * ns[sel].min():
        return None
    slope = np.polyfit(ts[sel], np.log2(ns[sel]), 1)[0]
    return float(1 / slope) if slope > 1e-6 else None


def shape_ratio(name):
    """R_eq / r_inscribed of the largest droplet in the last snapshot (disk ~1; ellipse 2:1 -> 1.41; worms >> 1).
    Covariance aspect fails for meandering worms (review finding), this does not. None if no snapshot/droplet."""
    from scipy import ndimage

    p = RUNS / f"{name}.npz"
    if not p.exists():
        return None
    z = np.load(p)
    k = max(z.files, key=lambda s: int(s[1:]))
    c = z[k].astype(float) > 0.5
    if not c.any():
        return None
    n = c.shape[0]
    big = np.tile(c, (3, 3))
    lab, _ = ndimage.label(big)
    dist = ndimage.distance_transform_edt(big)
    sub = lab[n:2 * n, n:2 * n]
    ids, counts = np.unique(sub[sub > 0], return_counts=True)
    if len(ids) == 0:
        return None
    best = ids[np.argmax(counts)]
    comp = lab == best
    box = ndimage.find_objects(comp.astype(int))[0]  # tuple of slices of the component's bounding box
    if any(sl.stop - sl.start >= 2 * n for sl in box):
        return float("inf")  # wraps the periodic box: stripe / labyrinth
    return float(np.sqrt(comp.sum() / np.pi) / max(dist[comp].max(), 1e-9))


def classify_ch(r):
    R0 = float(r.get("spec", {}).get("R0", 8.0))
    # every run starts from ONE droplet of radius R0 (virtual record at t = 0)
    h = [dict(t=0.0, n=1, rmax=R0, rmean=R0, perc=False, cfar=0.0, aspect=1.0)] + r["hist"]
    ns = [x["n"] for x in h]
    out = dict(n_max=max(ns), n_end=ns[-1], t_end=h[-1]["t"], perc_end=h[-1]["perc"], t_div=None, R_div=None,
               spinodal_birth=False, T_double=None, n_events=0, div_type=None, stop=r.get("stop"),
               aspect_max=max(x.get("aspect", 1.0) for x in h), shape_ratio=r.get("_shape_ratio"))
    for i in range(1, len(h)):
        cf0, cf1 = h[i - 1]["cfar"], h[i]["cfar"]
        jump = h[i]["n"] - h[i - 1]["n"]
        # birth from the medium (not division) if the far field was/is beyond the spinodal, if no far field is left
        # (NaN: droplets everywhere), or if the count more than triples in one record (each parent gives at most 3
        # daughters; synchronous divisions such as 4 -> 8 are still allowed)
        spin = (cf0 == cf0 and cf0 >= SPINODAL) or (cf1 == cf1 and cf1 >= SPINODAL) or cf1 != cf1 \
            or h[i]["n"] > 3 * max(h[i - 1]["n"], 1)
        if jump > 0:
            if spin:
                out["spinodal_birth"] = True
            elif not h[i - 1]["perc"] and h[i - 1]["n"] >= 1:
                out["n_events"] += 1
        if out["t_div"] is None and h[i - 1]["n"] == 1 and not h[i - 1]["perc"] and h[i]["n"] >= 2 and not spin:
            out["t_div"], out["R_div"] = h[i]["t"], h[i - 1]["rmax"]
            grew = max(x["rmax"] for x in h[:i]) > 1.1 * R0
            out["aspect_before_div"] = max(x.get("aspect", 1.0) for x in h[:i])
            if not grew:
                out["div_type"] = "büyük başlangıç damlası bölündü"
            elif out["aspect_before_div"] >= 4.0:
                out["div_type"] = "solucan/şerit → kopma"  # elongated worm broke, not a cell-like binary fission
            else:
                out["div_type"] = "büyüme→ikili bölünme"
    if out["t_div"] is not None:
        cls = "bölündü"
        out["T_double"] = doubling_time([x["t"] for x in h], ns)
    elif out["spinodal_birth"]:
        cls = "spinodal doğum"
    elif out["perc_end"]:
        cls = "süzüldü (şerit)"
    elif out["n_end"] == 0:
        cls = "çözündü"
    else:
        r_end = h[-1]["rmax"]
        r_prev = next((x["rmax"] for x in reversed(h) if x["t"] <= h[-1]["t"] - 2000), r_end)
        cls = "tek, kararlı" if abs(r_end - r_prev) <= 0.02 * max(r_end, 1e-9) else "tek, hâlâ değişiyor"
        sr = out.get("shape_ratio")
        if out["aspect_max"] >= 4.0 or (sr is not None and sr >= 1.6):
            cls = "solucan (uzadı, kopmadı)"
    out["class"] = cls
    out["persistent"] = bool(cls == "bölündü" and out["n_end"] >= 2)  # daughters still alive at the end
    return out


def classify_gs(r):
    h = r["hist"]
    ns = [x["n"] for x in h]
    last = h[-1]
    if last["area_frac"] == 0.0:
        cls = "söndü"
    elif last["area_frac"] > 0.95:
        cls = "düzgün"
    elif last["perc"]:
        cls = "şerit/labirent"  # review: a percolating pattern is a labyrinth even if some spots exist
    elif max(ns) >= 4 and last["n"] >= 4:
        cls = "çoğalan lekeler"
    else:
        cls = "az leke"
    return dict(n_max=max(ns), n_end=last["n"], cls_area=last["area_frac"], perc_end=last["perc"], **{"class": cls},
                T_double=doubling_time([x["t"] for x in h], ns))


def main():
    cl = {}
    for p in sorted(RUNS.glob("*.json")):
        if p.stem.endswith("__ext"):
            continue  # review extensions are reported separately
        r = json.loads(p.read_text(encoding="utf-8"))
        if not p.stem.startswith("GS_"):
            r["_shape_ratio"] = shape_ratio(p.stem)
        cl[p.stem] = classify_gs(r) if p.stem.startswith("GS_") else classify_ch(r)
    (RES / "classification.json").write_text(json.dumps(cl, indent=1, ensure_ascii=False), encoding="utf-8")

    # ------------------------------------------------ rho = R_div / l for dissipative runs (l = sqrt(D/k) = 10)
    rd = [v["R_div"] / 10.0 for k, v in cl.items()
          if k.startswith("B_e") and v["class"] == "bölündü" and v["div_type"] != "büyük başlangıç damlası bölündü"]
    if rd:
        (RES / "b_rho.json").write_text(json.dumps(dict(rho=float(np.median(rd)), values=rd), indent=1),
                                        encoding="utf-8")

    # ------------------------------------------------ table
    rows = []
    for k, v in cl.items():
        rows.append(dict(run=k, sinif=v["class"], bolunme_turu=v.get("div_type"), N_max=v["n_max"], N_son=v["n_end"],
                         bolunme_olayi=v.get("n_events"), t_bolunme=v.get("t_div"), R_bolunme=v.get("R_div"),
                         en_boy_max=v.get("aspect_max"), T_ikilenme=v.get("T_double"), durma=v.get("stop")))
    with open(RES / "comparison.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    def frac(prefix, ok, excl=None):
        sel = [v for k, v in cl.items() if k.startswith(prefix) and not (excl and k.startswith(excl))]
        return [sum(bool(ok(v)) for v in sel), len(sel)]

    pers = lambda v: v.get("persistent")  # noqa: E731
    stats = {
        "B 2B Zwicker: kalıcı bölünme": frac("B_e", pers),
        "B 2B Zwicker: büyüme→ikili bölünme": frac("B_e", lambda v: v.get("div_type") == "büyüme→ikili bölünme"),
        "B 2B Zwicker: solucan→kopma": frac("B_e", lambda v: v.get("div_type") == "solucan/şerit → kopma"),
        "B 2B Zwicker: büyük başlangıç bölündü (kalıcı)":
            frac("B_e", lambda v: v.get("div_type") == "büyük başlangıç damlası bölündü" and pers(v)),
        "B 2B Zwicker: solucan (uzadı, kopmadı)": frac("B_e", lambda v: v["class"] == "solucan (uzadı, kopmadı)"),
        "B 2B Zwicker: damla yaşadı": frac("B_e", lambda v: v["n_end"] >= 1),
        "A (A≥0,02): kalıcı bölünme": frac("A_", pers),
        "A (A≥0,02): hücre yaşadı": frac("A_", lambda v: v["n_end"] >= 1),
        "A2 (A≤0,01): kalıcı bölünme": frac("A2_", pers),
        "A2 (A≤0,01): hücre yaşadı": frac("A2_", lambda v: v["n_end"] >= 1),
        "A3 (şarj şiddeti J taraması): kalıcı bölünme": frac("A3_", pers),
        "A3 (şarj şiddeti J taraması): ortamdan doğum": frac("A3_", lambda v: v["class"] == "spinodal doğum"),
        "A+A2 spinodal doğum (bölünme değil)": frac("A", lambda v: v["class"] == "spinodal doğum", excl="A3_"),
        "GS: çoğalan lekeler": frac("GS_", lambda v: v["class"].startswith("çoğalan")),
    }
    (RES / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False), encoding="utf-8")
    for k, (a, b) in stats.items():
        print(f"{k:48s}: {a}/{b}")

    fig_sim(cl)
    fig_phys()
    fig_cells()


def fig_cells():
    """3D Zwicker division (mid-plane), Gray-Scott ramp (dividing while shrinking), Gray-Scott pilot."""
    panels = []
    p = RES / "b3d_slices.npz"
    if p.exists():
        z = np.load(p)
        for k in ("t1000", "t2000"):
            if k in z.files:
                panels.append((z[k].astype(float), f"B Zwicker 3B, orta kesit {k}", "viridis", (0, 1)))
    p = RES / "gs_ramp.npz"
    if p.exists():
        z = np.load(p)
        h = json.loads((RES / "gs_ramp.json").read_text(encoding="utf-8"))["hist"]
        for s in range(4):
            if f"stage{s}" in z.files:
                last = [x for x in h if x["stage"] == s][-1]
                panels.append((z[f"stage{s}"].astype(float), f"GS rampa evre {s}: D×{2.0**-s:g}, N={last['n']}",
                               "magma", (0, 0.5)))
    fig, axs = plt.subplots(1, len(panels), figsize=(3.3 * len(panels), 3.6))
    for ax, (img, title, cmap, lim) in zip(np.atleast_1d(axs), panels):
        ax.imshow(img.T, origin="lower", cmap=cmap, vmin=lim[0], vmax=lim[1])
        ax.set_title(title, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("Gerçek hücre benzeri bölünme: 3B aktif damla (ikiye) ve Gray–Scott (sürüş arttıkça küçülerek çoğalma)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(RES / "fig_exp11_cells.png", dpi=100)
    plt.close(fig)


MARK = {"büyüme→ikili bölünme": ("*", "tab:red"), "solucan/şerit → kopma": ("P", "tab:red"),
        "büyük başlangıç damlası bölündü": ("v", "salmon"), "solucan (uzadı, kopmadı)": ("D", "tab:brown"),
        "spinodal doğum": ("^", "tab:orange"), "süzüldü (şerit)": ("s", "tab:purple"),
        "çözündü": ("x", "0.5"), "tek, kararlı": ("o", "tab:blue"), "tek, hâlâ değişiyor": ("o", "tab:cyan")}


def mk(v):
    return MARK.get(v.get("div_type") or v["class"], ("?", "k"))


def snap(name, key=None):
    p = RUNS / f"{name}.npz"
    if not p.exists():
        return None, None
    z = np.load(p)
    keys = sorted(z.files, key=lambda s: int(s[1:]))
    k = keys[-1] if key is None else min(keys, key=lambda s: abs(int(s[1:]) - key))
    return z[k].astype(float), k


def fig_sim(cl):
    fig = plt.figure(figsize=(17, 10.5))
    gsp = fig.add_gridspec(2, 4, hspace=0.38, wspace=0.32)
    # (a) B phase diagram
    ax = fig.add_subplot(gsp[0, 0:2])
    for e in EPS:
        for d in DNU:
            v = cl.get(f"B_e{e:.2f}_d{d:+.3f}")
            if not v:
                continue
            m, c = mk(v)
            ax.scatter(e, d * 1e3, marker=m, c=c, s=160 if m == "*" else 70)
            ax.annotate(str(v["n_max"]), (e, d * 1e3), textcoords="offset points", xytext=(7, 5), fontsize=8)
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axvline(SPINODAL, color="tab:orange", lw=0.8, ls=":")
    for lab, (m, c) in MARK.items():
        ax.scatter([], [], marker=m, c=c, label=lab)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False)
    ax.set_xlim(0.07, 0.215)
    ax.set_xlabel("aşırı doygunluk ε = ν₊/k   (noktalı: spinodal 0,211)")
    ax.set_ylabel("Δν = ν₋+ν₊−k  (×10⁻³)")
    ax.set_title("B (disipatif, Zwicker 2B): tek damla, L=128w, T≤2·10⁴  (sayı: N_max)\n"
                 "kesikli çizgi Δν=0: sızıntılı şarj s=J−λc ≡ Ohta–Kawasaki (λ/q²)", fontsize=9)
    # (b) A grid
    ax = fig.add_subplot(gsp[0, 2:4])
    rowsA = [("denge", "eq"), ("şarj J", "charge"), ("sızıntılı J−λc", "leaky")]
    for i, (lab, drive) in enumerate(rowsA):
        for j, (A, d) in enumerate(A_COLS):
            v = cl.get(a_name(drive, A, d))
            if not v:
                ax.text(j, i, "·", ha="center", va="center", color="0.7")
                continue
            m, c = mk(v)
            ax.scatter(j, i, marker=m, c=c, s=200 if m == "*" else 80)
            ax.text(j, i - 0.32, f"{v['n_max']}", ha="center", fontsize=7)
    ax.set_xticks(range(len(A_COLS)), [f"{A:g}\n{d}" for A, d in A_COLS], fontsize=7)
    ax.set_xlabel("Coulomb A (üst) / kapı uzaklığı d (alt, w biriminde)", fontsize=8)
    ax.set_yticks(range(3), [r[0] for r in rowsA])
    ax.set_ylim(2.6, -0.6)
    ax.set_title("A (süperkapasitör: kapı-ekranlı Coulomb + elektrikle şarj)", fontsize=10)
    # (c) GS map
    ax = fig.add_subplot(gsp[1, 0])
    colors = {"çoğalan lekeler": "tab:red", "çoğalan lekeler (+şerit)": "salmon", "şerit/labirent": "tab:purple",
              "söndü": "0.6", "düzgün": "tab:green", "az leke": "tab:blue"}
    for k, v in cl.items():
        if not k.startswith("GS_"):
            continue
        F, kk = float(k.split("_F")[1].split("_k")[0]), float(k.split("_k")[1])
        ax.scatter(kk, F, c=colors.get(v["class"], "k"), s=60)
    for lab, c in colors.items():
        ax.scatter([], [], c=c, label=lab)
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False)
    ax.set_xlabel("k")
    ax.set_ylabel("F")
    ax.set_title("Gray–Scott (Pearson 1993), tek tohum", fontsize=10)
    # (d) N(t)
    ax = fig.add_subplot(gsp[1, 1])
    best_b = "B_e0.20_d-0.002"
    curves = [("GS_F0.026_k0.059", "B: Gray–Scott F=.026 k=.059"), (best_b, "B: Zwicker 2B ε=.20 Δν=−.002 (solucan)"),
              ("A2_charge_A0.003_dinf", "A: şarj, A=.003, d=∞"), ("A2_charge_A0.003_d1", "A: şarj, A=.003, d=1"),
              ("A_charge_A0.08_dinf", "A: şarj, A=.08, d=∞ (çözündü)")]
    for name, lab in curves:
        r = load(name)
        if r is None:
            continue
        t = np.array([x["t"] for x in r["hist"]], float)
        n = np.array([x["n"] for x in r["hist"]], float)
        ax.plot(t / t.max(), np.maximum(n, 0.5), label=lab)
    ax.set_yscale("log")
    ax.set_xlabel("t / t_son  (GS: 2·10⁴ adım; CH: 2·10⁴ t₀)")
    ax.set_ylabel("hücre (damla/leke) sayısı")
    ax.legend(fontsize=7)
    ax.set_title("Çoğalma eğrileri", fontsize=10)
    # (e)-(f) snapshots
    for j, (name, title) in enumerate([(best_b, "B Zwicker 2B: solucan koptu"),
                                       ("A2_charge_A0.003_dinf", "A şarjlı hücre, A=.003 d=∞")]):
        ax = fig.add_subplot(gsp[1, 2 + j])
        c, k = snap(name) if name else (None, None)
        if c is not None:
            ax.imshow(c.T, origin="lower", cmap="viridis", vmin=0, vmax=1)
            ax.set_title(f"{title}  [{k}]", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("Exp11: hangisi daha güçlü bölüyor? A = süperkapasitör fikri, B = disipatif yapı", fontsize=13)
    fig.savefig(RES / "fig_exp11_sim.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def fig_phys():
    a = json.loads((RES / "model_a_phys.json").read_text(encoding="utf-8"))
    ph = json.loads((RES / "phys.json").read_text(encoding="utf-8"))
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.2))
    ax = axs[0]
    for key, cur in a["curves"].items():
        n = np.array(cur["n_cm2"])
        c = np.array([np.nan if x is None else x for x in cur["C_star_uF_cm2"]])
        if np.all(np.isnan(c)):
            continue
        ax.plot(n, c, marker=".", label=key)
    ax.axhline(a["C_SKY130_uF_cm2"], color="k", ls="--", lw=1)
    ax.text(1.1e10, a["C_SKY130_uF_cm2"] * 1.1, "SKY130 C_ox = 0,83 µF/cm²", fontsize=8)
    ax.axhspan(1, 10, color="tab:green", alpha=0.12)
    ax.text(1.1e10, 3, "EDL / süperkapasitör kapısı 1–10 µF/cm²", fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e10, 1e13)
    ax.set_ylim(0.3, 3000)
    ax.set_xlabel("n_s (cm⁻²)")
    ax.set_ylabel("eşik kapasitans C* (µF/cm²)")
    ax.set_title("A fiziksel: düzgün Si elektron tabakası hangi C'de bölünür?\n(HF eğrileri hiç yok = hiçbir C'de bölünmez)",
                 fontsize=10)
    ax.legend(fontsize=7)
    ax = axs[1]
    qt = ph["quantum_target"]["4K"]
    for name, rows in ph["B"].items():
        for D in (0.1, 1.0, 10.0):
            rr = [r for r in rows if r["D_cm2_s"] == D]
            ax.plot([r["Pmax_W_cm2"] for r in rr], [r["R_min_nm"] for r in rr], marker="o",
                    ls="-" if "e-h" in name else "--", label=f"{name.split(' (')[0]}, D={D:g} cm²/s")
    for name, rows in ph.get("B_2D_layer", {}).items():
        if "|" in name:
            continue
        for D in (0.1, 1.0):
            rr = [r for r in rows if r["D_cm2_s"] == D]
            ax.plot([r["Pmax_W_cm2"] for r in rr], [r["R_min_nm"] for r in rr], marker="s", ls=":", lw=2,
                    label=f"{name.split(' (')[0]}, D={D:g} cm²/s")
    ax.axhspan(6, 20, color="tab:purple", alpha=0.15)
    ax.text(0.11, 7.5, "4 K kuantum hedefi R ≈ 6–20 nm (geometriye bağlı)", fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("4 K'de uzaklaştırılabilen güç P_max (W/cm²)  [varsayım]")
    ax.set_ylabel("en küçük sürdürülebilir hücre yarıçapı (nm)")
    asm = ph["assumptions"]
    extra = (f"\n2B solucan (ρ={asm['rho_2D_worm']:.2g}) için güç ×{asm['power_factor_2D_worm']:.1f} daha fazla"
             if asm.get("rho_2D_worm") else "")
    ax.set_title(f"B fiziksel: güç bütçesi + parçacık sayısı (N≥{asm['N_min']}), ρ=R_bölünme/ℓ={asm['rho']:.2g} (3B){extra}",
                 fontsize=9)
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(RES / "fig_exp11_phys.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
