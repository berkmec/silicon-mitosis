"""Exp14 figure (results/exp14_dimension/fig_exp14.png) and the outcome table outcomes.json collected from all
axisymmetric runs (results/exp14_dimension/axi/*.json): first event = division (>= 2 bodies) or hollowing (enclosed dilute
region = liquid shell), or none."""
import glob
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active_axi import AxiCH  # noqa: E402

OUT = ROOT / "results" / "exp14_dimension"
def dnu_of(nmk):
    """Chemical drive of a Z_nmk* run: dnu = (nu-/k + eps - 1) k with eps = 0.2, k = 0.01 (nmk = nu-/k)."""
    return round((nmk + 0.2 - 1.0) * 0.01, 6)


def first_event(path):
    r = json.loads(Path(path).read_text(encoding="utf-8"))
    j, h = r["job"], r["hist"]
    ts, th = r.get("t_split"), r.get("t_hollow")
    if "t_split" not in r:                                   # runs made before the event bookkeeping was added
        ts = next((x["t"] for x in h if x["n_axis"] + x["n_ring"] >= 2), None)
        z = np.load(str(path).replace(".json", ".npz"))
        m = AxiCH(j["d"], j["nz"], j["nr"], h=j["h"], dt=j["dt"])
        for k in sorted((k for k in z.files if k.startswith("t")), key=lambda s: int(s[1:])):
            if m.bubbles(z[k].astype(float)):
                th = float(k[1:])
                break
    kind = None
    if ts is not None or th is not None:
        kind = "split" if (th is None or (ts is not None and ts <= th)) else "hollow"
    t_ev = ts if kind == "split" else th
    n_pieces = None
    if kind == "split":
        x = next(x for x in h if x["t"] >= ts)
        n_pieces = len(x["bodies"])
    R_ev = None
    if t_ev is not None:
        x = min(h, key=lambda x: abs(x["t"] - t_ev))
        R_ev = max(b["R_eq"] for b in x["bodies"]) if kind == "hollow" else sum(b["R_eq"] ** j["d"] for b in x["bodies"]) ** (1 / j["d"])
    return dict(name=r["job"]["name"], d=j["d"], nmk=j["chem"]["nmk"], dnu=dnu_of(j["chem"]["nmk"]), kind=kind,
                t=t_ev, n=n_pieces, R=R_ev, T_end=h[-1]["t"])


def main():
    rows = [first_event(p) for p in sorted(glob.glob(str(OUT / "axi" / "Z_*.json")))]
    (OUT / "outcomes.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    for r in sorted(rows, key=lambda r: (r["nmk"], r["d"])):
        print(r)

    fig = plt.figure(figsize=(18, 11))
    gs = fig.add_gridspec(2, 3)
    # (a) passive Rayleigh-Plateau
    ax = fig.add_subplot(gs[0, 0])
    prp = json.loads((OUT / "passive_rp.json").read_text(encoding="utf-8"))
    for (d, v), col in zip(prp.items(), ("C0", "C1", "C2", "C3")):
        ax.plot(v["kR"], np.array(v["sigma"]) * 1e4, "o-", c=col, label=f"d={float(d):g}: ölçülen kR_c={v['kR_c']:.3f}")
        ax.axvline(v["theory"], c=col, ls="--", lw=1)
    ax.axhline(0, c="k", lw=0.8)
    ax.set_xlabel("k R (tüp, R = 8 w)")
    ax.set_ylabel("büyüme hızı σ (10⁻⁴ / t₀)")
    ax.set_title("(a) pasif Cahn–Hilliard: Rayleigh–Plateau bandı kR < √(d−2)\n(kesikli çizgi: kuram; 2B'de bant yok)", fontsize=10)
    ax.legend(fontsize=8)
    # (b) gallery: first event, dnu = -0.001
    sub = gs[0, 1].subgridspec(2, 3, wspace=0.05, hspace=0.25)
    pick = [("Z_nmk0.7_d2", None), ("Z_nmk0.7_d2.4", None), ("Z_nmk0.7_d3", None), ("Z_nmk0.7_d3.5", None),
            ("Z_nmk0.7_d4", None), ("Z_nmk0.7_d6", None)]
    byname = {r["name"]: r for r in rows}
    for i, (name, _) in enumerate(pick):
        a = fig.add_subplot(sub[i // 3, i % 3])
        r = byname.get(name)
        z = np.load(OUT / "axi" / f"{name}.npz")
        keys = [k for k in z.files if k.startswith(("split", "hollow"))]
        if not keys:
            tt = [int(k[1:]) for k in z.files if k.startswith("t")]
            target = r["t"] if r and r["t"] else max(tt)
            keys = [f"t{min(tt, key=lambda t: (t < target, abs(t - target)))}"]
        c = z[keys[0]].astype(float)
        full = np.concatenate([c[:, ::-1], c], axis=1)
        nz = full.shape[0]
        win = min(nz, int(nz * 0.8))
        a.imshow(full[(nz - win) // 2:(nz + win) // 2].T, origin="lower", cmap="viridis", vmin=0, vmax=1, aspect="equal")
        lab = {"split": f"bölünme ({r['n']} parça)", "hollow": "kabuk", None: "şerit, olay yok"}[r["kind"] if r else None]
        a.set_title(f"d = {name.split('_d')[1]}: {lab}\n{keys[0]}", fontsize=8)
        a.axis("off")
    fig.text(0.515, 0.925, "(b) Δν = −0,001, O(d−1) simetrik indirgeme, ilk olay (z yatay, ±r dikey)", fontsize=10, ha="center")
    # (c) phase diagram
    ax = fig.add_subplot(gs[0, 2])
    mk = {"split": ("o", "C2", "bölünme"), "hollow": ("s", "C3", "kabuk (içi boş)"), None: ("x", "0.4", "olay yok (şerit/damla)")}
    seen = set()
    for r in rows:
        m_, c_, lab = mk[r["kind"]]
        ax.scatter(r["d"], r["dnu"] * 1e3, marker=m_, c=c_, s=70, label=None if lab in seen else lab)
        seen.add(lab)
        if r["kind"] == "split":
            ax.annotate(str(r["n"]), (r["d"], r["dnu"] * 1e3), textcoords="offset points", xytext=(5, 4), fontsize=7)
    th = json.loads((OUT / "theory.json").read_text(encoding="utf-8"))
    for dnu, dc in th.get("d_c", {}).items():
        if isinstance(dc, (int, float)):
            ax.plot([dc], [float(dnu) * 1e3], "k|", ms=18, mew=2)
    ax.plot([], [], "k|", ms=12, mew=2, label="sonsuz tüp boğulma eşiği (doğrusal olmayan kuram)")
    ax.set_xlabel("boyut d (sürekli)")
    ax.set_ylabel("Δν (10⁻³ / t₀)")
    ax.set_title("(c) bölünme penceresi Δν ile yukarı kayıyor: Δν = −0,001'de 2,3 ≤ d ≤ 3,5;\n"
                 "6B ancak güçlü iç tüketimde (Δν ≥ +0,012) bölünüyor", fontsize=10)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0.3, 0.99))
    ax.grid(alpha=0.3)
    # (d) mechanism of the upper bound
    ax = fig.add_subplot(gs[1, 0])
    wt = json.loads((OUT / "window_theory.json").read_text(encoding="utf-8"))
    ds = [v["d"] for v in wt.values() if v["dnu"] == -0.001]
    rh = [v["R_hollow"] for v in wt.values() if v["dnu"] == -0.001]
    ax.plot(ds, rh, "k-o", label="R_boşalma(d): yoğun çekirdeğin var olabildiği en büyük R (kuram)")
    ax.set_xlabel("d")
    ax.set_ylabel("yarıçap (w)")
    pf = [r for r in rows if r["nmk"] == 0.7 and r["kind"]]
    ax.scatter([r["d"] for r in pf if r["kind"] == "split"], [r["R"] for r in pf if r["kind"] == "split"], c="C2", marker="o",
               s=60, label="faz alanı: bölünme anındaki R_eş (toplam hacimden)")
    ax.scatter([r["d"] for r in pf if r["kind"] == "hollow"], [r["R"] for r in pf if r["kind"] == "hollow"], c="C3", marker="s",
               s=60, label="faz alanı: içi boşalma anındaki R_eş")
    axr = ax.twinx()
    for R, col in ((10.0, "C1"), (12.0, "C4")):
        rel = [next((x.get("rel") for x in v["rows"] if x["R"] == R and not x.get("hollow")), np.nan) for v in wt.values()
               if v["dnu"] == -0.001]
        axr.plot(ds, np.array(rel, dtype=float) * 1e3, "--", c=col, label=f"σ₂ − v/R, R = {R:g} (sağ eksen)")
    axr.axhline(0, c="0.6", lw=0.8)
    axr.set_ylabel("göreli uzama hızı σ₂ − v/R (10⁻³ / t₀)")
    ax.set_title("(d) üst sınırın mekanizması (Δν = −0,001): boyut artınca uzama zayıflıyor,\n"
                 "damla bölünemeden boşalma yarıçapına ulaşıp kabuğa dönüşüyor", fontsize=10)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="lower right")
    ax.set_ylim(6, 22)
    # (e) 3D validation against the full Cartesian runs of exp11
    ax = fig.add_subplot(gs[1, 1])
    for tag, name, col, lab in (("b3d", "Z_nmk0.7_d3", "C0", "Δν = −0,001"), ("b3d_dnum0p004", "Z_nmk0.4_d3", "C3", "Δν = −0,004")):
        b = json.loads((ROOT / "results" / "exp11" / f"{tag}.json").read_text(encoding="utf-8"))
        ax.plot([x["t"] for x in b["hist"]], [max(x["radii"]) if x["radii"] else np.nan for x in b["hist"]], c=col, lw=2,
                label=f"Kartezyen 3B 64³ (exp11), {lab}: en büyük R_eş")
        a = json.loads((OUT / "axi" / f"{name}.json").read_text(encoding="utf-8"))
        ax.plot([x["t"] for x in a["hist"]], [max(bb["R_eq"] for bb in x["bodies"]) if x["bodies"] else np.nan for x in a["hist"]],
                c=col, ls="--", label=f"eksenel d = 3, {lab}")
    ax.axvspan(1750, 2000, color="C0", alpha=0.12)
    ax.text(1760, 18.8, "Kartezyen\nbölünme", fontsize=7, color="C0")
    ax.set_xlim(0, 4000)
    ax.set_xlabel("t (t₀)")
    ax.set_ylabel("R_eş (w)")
    ax.set_title("(e) 3B doğrulama: ilk bölünme t = 1850 vs 1750–2000 (kızlar 9,7 vs 9,9);\n"
                 "Δν = −0,004'te ikisi de kabuk (Kartezyen merkez c: 0,78 → 0,03)", fontsize=10)
    ax.legend(fontsize=7)
    # (f) 6D hypercube vs O(5) reduction
    ax = fig.add_subplot(gs[1, 2])
    h6 = json.loads((OUT / "hyper6d.json").read_text(encoding="utf-8"))
    ax.plot([x["t"] for x in h6["cart6d"]], [x["R_eq"] for x in h6["cart6d"]], "ko-", label="gerçek 6B hiperküp 16⁶ (16,8 M nokta), dx = 1")
    for k, st, lab in (("axi_big_h1.0", "C0--", "O(5) indirgeme, h = 1"), ("axi_big_h0.5", "C1:", "O(5) indirgeme, h = 0,5")):
        ax.plot([x["t"] for x in h6[k]], [x["R_eq"] for x in h6[k]], st, lw=2, label=lab)
    ax.set_xlabel("t (t₀)")
    ax.set_ylabel("R_eş (w)")
    ax.set_title("(f) 6B hiperküpte Cahn–Hilliard (aktif damla büyümesi)\nO(5) indirgeme ile her noktada fark ≤ %2,5", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    ins = ax.inset_axes([0.06, 0.55, 0.36, 0.36])
    ins.imshow(np.load(OUT / "hyper6d_slice.npy").astype(float).T, origin="lower", cmap="viridis", vmin=0, vmax=1)
    ins.set_xticks([])
    ins.set_yticks([])
    ins.set_title("6B alanın 2B kesiti, t = 40", fontsize=7)
    fig.suptitle("Exp14: aktif Cahn–Hilliard damlasında bölünme boyuta bağlı — bölünme yalnız bir boyut penceresinde", fontsize=13,
                 y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "fig_exp14.png", dpi=95)
    print("ok")


if __name__ == "__main__":
    main()
