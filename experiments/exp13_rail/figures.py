"""Exp13 figure (fig_exp13.png): robustness of the single filament, end-fed and distributed-ballast inhibitors
(kymographs of V_b(y) vs V_DD), inhibitor profile at the split, convergence of the R_sq = 1e10 case."""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))
import spatial as sp  # noqa: E402

OUT = ROOT / "results" / "exp13"
BASE = dict(ss=0.02, vs=0.005, fii=10, vov0=-0.1, g_sub=1e-13, von=1.05)


def kymo(kw, vdd=(0.6, 3.4), dv=0.02):
    vlist = list(np.round(np.arange(vdd[0], vdd[1] + 1e-9, dv), 4))
    r = sp.simulate(vdd_up=vdd, dvdd=dv, store_at=[("up", v) for v in vlist] + [("down", v) for v in vlist], **BASE, **kw)
    up = np.array([r["profiles"][f"up_{v:.2f}"] for v in vlist])
    dn = np.array([r["profiles"][f"down_{v:.2f}"] for v in vlist])
    vd = np.array([r["profiles"][f"up_{v:.2f}_VD"] for v in vlist])
    return r, np.array(vlist), up, dn, vd


fig, axs = plt.subplots(2, 3, figsize=(17, 10))

# (a) robustness of the single filament (global coupling, 5 parameter sets)
ax = axs[0, 0]
for p in sorted(OUT.glob("A*_global.json")):
    d = json.loads(p.read_text(encoding="utf-8"))
    h = [x for x in d["result"]["hist"] if x["dir"] == "up"]
    s = d["job"]["set"]
    ax.plot([x["vdd"] for x in h], [x["frac_high"] * 100 for x in h], lw=1.2,
            label=f"SS {s['ss']*1e3:.0f}, V_s {s['vs']*1e3:.0f} mV, F_ii {s['fii']}, V_ov0 {s['vov0']}, g {s['g_sub']:.0e}")
ax.set_xlabel("V_DD (V), R_L = 5 kΩ")
ax.set_ylabel("iletken (V_b > 0,5 V) genişlik %")
ax.set_title("(a) global bağlaşım: 5 parametre setinin 5'inde tek filaman\n(hepsinde doğum 1, bölünme 0)", fontsize=10)
ax.legend(fontsize=6)

# (b) end-fed rail kymograph (W = 100, r_D = 300)
r_b, vl, up_b, dn_b, _ = kymo(dict(R_sq=1e12, W=100.0, dy=0.5, r_D=300.0, contacts="both"))
ax = axs[0, 1]
ax.imshow(np.vstack([up_b, dn_b[::-1]]).T, aspect="auto", origin="lower", cmap="magma", vmin=0, vmax=1.05,
          extent=[0, 2 * len(vl), 0, 100])
ax.axvline(len(vl), color="w", lw=1)
ax.set_xticks([0, len(vl) // 2, len(vl), 3 * len(vl) // 2, 2 * len(vl)],
              [f"{vl[0]:.1f}", f"{vl[len(vl)//2]:.1f}", f"{vl[-1]:.1f}", f"{vl[len(vl)//2]:.1f}", f"{vl[0]:.1f}"])
ax.set_xlabel("V_DD (V): yukarı | aşağı")
ax.set_ylabel("y (µm)")
ax.set_title("(b) uçlardan beslenen dirençli drenaj (r_D = 300 Ω/µm, W = 100 µm):\nfilamanlar kontaklarda doğuyor; aşağıda ortadan bölünme", fontsize=10)

# (c) distributed ballast kymograph (lambda = 40 um, R_L = 5 kOhm, W = 300)
r_c, vl_c, up_c, dn_c, vd_c = kymo(dict(R_sq=1e12, W=300.0, dy=0.5, r_D=10.0, contacts="distributed", r_v=10.0 * 40.0**2,
                                         R_L=5e3))
ax = axs[0, 2]
ax.imshow(up_c.T, aspect="auto", origin="lower", cmap="magma", vmin=0, vmax=1.05, extent=[vl_c[0], vl_c[-1], 0, 300])
ax.set_xlabel("V_DD (V), yukarı tarama")
ax.set_ylabel("y (µm)")
ax.set_title("(c) dağıtık balast (λ = 40 µm, R_L = 5 kΩ, W = 300 µm):\nbüyüme → ikiye bölünme (3/3 tohumda)", fontsize=10)

# (d) inhibitor profile around the split
ax = axs[1, 0]
k_split = next((i for i, h in enumerate([x for x in r_c["hist"] if x["dir"] == "up"]) if h["splits"]), None)
y = np.arange(up_c.shape[1]) * 0.5
for i, c in ((max(k_split - 3, 0), "tab:blue"), (k_split, "tab:red")) if k_split else ():
    ax.plot(y, up_c[i], c=c, label=f"V_b, V_DD = {vl_c[i]:.2f} V")
ax3 = ax.twinx()
if k_split:
    ax3.plot(y, vd_c[max(k_split - 3, 0)], c="tab:blue", ls=":", lw=1)
    ax3.plot(y, vd_c[k_split], c="tab:red", ls=":", lw=1)
ax3.set_ylabel("yerel drenaj gerilimi V_D(y) (V, noktalı)")
ax.set_xlabel("y (µm)")
ax.set_ylabel("V_b (V)")
ax.set_title("(d) bölünme anı: filamanın ortasında V_D en düşük\n(darbe iyonlaşması orada en çok bastırılıyor)", fontsize=10)
ax.legend(fontsize=7, loc="lower left")

# (e) convergence of R_sq = 1e10
ax = axs[1, 1]
conv = json.loads((ROOT / "results" / "exp12" / "spatial_convergence.json").read_text(encoding="utf-8"))
for dy in (0.5, 0.25, 0.1):
    k = f"Rsq1e+10_dy{dy}_nsub12"
    if k in conv:
        vs_ = sorted(float(v) for v in conv[k]["up"])
        ax.plot(vs_, [conv[k]["up"][str(v)][1] * 100 for v in vs_], "o-", label=f"dy = {dy} µm (exp12)")
for p in sorted(OUT.glob("D_Rsq1e10_dy*.json")):
    d = json.loads(p.read_text(encoding="utf-8"))
    h = [x for x in d["result"]["hist"] if x["dir"] == "up" and 1.2 <= x["vdd"] <= 2.2]
    ax.plot([x["vdd"] for x in h], [x["frac_high"] * 100 for x in h], lw=1.5, label=f"dy = {d['job']['kw']['dy']} µm (exp13)")
ax.set_xlabel("V_DD (V)")
ax.set_ylabel("iletken genişlik %")
ax.set_title("(e) R_□ = 10¹⁰ yakınsaması: dy ≤ 0,1 µm'de aynı\n(kaba 0,5 µm'deki 'tümü birden' artefakttı)", fontsize=10)
ax.legend(fontsize=7)

# (f) event summary
ax = axs[1, 2]
rows = []
for p in sorted(OUT.glob("[ABCD]*.json")):
    if p.stem.endswith("_fine"):
        continue
    d = json.loads(p.read_text(encoding="utf-8"))
    ev = d["summary"]["events"]["up"]
    rows.append((p.stem, d["summary"]["n_max_up"], ev["births"], ev["splits"]))
ax.axis("off")
txt = "koşu                          N_max  doğum  bölünme(yukarı)\n" + "\n".join(
    f"{n:30s} {m:5d} {b:6d} {s:8d}" for n, m, b, s in rows)
ax.text(0.0, 1.0, txt, family="monospace", fontsize=7, va="top")
ax.set_title("(f) tüm koşular, dV_DD = 0,04\n(C_lam5_RL0'daki 3 filaman ince adımda kayboldu)", fontsize=9)
fig.suptitle("Exp13: yerel inhibitörlü (dirençli drenaj) SKY130 nFET, 4 K eşik altı kilitlenme — çoklu filaman ve bölünme", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "fig_exp13.png", dpi=100)
print("ok", "bölünme adımı:", k_split, vl_c[k_split] if k_split else None)
