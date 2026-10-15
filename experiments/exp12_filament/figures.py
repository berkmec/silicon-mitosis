"""Exp12 figure: fig_exp12.png (6 panels) from the BSIM4 tables, reduced-model maps and spatial runs."""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import local_bistability as lb  # noqa: E402
import reduced_model as rm  # noqa: E402
import spatial as sp  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp12"
KW = dict(vov0=-0.1, ss=0.02, fii=10, von=1.05, vs=0.005, g_sub=1e-13)

fig, axs = plt.subplots(2, 3, figsize=(17, 10))

# (a) BSIM4 feedback slope vs V_gs
ax = axs[0, 0]
vb = np.linspace(0, 1.0, 101)
vgs = np.round(np.arange(0.4, 1.81, 0.05), 3)
for vd, c in ((1.2, "tab:blue"), (1.8, "tab:red")):
    s = [np.max(np.gradient(np.log(lb.j_ii(v, vd, vb)), vb)) for v in vgs]
    ax.plot(vgs, s, "o-", c=c, ms=3, label=f"V_ds = {vd} V")
ax.axvline(rm.AN["vth0"], color="k", ls=":", lw=1)
ax.text(rm.AN["vth0"] + 0.02, 6.5, "V_T(V_b=0)", fontsize=8)
ax.set_xlabel("V_gs (V)")
ax.set_ylabel("maks d ln I_ii / dV_b (1/V)")
ax.set_title("(a) BSIM4 (SKY130, 300 K): yüzer gövde geri beslemesi\neşik altında güçlü, eşik üstünde zayıf", fontsize=10)
ax.legend(fontsize=8)

# (b) bistability (kink + hysteresis) vs subthreshold slope
ax = axs[0, 1]
rb = json.loads((OUT / "reduced_bistability.json").read_text(encoding="utf-8"))
ks = list(rb["by_ss"].keys())
ax.bar(range(len(ks)), [rb["by_ss"][k]["hysteresis"] / rb["by_ss"][k]["n"] for k in ks], color="tab:purple")
for i, k in enumerate(ks):
    w = rb["by_ss"][k]["width_median_V"]
    ax.text(i, rb["by_ss"][k]["hysteresis"] / rb["by_ss"][k]["n"] + 0.02, f"{rb['by_ss'][k]['hysteresis']}/108\nmedyan {w:.2f} V",
            ha="center", fontsize=7)
ax.set_xticks(range(len(ks)), ks, fontsize=8)
ax.set_ylim(0, 1.15)
ax.set_ylabel("kink + histerezis görülen düğme oranı")
ax.set_title("(b) indirgenmiş 4 K modeli: eşik altı eğimine göre kilitlenme\n(V_ov0, V_s, F_ii, g_sub taraması)", fontsize=10)

# (c) single-device S-curve
ax = axs[0, 2]
vds = np.round(np.arange(0.6, 2.6001, 0.01), 3)
for v in vds[::2]:
    for r in rm.roots_at(v, **KW):
        ax.plot(v, r, "k.", ms=2)
V_up, V_down, up, dn = rm.sweep(vds, **KW)
ax.plot(vds, up, "tab:blue", lw=1.5, label=f"yukarı tarama (kink V_up = {V_up} V)")
ax.plot(vds, dn, "tab:red", lw=1.5, ls="--", label=f"aşağı tarama (dönüş V_down = {V_down} V)")
ax.set_xlabel("V_ds (V)")
ax.set_ylabel("gövde gerilimi V_b (V)")
ax.set_title("(c) tek cihaz (SS 20 mV/dek, V_s 5 mV, F_ii 10, V_ov0 −0,1 V):\nkararlı dallar ve histerezis", fontsize=10)
ax.legend(fontsize=8)

# (d)-(e) spatial
runs = {}
for R_sq in (1e8, 1e12):
    runs[R_sq] = sp.simulate(R_sq=R_sq, dy=0.25, mismatch=True, vdd_up=(0.6, 2.4), dvdd=0.02,
                             store_at=[("up", v) for v in (1.4, 1.6, 1.8, 2.0)])
ax = axs[1, 0]
y = np.arange(len(runs[1e12]["profiles"]["up_1.40"])) * 0.25
for k, c in zip(("up_1.40", "up_1.60", "up_1.80", "up_2.00"), ("tab:blue", "tab:green", "tab:orange", "tab:red")):
    ax.plot(y, runs[1e12]["profiles"][k], c=c, label=f"V_DD = {k.split('_')[1]} V")
ax.set_xlabel("kanal genişliği boyunca y (µm), W = 100 µm")
ax.set_ylabel("V_b (V)")
ax.set_title("(d) donmuş kuyu R_□ = 10¹² Ω/□: akım filamanı\nsürüş arttıkça bölünmüyor, genişliyor", fontsize=10)
ax.legend(fontsize=8)
ax = axs[1, 1]
for R_sq, c in ((1e8, "tab:gray"), (1e12, "tab:red")):
    h = runs[R_sq]["hist"]
    for d_, ls in (("up", "-"), ("down", "--")):
        hh = [x for x in h if x["dir"] == d_]
        ax.plot([x["I_uA"] for x in hh], [x["vds"] for x in hh], c=c, ls=ls, lw=1.2,
                label=f"R_□ = {R_sq:.0e} Ω/□, {('yukarı' if d_ == 'up' else 'aşağı')}")
vM = json.loads((OUT / "maxwell.json").read_text(encoding="utf-8"))["V_ds_maxwell"]
ax.axhline(vM, color="k", ls=":", lw=1)
ax.text(5, vM + 0.03, f"Maxwell (eşit alan) V_ds = {vM:.3f} V", fontsize=8)
ax.set_xlabel("toplam drenaj akımı (µA), R_L = 5 kΩ")
ax.set_ylabel("V_ds (V)")
ax.set_title("(e) yük doğrusuyla I–V: filamanda V_ds Maxwell geriliminde kenetli\n(dikey bölüm), güçlü bağlaşımda tüm kanal birlikte geçiyor", fontsize=10)
ax.legend(fontsize=7)

# (f) ngspice single device above threshold: no hysteresis
ax = axs[1, 2]
sk = json.loads((OUT / "single_kink.json").read_text(encoding="utf-8"))
SP = Path(__file__).resolve().parent / "spice"
for name, c in (("warm_300K_contacted", "tab:blue"), ("frozen_von1.0_fii1", "tab:red"), ("frozen_von1.0_fii10", "tab:orange")):
    for d_, ls in (("up", "-"), ("dn", "--")):
        dd = np.loadtxt(SP / f"{name}_{d_}.dat")
        ax.plot(dd[:, 0], dd[:, 1] * 1e3, c=c, ls=ls, lw=1.1, label=name if d_ == "up" else None)
ax.set_xlabel("V_ds (V)")
ax.set_ylabel("I_d (mA), W = 10 µm")
ax.set_title("(f) ngspice + SKY130 BSIM4, V_gs = 1,2 V (eşik üstü):\nyüzer gövdeyle bile histerezis yok (düz = yukarı, kesik = aşağı)", fontsize=10)
ax.legend(fontsize=7)
fig.suptitle("Exp12: 4 K'de SKY130 nFET — yüzer gövde kilitlenmesi ve akım filamanı (disipatif desen)", fontsize=13)
fig.tight_layout()
fig.savefig(OUT / "fig_exp12.png", dpi=105)
json.dump({str(k): v["hist"] for k, v in runs.items()}, open(OUT / "spatial_fig_runs.json", "w"))
print("ok")
