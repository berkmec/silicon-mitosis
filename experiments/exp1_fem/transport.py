"""Combine FEM capacitances (cap_cell_A.csv) and FEM level spacings (schrodinger_wells.csv):
charging-energy thresholds per generation and Coulomb diamonds (orthodox theory on FEM C).
"""
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import model as m  # noqa: E402

OUT = ROOT / "results" / "exp1"
KB_EV = m.KB / m.E
rowsA = list(csv.DictReader(open(OUT / "cap_cell_A.csv")))
wells = list(csv.DictReader(open(OUT / "schrodinger_wells.csv")))
Lw = np.array([float(w["L_nm"]) for w in wells])
dEw = np.array([float(w["dE_eV"]) for w in wells])
ratio_edge = float(wells[-1]["ratio_dE_fem_over_hw"])


def cap_row(r, level=3):
    return next(x for x in rowsA if float(x["r"]) == r and int(x["level"]) == level)


def cap_unc(r):
    c3, c2 = float(cap_row(r, 3)["c_sigma"]), float(cap_row(r, 2)["c_sigma"])
    return abs(c3 - c2) / c3


def dE_eV(L_nm):
    if L_nm < Lw.min():
        return None  # outside FEM range
    if L_nm <= Lw.max():
        return float(np.exp(np.interp(np.log(L_nm), np.log(Lw), np.log(dEw))))
    return 3 * m.HBAR**2 * np.pi**2 / (2 * 0.26 * m.M0 * (L_nm * 1e-9) ** 2) / m.E * ratio_edge


res = {"eps": {}, "diamond_gen": {}, "cap_rel_uncertainty": {}}
for r in (0.5, 0.4, 0.3, 0.2):
    res["cap_rel_uncertainty"][str(r)] = cap_unc(r)
MARGIN = 10.0
for eps_r in (3.9, 11.7):
    r = 0.5
    c = cap_row(r)
    table = []
    for g in range(0, 9):
        a = 100e-9 * r**g
        C = eps_r * m.EPS0 * a * float(c["c_sigma"])
        Ec = m.E**2 / (2 * C)
        dE = dE_eV(a * 1e9)
        table.append(dict(gen=g, L_nm=a * 1e9, C_aF=C * 1e18, EC_over_kB_K=Ec / m.KB,
                          dE_over_kB_K=None if dE is None else dE / KB_EV))
    firsts = {}
    for T in (300, 77, 4):
        gc = next((t["gen"] for t in table if t["EC_over_kB_K"] >= MARGIN * T), None)
        gq = next((t["gen"] for t in table if t["dE_over_kB_K"] is not None and t["dE_over_kB_K"] >= MARGIN * T), None)
        firsts[str(T)] = dict(coulomb_gen=gc, coulomb_L_nm=None if gc is None else table[gc]["L_nm"],
                              confinement_gen=gq, confinement_L_nm=None if gq is None else table[gq]["L_nm"])
    res["eps"][str(eps_r)] = dict(table=table, first_regime=firsts)

# diamonds: eps=3.9 (SiO2 fill), T=4 K and 77 K, first Coulomb generation
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
eps_r = 3.9
r = 0.5
c = cap_row(r)
for ax, T in zip(axs, (4.0, 77.0)):
    g = res["eps"][str(eps_r)]["first_regime"][str(int(T))]["coulomb_gen"]
    if g is None:
        ax.set_title(f"{T:.0f} K: Coulomb eşiğine ulaşılmadı")
        continue
    a = 100e-9 * r**g
    k = eps_r * m.EPS0 * a
    Cs, Cd, Cg = float(c["c_s"]) * k, float(c["c_d"]) * k, float(c["c_g"]) * k
    Ec = m.E**2 / (2 * (Cs + Cd + Cg))
    Vg = np.linspace(0, 3 * m.E / Cg, 240)
    Vsd = np.linspace(-4 * Ec / m.E, 4 * Ec / m.E, 240)
    I = m.single_dot_current(Cs, Cd, Cg, Vg, Vsd, T)
    dI = np.gradient(I, Vsd, axis=1)
    v = np.percentile(np.abs(dI), 99)
    im = ax.pcolormesh(Vg * 1e3, Vsd * 1e3, dI.T, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto")
    invalid = " [<1 nm: model geçersiz]" if a < 1e-9 else ""
    ax.set(title=f"{T:.0f} K, nesil {g} (L={a*1e9:.1f} nm){invalid}\nFEM C: E_C/kT={Ec/(m.KB*T):.1f}",
           xlabel="V_gate (mV)", ylabel="V_sd (mV)")
    fig.colorbar(im, ax=ax, label="dI/dV (S)")
    res["diamond_gen"][str(int(T))] = dict(gen=g, L_nm=a * 1e9, Cs_aF=Cs * 1e18, Cd_aF=Cd * 1e18, Cg_aF=Cg * 1e18,
                                           EC_meV=Ec / m.E * 1e3)
fig.suptitle("Coulomb elmasları: FEM kapasitans matrisi üstünde ortodoks teori (sabit etkileşim)")
fig.tight_layout()
fig.savefig(OUT / "fig_diamonds_fem.png", dpi=140)
plt.close(fig)
(OUT / "transport.json").write_text(json.dumps(res, indent=2))
for e, d in res["eps"].items():
    print("eps_r", e, d["first_regime"])
print("diamond gens", res["diamond_gen"])
print("cap rel. belirsizlik (|c3-c2|/c3):", res["cap_rel_uncertainty"])
