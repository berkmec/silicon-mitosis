"""Experiment 0: recursive cell scaling, quantum-regime thresholds, cascade limits,
Coulomb diamonds. Outputs to E:/quantum_cell/results/exp0 (+ desktop copy by caller)."""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from qcell import model as m  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "exp0"
OUT.mkdir(parents=True, exist_ok=True)
TEMPS = {"300 K": 300.0, "77 K": 77.0, "4 K": 4.0}
KBT = {k: m.KB * v for k, v in TEMPS.items()}

# ---- 1. generation table for the baseline cell -----------------------------
base = m.CellParams()
G = np.arange(0, 12)
rows = []
for g in G:
    rows.append(dict(
        gen=int(g), L_nm=m.size(base, g) * 1e9,
        C_aF=m.c_total(base, g) * 1e18,
        EC_meV=m.charging_energy(base, g) / m.E * 1e3,
        EC_over_kB_K=m.charging_energy(base, g) / m.KB,
        dE_meV=m.level_spacing(base, g) / m.E * 1e3,
        dE_over_kB_K=m.level_spacing(base, g) / m.KB,
        V_spawn_V=m.spawn_voltage(base, g), V_break_V=m.breakdown_voltage(base, g)))
with open(OUT / "generations.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

# ---- 2. parameter sweep: first generation reaching each regime -------------
sweep = []
for r in (0.7, 0.5, 0.35):
    for boost in (1, 10, 100):
        p = m.CellParams(r=r, boost=boost)
        for tname, T in TEMPS.items():
            gc, gq, gb = m.regime_generation(p, T)
            # deepest cascade reachable if pulse amplitude is ramped to the largest allowed
            sweep.append(dict(r=r, boost=boost, T=tname, g_coulomb=gc, g_confine=gq, g_both=gb,
                              L_both_nm=None if gb is None else m.size(p, gb) * 1e9,
                              spawn_ok_V=None if gb is None else m.spawn_voltage(p, gb - 1 if gb > 0 else 0),
                              break_V=None if gb is None else m.breakdown_voltage(p, gb - 1 if gb > 0 else 0)))
with open(OUT / "sweep.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(sweep[0]))
    w.writeheader()
    w.writerows(sweep)

# ---- fig 1: scaling ---------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
L = np.array([r["L_nm"] for r in rows])
ax[0].loglog(L, [r["EC_over_kB_K"] for r in rows], "o-", label="E_C / k_B  (Coulomb)")
ax[0].loglog(L, [r["dE_over_kB_K"] for r in rows], "s-", label="ΔE / k_B  (confinement)")
for tname, T in TEMPS.items():
    ax[0].axhline(base.margin * T, ls="--", c="gray", lw=0.8)
    ax[0].text(L.max(), base.margin * T * 1.1, f"10·T = {tname}", ha="right", fontsize=8)
ax[0].invert_xaxis()
ax[0].set(xlabel="cell edge L (nm)  →  smaller", ylabel="temperature scale (K)",
          title="Quantum scales vs generation (r=0.5, boost=1)")
ax[0].legend()
gens = np.array([r["gen"] for r in rows])
ax[1].semilogy(gens, [r["V_spawn_V"] for r in rows], "o-", label="min pulse to spawn child")
ax[1].semilogy(gens, [r["V_break_V"] for r in rows], "x--", label="dielectric breakdown limit")
ax[1].set(xlabel="generation", ylabel="pulse amplitude (V)", title="Cascade window")
ax[1].legend()
fig.tight_layout()
fig.savefig(OUT / "fig1_scaling.png", dpi=140)
plt.close(fig)

# ---- fig 2: Coulomb diamonds at the first Coulomb generation per T ---------
fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))
for j, (tname, T) in enumerate(TEMPS.items()):
    gc = m.regime_generation(base, T)[0]
    if gc is None:
        for i in range(2):
            axs[i, j].set_title(f"{tname}: not reached")
        continue
    for i, g in enumerate((max(gc - 2, 0), gc)):
        Cg = m.c_gate(base, g)
        Cs = Cd = 0.5 * 2 * np.pi * m.EPS0 * base.eps_r * m.size(base, g)
        Ec = m.charging_energy(base, g)
        Vg = np.linspace(0, 3 * m.E / Cg, 220)
        Vsd = np.linspace(-4 * Ec / m.E, 4 * Ec / m.E, 220)
        I = m.single_dot_current(Cs, Cd, Cg, Vg, Vsd, T)
        # dI/dV in log-friendly form
        dI = np.gradient(I, Vsd, axis=1)
        v = np.percentile(np.abs(dI), 99)
        im = axs[i, j].pcolormesh(Vg * 1e3, Vsd * 1e3, dI.T, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto")
        axs[i, j].set(title=f"{tname}, gen {g} (L={m.size(base, g)*1e9:.1f} nm, E_C/kT={Ec/(m.KB*T):.1f})",
                      xlabel="V_gate (mV)", ylabel="V_sd (mV)")
        fig.colorbar(im, ax=axs[i, j], label="dI/dV (S)")
fig.suptitle("Coulomb diamonds: two generations before the threshold (top) vs at it (bottom)")
fig.tight_layout()
fig.savefig(OUT / "fig2_diamonds.png", dpi=140)
plt.close(fig)

# ---- fig 3: cascade depth vs pulse amplitude -------------------------------
vp = np.logspace(-3, 1, 400)
fig, ax = plt.subplots(figsize=(7, 4.5))
for boost in (1, 10, 100):
    p = m.CellParams(boost=boost)
    depth = [m.cascade_depth(p, v)[0] for v in vp]
    ax.semilogx(vp, depth, label=f"boost={boost}")
ax.set(xlabel="fixed pulse amplitude (V)", ylabel="deepest generation reached",
       title="Pulse amplitude sets cascade depth (r=0.5)")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "fig3_cascade_depth.png", dpi=140)
plt.close(fig)

# ---- summary ---------------------------------------------------------------
def fmt(x, d=1):
    return "—" if x is None else f"{x:.{d}f}"

lines = ["# Deney 0 — Özyinelemeli hücre: özet", "",
         "Model: küresel öz-kapasitans + kapı kapasitansı, kübik kutu hapsi, sabit etkileşim modeli.",
         "Eşik ölçütü: enerji ≥ 10·k_B·T.", "",
         "## Temel hücre (r=0.5, boost=1, L0=100 nm)", "",
         "| nesil | L (nm) | E_C/k_B (K) | ΔE/k_B (K) | spawn V (V) | kırılma V (V) |", "|---|---|---|---|---|---|"]
for r in rows:
    lines.append(f"| {r['gen']} | {r['L_nm']:.2f} | {r['EC_over_kB_K']:.0f} | {r['dE_over_kB_K']:.0f} | "
                 f"{r['V_spawn_V']:.3f} | {r['V_break_V']:.3f} |")
lines += ["", "## Eşiğe ilk ulaşılan nesil (Coulomb | hapsolma | ikisi birden) ve o nesildeki kenar (nm)", "",
          "| r | boost | T | Coulomb | hapsolma | ikisi | L (nm) |", "|---|---|---|---|---|---|---|"]
for s in sweep:
    lines.append(f"| {s['r']} | {s['boost']} | {s['T']} | {s['g_coulomb']} | {s['g_confine']} | {s['g_both']} | {fmt(s['L_both_nm'])} |")
(OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
print("done", OUT)
