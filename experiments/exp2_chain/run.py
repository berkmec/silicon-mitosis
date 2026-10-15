"""Step A: graded dot chain S - I0(a) - I1(r a) - I2(r^2 a) - D with a common gate box.
FEM Maxwell capacitance matrix (scikit-fem) + orthodox master equation (qcell.chain)."""
import json
import sys
from pathlib import Path

import gmsh
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import chain, fem  # noqa: E402
from qcell import model as m  # noqa: E402

OUT = ROOT / "results" / "exp2"
OUT.mkdir(parents=True, exist_ok=True)
R_RATIO = 0.5
A0_NM = 12.0
EPS_R = 3.9
LEVELS = {1: (0.02, 0.4), 2: (0.0125, 0.3)}


def chain_caps(level):
    s = [R_RATIO**i for i in range(3)]
    gaps = [0.05, 0.025]
    x0 = 0.0
    x1 = x0 + s[0] / 2 + gaps[0] + s[1] / 2
    x2 = x1 + s[1] / 2 + gaps[1] + s[2] / 2
    sD, sS = 0.125, 0.5
    xD = x2 + s[2] / 2 + gaps[1] + sD / 2
    xS = x0 - s[0] / 2 - 0.1 - sS / 2
    boxes = {
        "I0": fem.box_bounds(x0, 0, 0, s[0]), "I1": fem.box_bounds(x1, 0, 0, s[1]),
        "I2": fem.box_bounds(x2, 0, 0, s[2]), "S": fem.box_bounds(xS, 0, 0, sS), "D": fem.box_bounds(xD, 0, 0, sD)}
    fem.gmsh_begin("chain")
    G = fem.box_bounds(0.3, 0, 0, 3.4)
    g = fem.add_box(G)
    tools = [(3, fem.add_box(b)) for b in boxes.values()]
    gmsh.model.occ.cut([(3, g)], tools)
    gmsh.model.occ.synchronize()
    surf = {k: fem.surfaces_in_box(b) for k, b in boxes.items()}
    allc = sum(surf.values(), [])
    surf["G"] = [x for _, x in gmsh.model.getEntities(2) if x not in allc]
    h_fine, h_max = LEVELS[level]
    fem.set_mesh_size(h_max, h_fine, allc, dist_min=0.03, dist_max=0.6)
    gmsh.model.mesh.generate(3)
    gm = fem.GmshMesh()
    ns = {k: gm.surface_nodes(v) for k, v in surf.items()}
    names, C, _ = fem.capacitance_matrix(gm.mesh, np.ones(gm.mesh.t.shape[1]), ns)
    nn = int(gm.mesh.p.shape[1])
    fem.gmsh_end()
    ix = {n: i for i, n in enumerate(names)}
    I = [ix["I0"], ix["I1"], ix["I2"]]
    return dict(nodes=nn, C=C, ix=ix, I=I, row_sum_max=float(np.abs(C.sum(axis=1)).max()), sym=float(np.abs(C - C.T).max()))


res = {}
caps = {}
CACHE = OUT / "caps_cache.npz"
if CACHE.exists():
    z = np.load(CACHE, allow_pickle=True)
    caps = z["caps"].item()
else:
    for lv in (1, 2):
        caps[lv] = chain_caps(lv)
    np.savez(CACHE, caps=np.array(caps, dtype=object))
for lv in (1, 2):
    c = caps[lv]
    print(f"level {lv}: nodes={c['nodes']}  row_sum_max={c['row_sum_max']:.1e}  sym={c['sym']:.1e}", flush=True)
C = caps[2]["C"]
ix, Iidx = caps[2]["ix"], caps[2]["I"]
dev = float(np.max(np.abs(caps[2]["C"][np.ix_(Iidx, Iidx)] - caps[1]["C"][np.ix_(Iidx, Iidx)])
                   / np.abs(caps[2]["C"][np.ix_(Iidx, Iidx)])))
res["level_change_max_rel"] = dev
print("island-block max relative change level1->2:", f"{dev:.3f}")

k = EPS_R * m.EPS0 * A0_NM * 1e-9  # farads per unit of the FEM result
C_II = C[np.ix_(Iidx, Iidx)] * k
cS = -C[Iidx, ix["S"]] * k
cD = -C[Iidx, ix["D"]] * k
cG = -C[Iidx, ix["G"]] * k
Cinv = np.linalg.inv(C_II)
Ec = m.E**2 / 2 * np.diag(Cinv)  # charging energy of each island with the others neutral-coupled
res["island_sizes_nm"] = [A0_NM * R_RATIO**i for i in range(3)]
res["C_II_aF"] = (C_II * 1e18).tolist()
res["cS_aF"], res["cD_aF"], res["cG_aF"] = [(x * 1e18).tolist() for x in (cS, cD, cG)]
res["EC_over_kB_K"] = (Ec / m.KB).tolist()
print("island sizes (nm):", res["island_sizes_nm"])
print("E_C/k_B per island (K):", np.round(res["EC_over_kB_K"], 1))

NMAX = 6
ch = chain.Chain(C_II, cS, cD, cG, R=1e6, nmax=NMAX)
Vsd = 0.004
temps = (4.0, 77.0, 300.0)
# gate window: stay where no island is near the state-space truncation (mean n < NMAX-2) at 4 K
coarse = np.arange(0.0, 1.5, 0.02)
Vmax = coarse[-1]
for v in coarse:
    if ch.stationary(v, Vsd, 4.0)[0].max() >= NMAX - 2:
        Vmax = v
        break
Vg = np.linspace(0, Vmax, 181)
res["gate_window_V"] = float(Vmax)
print(f"gate penceresi: 0..{Vmax*1e3:.0f} mV (kesme sınırı NMAX={NMAX})", flush=True)
fig, axs = plt.subplots(2, 3, figsize=(15, 7), sharex=True)
res["oscillations"] = {}
for j, T in enumerate(temps):
    occ = np.zeros((len(Vg), 3))
    cur = np.zeros(len(Vg))
    for i, v in enumerate(Vg):
        occ[i], cur[i] = ch.stationary(v, Vsd, T)
    cur_abs = np.abs(cur)
    assert occ.max() < NMAX - 0.5, f"truncation reached at T={T}: max n = {occ.max():.2f}"
    # detrended oscillation amplitude: remove the slow trend (moving average, 31 pts), compare to the mean current
    win = 31
    ker = np.ones(win) / win
    trend = np.convolve(cur_abs, ker, mode="same")
    dI = (cur_abs - trend)[win // 2:-(win // 2)]
    osc = float((np.percentile(dI, 95) - np.percentile(dI, 5)) / cur_abs.mean()) if cur_abs.mean() > 0 else 0.0
    quant = [float(1 - 4 * np.mean(np.abs(occ[:, i] - np.round(occ[:, i])))) for i in range(3)]
    res["oscillations"][str(int(T))] = dict(
        detrended_osc_amp_over_mean=osc, charge_quantisation_index=quant, I_max_A=float(cur_abs.max()),
        I_min_A=float(cur_abs.min()), n_range_per_island=[[float(occ[:, i].min()), float(occ[:, i].max())] for i in range(3)])
    axs[0, j].plot(Vg * 1e3, cur_abs * 1e12)
    axs[0, j].set(title=f"{T:.0f} K: salınım genliği/ort. = {osc:.3f}", ylabel="|I| (pA)")
    vis = osc
    for i in range(3):
        axs[1, j].plot(Vg * 1e3, occ[:, i], label=f"ada {i} ({A0_NM * R_RATIO**i:.0f} nm)")
    axs[1, j].set(xlabel="V_gate (mV)", ylabel="ortalama elektron sayısı")
    if j == 0:
        axs[1, j].legend()
    print(f"T={T:.0f} K: detrended salınım/ort={osc:.4f}  yük-kuantizasyon indeksi (ada0,1,2)={np.round(quant, 3)}  "
          f"n_ada2 aralığı={occ[:,2].min():.3f}..{occ[:,2].max():.3f}", flush=True)
fig.suptitle("Kademeli küçülen üç adalı zincir (12/6/3 nm, ε=3,9), FEM kapasitans + ana denklem, V_sd = 4 mV")
fig.tight_layout()
fig.savefig(OUT / "fig_chain.png", dpi=140)
plt.close(fig)
(OUT / "chain.json").write_text(json.dumps(res, indent=2))
