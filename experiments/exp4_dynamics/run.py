"""Stage 1 of the roadmap:
 (a) inter-cell tunnel coupling from the FEM double dot, validated against WKB decay,
 (b) dynamic loading = Landau-Zener passage, numerically validated, with thermal relaxation (Lindblad) for the best geometry,
 (c) SKY130 scale gap: charging energy / tunnel coupling vs the 150 nm gate-length floor."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import model as m  # noqa: E402
from qcell import schro  # noqa: E402

OUT = ROOT / "results" / "exp4"
OUT.mkdir(parents=True, exist_ok=True)
rows = json.loads((ROOT / "results" / "exp3" / "dqd.json").read_text())
HBAR_EV = m.HBAR / m.E  # eV s
res = {}

# ---- (a) tunnel splitting vs barrier thickness, WKB check -----------------
M_OUT, V0 = 0.5, 3.1
tbl = []
for L in sorted({r["L_nm"] for r in rows}):
    rr = sorted([r for r in rows if r["L_nm"] == L], key=lambda r: r["d_nm"])
    d = np.array([r["d_nm"] for r in rr])
    delta = np.array([r["delta_eV"] for r in rr])
    E1 = np.mean([r["E1_eV"] for r in rr])
    slope = -np.polyfit(d, np.log(delta), 1)[0]
    kappa = np.sqrt(M_OUT * (V0 - E1) / schro.HB)  # 1/nm
    # local slopes between adjacent d
    loc = [float(-(np.log(delta[i + 1]) - np.log(delta[i])) / (d[i + 1] - d[i])) for i in range(len(d) - 1)]
    tbl.append(dict(L_nm=L, E1_eV=float(E1), fit_decay_per_nm=float(slope), local_decay_per_nm=loc,
                    kappa_WKB_per_nm=float(kappa), ratio_fit_over_WKB=float(slope / kappa)))
    print(f"L={L}: fit decay={slope:.2f}/nm  local={np.round(loc,2)}  WKB kappa={kappa:.2f}/nm  ratio={slope/kappa:.2f}")
res["tunnel_coupling"] = tbl

# ---- (b1) Landau-Zener numerical validation (hbar = 1) --------------------
def lz_prob_diabatic(tc, v, span):
    # H = (v t / 2) sz + tc sx ; start in the lower adiabatic state at t=-span, measure staying diabatic at +span
    def rhs(t, y):
        psi = y[:2] + 1j * y[2:]
        H = np.array([[v * t / 2, tc], [tc, -v * t / 2]])
        d = -1j * (H @ psi)
        return np.concatenate([d.real, d.imag])
    t0 = -span
    H0 = np.array([[v * t0 / 2, tc], [tc, -v * t0 / 2]])
    w, U = np.linalg.eigh(H0)
    psi0 = U[:, 0]  # ground state at t0
    sol = solve_ivp(rhs, [t0, span], np.concatenate([psi0.real, psi0.imag]), rtol=1e-9, atol=1e-11, method="DOP853")
    psi = sol.y[:2, -1] + 1j * sol.y[2:, -1]
    # at t=+span the ground state is mostly |down> (index 1); diabatic = stays in index 0 population
    return float(abs(psi[0]) ** 2)


lz = []
for gamma in (0.1, 0.5, 1.0, 2.0, 4.0):
    tc = 1.0
    v = 2 * np.pi * tc**2 / gamma  # gamma = 2 pi tc^2 / v   (hbar=1)
    p_num = lz_prob_diabatic(tc, v, span=60.0 / np.sqrt(v) * 4 + 40)
    p_th = np.exp(-gamma)
    lz.append(dict(gamma=gamma, P_num=p_num, P_theory=p_th, abs_err=abs(p_num - p_th)))
    print(f"LZ gamma={gamma}: numeric={p_num:.5f} theory={p_th:.5f} err={abs(p_num-p_th):.1e}")
res["lz_validation"] = lz

# ---- (b2) loading map: adiabaticity vs electrical pulse speed -------------
swing_eV = 0.100  # detuning sweep range through the anticrossing
rise_ps = np.array([1, 3, 10, 30, 100, 1000, 10000], float)
load = []
for r in rows:
    t_eV = r["delta_eV"] / 2
    entry = dict(L_nm=r["L_nm"], d_nm=r["d_nm"], t_meV=t_eV * 1e3, gamma_vs_rise_ps={})
    for rp in rise_ps:
        v_eV_per_s = swing_eV / (rp * 1e-12)
        gamma = 2 * np.pi * t_eV**2 / (HBAR_EV * v_eV_per_s)
        entry["gamma_vs_rise_ps"][str(int(rp))] = float(gamma)
    load.append(entry)
res["loading_map"] = load

# ---- (b3) Lindblad sweep at 300 K for the best geometry (L=1,d=0.5) -------
best = max(rows, key=lambda r: r["delta_eV"])
dJ = best["delta_eV"]  # eV
t_e = dJ / 2
c300 = best["temps"]["300"]
Gamma_rel = c300["Gamma1"]  # 1/s total relaxation rate (single-phonon)
n = c300["n_bose"]
g_down, g_up = Gamma_rel * (n + 1) / (2 * n + 1), Gamma_rel * n / (2 * n + 1)
phi = {tp: 1.0 / (c300["scenarios"][tp]["T2_s"]) - 0.5 * Gamma_rel for tp in ("10fs", "100fs", "1000fs")}


def lindblad_transfer(rise_ps, tau_phi_rate):
    """Sweep detuning from -X to +X (eV) in time rise; electron starts in left dot (ground state at -X).
    Returns probability to be in the right dot at the end (successful transfer)."""
    X = swing_eV / 2
    T = rise_ps * 1e-12
    def rhs(t, y):
        rho = (y[:4] + 1j * y[4:]).reshape(2, 2)
        eps = -X + 2 * X * t / T
        H = np.array([[eps / 2, t_e], [t_e, -eps / 2]]) / HBAR_EV  # rad/s (basis: left, right)
        d = -1j * (H @ rho - rho @ H)
        # relaxation/excitation between instantaneous eigenstates; total rate held at its anticrossing value
        # (F(q) at other energies not computed), detailed balance enforced at the INSTANTANEOUS gap
        w, U = np.linalg.eigh(H * HBAR_EV)
        gap = w[1] - w[0]
        bz = np.exp(-gap / (m.KB * 300.0 / m.E))
        g_dn, g_upi = Gamma_rel / (1 + bz), Gamma_rel * bz / (1 + bz)
        L_dn = U[:, 0][:, None] @ U[:, 1].conj()[None, :]
        L_up = L_dn.conj().T
        for L, g in ((L_dn, g_dn), (L_up, g_upi)):
            d += g * (L @ rho @ L.conj().T - 0.5 * (L.conj().T @ L @ rho + rho @ L.conj().T @ L))
        # detuning (charge) noise: pure dephasing in the left/right basis, sigma_z
        sz = np.diag([1.0, -1.0])
        d += 0.5 * tau_phi_rate * (sz @ rho @ sz - rho)
        return np.concatenate([d.real.reshape(-1), d.imag.reshape(-1)])
    rho0 = np.zeros((2, 2), complex)
    rho0[0, 0] = 1.0  # in the left dot at t = 0 (eps = -X: left dot lower in energy)
    sol = solve_ivp(rhs, [0, T], np.concatenate([rho0.real.reshape(-1), rho0.imag.reshape(-1)]), rtol=1e-8, atol=1e-10, method="LSODA")
    rho = (sol.y[:4, -1] + 1j * sol.y[4:, -1]).reshape(2, 2)
    return float(rho[1, 1].real)


lind = []
for rp in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0):
    row = dict(rise_ps=rp, P_transfer={})
    for tp, rate in phi.items():
        row["P_transfer"][tp] = lindblad_transfer(rp, max(rate, 0.0))
    lind.append(row)
    print("Lindblad 300 K L=1,d=0.5:", rp, "ps", {k: round(v, 3) for k, v in row["P_transfer"].items()}, flush=True)
gap_end = np.hypot(swing_eV / 2, 2 * t_e)  # eV, instantaneous gap at the end of the sweep (eps=+X)
bz_end = np.exp(-gap_end / (m.KB * 300.0 / m.E))
p_lo = 1 / (1 + bz_end)
w_right_lo = 0.5 * (1 + (swing_eV / 2) / gap_end)  # right-dot weight of the lower adiabatic state
P_eq_right = p_lo * w_right_lo + (1 - p_lo) * (1 - w_right_lo)
print(f"termal denge P_right(300 K, eps=+X) = {P_eq_right:.4f}")
res["lindblad_300K_best"] = dict(L_nm=best["L_nm"], d_nm=best["d_nm"], delta_meV=dJ * 1e3, T1_s=1 / Gamma_rel,
                                 P_right_thermal_equilibrium=float(P_eq_right), rows=lind)

# ---- (c) SKY130 scale gap -------------------------------------------------
import csv  # noqa: E402

A = list(csv.DictReader(open(ROOT / "results" / "exp1" / "cap_cell_A.csv")))
c_sigma = float(next(x for x in A if float(x["r"]) == 0.5 and int(x["level"]) == 3)["c_sigma"])
K_EC = {}  # E_C/kB = const / a(nm), for eps_r = 3.9 and 11.7
for eps in (3.9, 11.7):
    C_per_nm = eps * m.EPS0 * c_sigma * 1e-9  # F per nm of island edge
    K_EC[str(eps)] = m.E**2 / (2 * C_per_nm) / m.KB  # K*nm
res["EC_const_Knm"] = K_EC
gap = []
for T in (0.05, 0.1, 0.3, 1.0, 4.0, 20.0, 77.0, 300.0):
    gap.append(dict(T_K=T, a_max_nm_eps3p9=K_EC["3.9"] / (10 * T), a_max_nm_eps11p7=K_EC["11.7"] / (10 * T)))
res["required_island_size"] = gap
for a in (150, 210):
    res[f"EC_over_kB_K_at_{a}nm_eps3p9"] = K_EC["3.9"] / a
print("E_C/kB at 150 nm (eps 3.9):", K_EC["3.9"] / 150, "K; 10 kT limit T <=", K_EC["3.9"] / 150 / 10, "K")

# electrostatic tunnel barrier between dots: WKB decay for realistic barrier heights in silicon (m*=0.19)
wk = []
for dV_meV in (5, 20, 50):
    kap = np.sqrt(0.19 * (dV_meV * 1e-3) / schro.HB)  # 1/nm
    wk.append(dict(barrier_meV=dV_meV, kappa_per_nm=float(kap), decay_length_nm=float(1 / kap),
                   t_ratio_at_30nm=float(np.exp(-kap * 30)), t_ratio_at_210nm=float(np.exp(-kap * 210))))
res["silicon_barrier_wkb"] = wk
(OUT / "stage1.json").write_text(json.dumps(res, indent=2))

# ---- figure ---------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
for r_ in tbl:
    pass
Ls = sorted({r["L_nm"] for r in rows})
for L in Ls:
    rr = sorted([r for r in rows if r["L_nm"] == L], key=lambda r: r["d_nm"])
    ax[0].semilogy([r["d_nm"] for r in rr], [r["delta_eV"] * 1e3 for r in rr], "o-", label=f"L={L} nm")
ax[0].set(xlabel="bariyer kalınlığı d (nm)", ylabel="tünel yarılması Δ (meV)", title="FEM: Δ ∝ exp(−κ d)  (WKB ile uyum)")
ax[0].legend()
g = [x["gamma"] for x in lz]
ax[1].semilogx(g, [x["P_theory"] for x in lz], "k-", label="Landau-Zener (analitik)")
ax[1].semilogx(g, [x["P_num"] for x in lz], "ro", label="sayısal (DOP853)")
ax[1].set(xlabel="γ = 2π t²/(ħ v)", ylabel="diyabatik geçiş olasılığı", title="Landau-Zener doğrulaması")
ax[1].legend()
Tg = [x["T_K"] for x in gap]
ax[2].loglog(Tg, [x["a_max_nm_eps3p9"] for x in gap], "o-", label="ε=3,9")
ax[2].loglog(Tg, [x["a_max_nm_eps11p7"] for x in gap], "s-", label="ε=11,7")
ax[2].axhline(150, c="r", ls="--", label="SKY130 poly genişliği 150 nm")
ax[2].set(xlabel="sıcaklık T (K)", ylabel="E_C ≥ 10 kT için en büyük ada kenarı (nm)", title="Süreç ölçeği boşluğu")
ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "fig_stage1.png", dpi=140)
plt.close(fig)
print("ok")
