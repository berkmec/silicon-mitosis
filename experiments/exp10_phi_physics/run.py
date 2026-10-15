"""Experiment 10: physical meaning of phi. phi ~ (n - n0)/n0 in the SKY130 inversion layer.
 (1) local inverse compressibility K(n) = d^2(E/A)/dn^2 with kinetic + exchange (+ correlation) -> sign = phi double well?
 (2) linear stability of the uniform layer including the gate-screened Hartree term, restricted to q <= k_F
     (LDA-type local energy is not valid at shorter wavelengths): S(q) = K + V_H(q) F(q).  S < 0 somewhere -> division.
Cases (all SKY130 numbers from the official model file: t_ox = 4.148 nm; m* = 0.19, eps_Si = 11.7, eps_ox = 3.9)."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import egas2d as G  # noqa: E402

OUT = ROOT / "results" / "exp10"
OUT.mkdir(parents=True, exist_ok=True)
M, EPS_SI, EPS_OX, TOX = 0.19, 11.7, 3.9, 4.148e-9
EPS_BAR = (EPS_SI + EPS_OX) / 2
Z_AVG = 2.45e-9  # inversion-layer centroid at n_s=1e12 cm^-2 (exp5)
VG = G.V_gate(EPS_SI, EPS_OX, TOX)
FH = G.fang_howard(Z_AVG)

CASES = {
    "A_tek-vadi_g2_ε̄_xc(DMC)": dict(g=2, m_eff=M, V_x=G.V_uniform(EPS_BAR), eps_c=EPS_BAR, with_corr=True),
    "B_Si_g4_HF_ε̄": dict(g=4, m_eff=M, V_x=G.V_uniform(EPS_BAR), eps_c=EPS_BAR, with_corr=False),
    "C_Si_g4_HF_kapı": dict(g=4, m_eff=M, V_x=VG, eps_c=EPS_BAR, with_corr=False),
    "D_Si_g4_HF_kapı+korelasyon(vekil)": dict(g=4, m_eff=M, V_x=VG, eps_c=EPS_BAR, with_corr=True),
    "E_Si_g4_HF_kapı+korelasyon+FangHoward": dict(g=4, m_eff=M, V_x=VG, eps_c=EPS_BAR, with_corr=True, form=FH),
}

ns_cm2 = np.logspace(9, 13, 49)
ns = ns_cm2 * 1e4
res = {"params": dict(m_eff=M, eps_si=EPS_SI, eps_ox=EPS_OX, t_ox_nm=TOX * 1e9, eps_bar=EPS_BAR, z_avg_nm=Z_AVG * 1e9),
       "cases": {}}
for name, kw in CASES.items():
    K = np.array([G.inverse_compressibility(n, **kw) for n in ns])
    Kkin = 2 * np.pi * G.HBAR**2 / (kw["g"] * kw["m_eff"] * G.M0)
    # sign change
    idx = np.where(np.sign(K[:-1]) != np.sign(K[1:]))[0]
    n_c = None
    if len(idx):
        i = idx[-1]
        n_c = float(np.exp(np.interp(0.0, [K[i], K[i + 1]] if K[i] < K[i + 1] else [K[i + 1], K[i]],
                                     np.log([ns_cm2[i], ns_cm2[i + 1]]) if K[i] < K[i + 1] else np.log([ns_cm2[i + 1], ns_cm2[i]]))))
    # stability with Hartree for q <= kF (form factor on Hartree only in case E)
    form = kw.get("form")
    Smin, Rmin = [], []
    for n, k in zip(ns, K):
        kF = np.sqrt(4 * np.pi * n / kw["g"])
        qs = np.linspace(kF / 200, kF, 200)
        VH = np.array([VG(q) * (1.0 if form is None else form(q)) for q in qs])
        S = k + VH
        Smin.append(float(S.min()))
        Rmin.append(float(VH.min() / abs(k)) if k < 0 else None)
    unstable = [float(ns_cm2[i]) for i, s in enumerate(Smin) if s < 0]
    res["cases"][name] = dict(n_cm2=ns_cm2.tolist(), K_over_Kkin=(K / Kkin).tolist(),
                             rs=[float(G.rs_of(n, M, EPS_BAR)) for n in ns],
                             n_c_cm2=n_c, rs_c=None if n_c is None else float(G.rs_of(n_c * 1e4, M, EPS_BAR)),
                             S_min=Smin, Hartree_over_absK_min=Rmin, unstable_n_cm2=unstable)
    print(f"{name}: K<0 below n_c = {('—' if n_c is None else f'{n_c:.3g}')} cm^-2 (rs_c = "
          f"{('—' if n_c is None else f'{G.rs_of(n_c*1e4, M, EPS_BAR):.2f}')});  uniform layer unstable (q<=kF) at "
          f"{len(unstable)} of {len(ns)} densities{'' if not unstable else f' (max n {max(unstable):.3g})'};  "
          f"min Hartree/|K| = {min([r for r in Rmin if r is not None], default=float('nan')):.2f}", flush=True)

# SKY130 operating point
n1 = 1e12 * 1e4
op = {}
for name, kw in CASES.items():
    K = G.inverse_compressibility(n1, **kw)
    kF = np.sqrt(4 * np.pi * n1 / kw["g"])
    op[name] = dict(K_over_Kkin=float(K / (2 * np.pi * G.HBAR**2 / (kw["g"] * M * G.M0))),
                    VH_at_kF_over_absK=float(VG(kF) / abs(K)))
res["operating_point_1e12"] = dict(rs=float(G.rs_of(n1, M, EPS_BAR)), cases=op)
print("n_s = 1e12 cm^-2: rs =", round(G.rs_of(n1, M, EPS_BAR), 2), {k: {kk: round(vv, 2) for kk, vv in v.items()} for k, v in op.items()})
(OUT / "phi_physics.json").write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")

fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
for name, c in res["cases"].items():
    ax[0].semilogx(c["n_cm2"], c["K_over_Kkin"], label=name)
    ax[1].semilogx(c["n_cm2"], np.array(c["S_min"]) / np.abs(np.array(c["K_over_Kkin"]) * 0 + 1) * 0 + np.array(c["S_min"]) / (2 * np.pi * G.HBAR**2 / (4 * M * G.M0)), label=name)
for a in ax:
    a.axhline(0, c="k", lw=0.8)
    a.axvspan(1e11, 1e13, color="0.92", zorder=0)
ax[0].set(xlabel="n_s (cm⁻²)", ylabel="K / K_kinetik", title="Yerel ters kompresibilite (<0: φ çift kuyusu)")
ax[1].set(xlabel="n_s (cm⁻²)", ylabel="min_q≤kF [K + V_H(q)] / K_kin(g=4)", title="Kapı-ekranlı Hartree dahil kararlılık (<0: bölünme)")
ax[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT / "fig_phi_physics.png", dpi=140)
print("ok")
