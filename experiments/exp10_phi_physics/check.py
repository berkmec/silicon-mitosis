"""Re-check of exp10 ("no division"):
 (1) independent analytic check: V_H(k_F)/K_kin must equal the Thomas-Fermi ratio g/(a_B* k_F) for a uniform dielectric;
 (2) is the result an artefact of restricting q <= k_F?  -> scan q_max = k_F, 2 k_F, 4 k_F;
 (3) where does division appear? -> scan gate distance t_ox and gate dielectric (SiO2 3.9 vs high-k 20)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import egas2d as G  # noqa: E402

OUT = ROOT / "results" / "exp10"
M, EPS_SI = 0.19, 11.7
res = {}

# (1) analytic Thomas-Fermi check (uniform eps_bar, g=4, n=1e12 cm^-2)
eps_bar, g, n = 7.8, 4, 1e16
_, aB = G.units(M, eps_bar)
kF = np.sqrt(4 * np.pi * n / g)
ratio_code = G.V_uniform(eps_bar)(kF) / (2 * np.pi * G.HBAR**2 / (g * M * G.M0))
ratio_tf = g / (aB * kF)
res["thomas_fermi_check"] = dict(code=float(ratio_code), analytic=float(ratio_tf))
print(f"(1) V_H(kF)/K_kin: kod {ratio_code:.4f}, analitik g/(a_B* k_F) {ratio_tf:.4f}")


def unstable_window(t_ox, eps_ox, qmax_factor=1.0, with_corr=True, ns_cm2=np.logspace(10, 13, 31)):
    VG = G.V_gate(EPS_SI, eps_ox, t_ox)
    eps_c = (EPS_SI + eps_ox) / 2
    out = []
    for nc in ns_cm2:
        n = nc * 1e4
        K = G.inverse_compressibility(n, g=4, m_eff=M, V_x=VG, eps_c=eps_c, with_corr=with_corr)
        kF = np.sqrt(np.pi * n)
        qs = np.linspace(kF / 200, qmax_factor * kF, 300)
        S = K + np.array([VG(q) for q in qs])
        out.append((float(nc), float(S.min()), float(qs[np.argmin(S)] / kF)))
    uns = [o for o in out if o[1] < 0]
    return dict(unstable_n_cm2=[o[0] for o in uns], n_max_unstable=max([o[0] for o in uns], default=None),
                q_at_min_over_kF=[o[2] for o in uns])


# (2) q-range sensitivity at the SKY130 stack
res["q_range"] = {}
for qf in (1.0, 2.0, 4.0):
    r = unstable_window(4.148e-9, 3.9, qmax_factor=qf)
    res["q_range"][str(qf)] = r
    print(f"(2) SKY130 (SiO2 4.148 nm), q <= {qf:g} k_F: bölünme n aralığı üst sınırı = {r['n_max_unstable']}")

# (3) gate distance / dielectric map
res["gate_map"] = {}
for eps_ox in (3.9, 20.0):
    for t in (4.148e-9, 2e-9, 1e-9, 0.5e-9):
        r = unstable_window(t, eps_ox)
        key = f"eps_ox={eps_ox},t_ox={t*1e9:g}nm"
        res["gate_map"][key] = r
        print(f"(3) {key}: bölünme var mı (q<=kF, 1e10–1e13 cm^-2)? "
              f"{'EVET, n ≤ ' + format(r['n_max_unstable'], '.2g') + ' cm^-2' if r['n_max_unstable'] else 'hayır'}", flush=True)
(OUT / "check.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
