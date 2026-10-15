"""Stage 4 preparation: statistical power for a device-array measurement, tunnel resistance / current scale of an accidental dot,
and Coulomb-oscillation visibility at the measurement temperatures (orthodox theory on the predicted E_C)."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import model as m  # noqa: E402
from qcell import schro  # noqa: E402

OUT = ROOT / "results" / "exp7"
OUT.mkdir(parents=True, exist_ok=True)
res = {}

# ---- (1) statistical power ------------------------------------------------------
# H10 predicts a fraction p of min-size devices showing CB-like oscillations: 0.04 (strict) .. 0.12 (optimistic).
# Null: p = 0 (no accidental dots) -> any positive device contradicts it only if false positives are controlled.
rows = []
for p in (0.02, 0.04, 0.08, 0.12):
    for n in (24, 48, 96, 192, 384):
        p_ge1 = 1 - (1 - p) ** n
        p_ge3 = 1 - stats.binom.cdf(2, n, p)
        exp_pos = n * p
        # Clopper-Pearson 95% CI width at the expected count
        k = int(round(exp_pos))
        lo, hi = (stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0, stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0)
        rows.append(dict(p_true=p, n_devices=n, P_at_least_1=float(p_ge1), P_at_least_3=float(p_ge3),
                         expected_positives=float(exp_pos), ci95_at_expected=[float(lo), float(hi)]))
res["power"] = rows
# smallest n for P(>=3 positives) >= 0.95 at each p
need = {}
for p in (0.02, 0.04, 0.08, 0.12):
    n = 3
    while 1 - stats.binom.cdf(2, n, p) < 0.95:
        n += 1
    need[str(p)] = n
res["n_for_95pct_to_see_ge3_positives"] = need
print("n for >=3 positives with 95% power:", need)

# ---- (2) tunnel resistance of an accidental dot ------------------------------------------
# barriers formed by dopant potential: height Vb, width w; effective mass 0.19; WKB transmission T = exp(-2 kappa w)
# R_T ~ (h/e^2)/T (single channel).  Current at V_sd = 100 uV through two equal barriers (series -> R_T total ~ 2 R_T).
H_E2 = m.HBAR * 2 * np.pi / m.E**2  # h/e^2 (ohm)
tun = []
for Vb_meV in (10, 20, 30, 50):
    for w_nm in (5, 10, 15, 20):
        kap = np.sqrt(0.19 * Vb_meV * 1e-3 / schro.HB)
        T = float(np.exp(-2 * kap * w_nm))
        R = H_E2 / T
        I_100uV = 100e-6 / (2 * R)
        tun.append(dict(barrier_meV=Vb_meV, width_nm=w_nm, T=T, R_T_ohm=float(R), I_at_100uV_A=float(I_100uV),
                        CB_ok=bool(R > 3 * H_E2)))
res["tunnel"] = tun
res["h_over_e2_ohm"] = float(H_E2)

# ---- (3) CB visibility vs temperature for the predicted dot ------------------------------------
# single dot, C_sigma from E_C/kB = 30 K, symmetric barriers, orthodox theory; peak-to-valley of G(Vg) at small bias
vis = []
for EC_K in (25.0, 34.0, 40.0):
    Cs = m.E**2 / (2 * EC_K * m.KB)
    Cg_, Cs_, Cd_ = 0.5 * Cs, 0.25 * Cs, 0.25 * Cs
    Vg = np.linspace(0, 3 * m.E / Cg_, 400)
    for T in (0.1, 0.3, 1.0, 4.0, 10.0, 20.0):
        Vsd = 20e-6
        I = np.abs(m.single_dot_current(Cs_, Cd_, Cg_, Vg, [Vsd], T, R_s=1e9, R_d=1e9)[:, 0])
        # peak-to-valley of Coulomb oscillations after removing mean
        pv = float((I.max() - I.min()) / (I.max() + I.min())) if I.max() + I.min() > 0 else 0.0
        vis.append(dict(EC_over_kB_K=EC_K, T_K=T, EC_over_kT=EC_K / T, visibility=pv, I_peak_A=float(I.max())))
res["visibility"] = vis
for v in vis:
    if v["EC_over_kB_K"] == 34.0:
        print(f"E_C=34 K, T={v['T_K']:>5} K: E_C/kT={v['EC_over_kT']:6.1f}  peak-valley visibility={v['visibility']:.3f}  I_peak={v['I_peak_A']:.2e} A")
(OUT / "testplan.json").write_text(json.dumps(res, indent=2))
print("ok")
