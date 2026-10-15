"""Exp16 step 1b-2: follow-up of s1b (which FAILED: no curve gave the Chynoweth support, body-effect model collinear).
Diagnosis: (i) V0 was scanned on a 0.02 V grid -> x carried up to 10 mV systematic error that dense models absorbed;
(ii) BSIM4 clamps V_bseff smoothly below 0.95 Phi_s, which shifts V_bseff by -3...-8 mV at V_bs = 0.6...0.7 V (the
validation points). Fix, using BSIM4's own definitions of the variables (this is a test of the regression, not of
BSIM4): (i) x = V_ds - V_dseff(V_ds; V_dsat), V_dseff = V_dsat - (V_dsat - V_ds - d + sqrt((V_dsat - V_ds - d)^2
+ 4 d V_dsat))/2, d = DELTA = 0.01, V_dsat fitted CONTINUOUSLY per curve; then the sparse path on
{1, ln x, 1/x, x, x^2, 1/x^2, sqrt x}. (ii) V_bseff with the BSIM4 upper clamp (T0 = 0.95 phi, delta 1e-3), phi fitted
continuously; sparse path on {1, s, Vb, Vb^2, Vds, Vds^2, Vds Vb, Vds s, s Vb} (Vb = V_bseff, s = sqrt(phi - Vb)).
Output: results/exp16_equations/s1b2_bsim4.json"""
import common as C  # noqa: I001

import json
import time

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

C.low_priority()
t_start = time.time()
SP = C.ROOT / "experiments" / "exp12_filament" / "spice"
VB = np.round(np.arange(0.0, 1.0001, 0.1), 2)
VG = np.round(np.arange(0.3, 1.8001, 0.05), 3)
VD = np.round(np.arange(0.0, 2.6001, 0.05), 3)
ID = np.zeros((len(VB), len(VG), len(VD)))
II, VTH = np.zeros_like(ID), np.zeros_like(ID)
for i, vb in enumerate(VB):
    d = np.loadtxt(SP / f"tab_vb{vb:.1f}.dat")
    ID[i], II[i], VTH[i] = (d[:, c].reshape(len(VG), len(VD)) for c in (1, 3, 9))
EXP = json.loads((C.OUT / "s1b_bsim4.json").read_text(encoding="utf-8"))["expected"]
DELTA = 0.01
NAMES1 = ["1", "ln x", "1/x", "x", "x^2", "1/x^2", "sqrt x"]
TAUS = np.logspace(-7, -1, 19)


def xeff(vd, vdsat):
    a = vdsat - vd - DELTA
    return vd - (vdsat - 0.5 * (a + np.sqrt(a * a + 4 * DELTA * vdsat)))


def lib1(x):
    return np.column_stack([np.ones_like(x), np.log(x), 1 / x, x, x**2, 1 / x**2, np.sqrt(x)])


curves = []
for ib in range(0, len(VB), 2):
    for ig in range(2, len(VG), 2):
        ratio = II[ib, ig] / np.maximum(ID[ib, ig], 1e-30)
        ok = (ratio > 1e-13) & (ID[ib, ig] > 1e-12) & (VD > 0)
        if ok.sum() < 14:
            continue
        vd, y = VD[ok], np.log(ratio[ok])
        n_tr = int(round(0.7 * len(vd)))

        def resid(p, sl=slice(None)):
            x = xeff(vd[sl], p[2])
            return p[0] + np.log(np.maximum(x, 1e-9)) - p[1] / np.maximum(x, 1e-9) - y[sl]

        best = None
        for v0 in np.linspace(0.02, max(vd[0] - 0.05, 0.03), 8):
            try:
                f = least_squares(resid, x0=[0.0, 13.0, v0], bounds=([-20, 0.1, 1e-3], [20, 60, vd[0] + 0.3]),
                                  kwargs=dict(sl=slice(0, n_tr)))
            except ValueError:
                continue
            if best is None or f.cost < best.cost:
                best = f
        c0, b, vdsat = best.x
        x = xeff(vd, vdsat)
        val_rmse_fixed = float(np.sqrt(np.mean(resid(best.x, slice(n_tr, None)) ** 2)))
        # free ln-coefficient check (linear LS at the fitted V_dsat)
        A = np.column_stack([np.ones_like(x), np.log(x), -1 / x])
        cfree, *_ = np.linalg.lstsq(A[:n_tr], y[:n_tr], rcond=None)
        path = C.threshold_path(lib1(x[:n_tr]), y[:n_tr], lib1(x[n_tr:]), y[n_tr:], TAUS)
        sel = C.select_model_r2(path, r2_tol=1e-6)
        supp = [NAMES1[i] for i in np.flatnonzero(sel["coef"][0])]
        chyn_taus = [r["tau"] for r in path if set(NAMES1[i] for i in np.flatnonzero(r["coef"][0])) == {"1", "ln x", "1/x"}]
        curves.append(dict(vbs=float(VB[ib]), vgs=float(VG[ig]), n=int(ok.sum()), c0=float(c0), beta0=float(b),
                           vdsat=float(vdsat), val_rmse_ln_fixed_form=val_rmse_fixed,
                           free_ln_coef=float(cfree[1]), free_beta0=float(cfree[2]), free_c0=float(cfree[0]),
                           sparse_support=supp, sparse_r2_val=sel["r2_val"][0],
                           tau_range_chynoweth=[min(chyn_taus), max(chyn_taus)] if chyn_taus else None))
print(f"(i) {len(curves)} curves ({time.time() - t_start:.0f}s)", flush=True)


def st(v):
    v = np.asarray(v, float)
    return dict(median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)))


sup_counts = {}
for c in curves:
    k = " + ".join(c["sparse_support"])
    sup_counts[k] = sup_counts.get(k, 0) + 1
part1 = dict(n_curves=len(curves), beta0_fixed_form=st([c["beta0"] for c in curves]),
             c0_fixed_form=st([c["c0"] for c in curves]), expected_beta0=EXP["beta0"], expected_ln_a=EXP["ln_a"],
             val_rmse_ln_fixed_form=st([c["val_rmse_ln_fixed_form"] for c in curves]),
             free_ln_coef=st([c["free_ln_coef"] for c in curves]), free_beta0=st([c["free_beta0"] for c in curves]),
             sparse_support_counts=dict(sorted(sup_counts.items(), key=lambda kv: -kv[1])),
             frac_curves_with_chynoweth_on_path=float(np.mean([c["tau_range_chynoweth"] is not None for c in curves])),
             frac_selected_chynoweth=float(np.mean([set(c["sparse_support"]) == {"1", "ln x", "1/x"} for c in curves])),
             curves=curves)
print("(i)", {k: v for k, v in part1.items() if k != "curves"}, flush=True)

# ---------------------------------------------------------------- (ii) body effect with the V_bseff clamp
NAMES2 = ["1", "s", "Vb", "Vb^2", "Vds", "Vds^2", "Vds·Vb", "Vds·s", "s·Vb"]
VBg, VDg = np.meshgrid(VB, VD, indexing="ij")
vth = VTH[:, 0, :]
m_tr, m_va = VBg <= 0.5 + 1e-9, (VBg > 0.55) & (VBg <= 0.7 + 1e-9)


def vbseff(vbs, phi):
    t0 = 0.95 * phi
    t1 = t0 - vbs - 1e-3
    return t0 - 0.5 * (t1 + np.sqrt(t1 * t1 + 4e-3 * t0))


def lib2(vbs, vds, phi):
    vb = vbseff(vbs, phi)
    s = np.sqrt(phi - vb)
    return np.column_stack([np.ones_like(vb), s, vb, vb**2, vds, vds**2, vds * vb, vds * s, s * vb])


TRUE_IDX = [0, 1, 2, 4, 6]


def fixed_err(phi, mask):
    Th = lib2(VBg[mask], VDg[mask], phi)[:, TRUE_IDX]
    cf, *_ = np.linalg.lstsq(Th, vth[mask], rcond=None)
    return float(np.sqrt(np.mean((Th @ cf - vth[mask]) ** 2))), cf


opt = minimize_scalar(lambda p: fixed_err(p, m_tr)[0], bounds=(0.72, 1.3), method="bounded",
                      options=dict(xatol=1e-6))
phi_f = float(opt.x)
e_tr, cf = fixed_err(phi_f, m_tr)
Thv = lib2(VBg[m_va], VDg[m_va], phi_f)[:, TRUE_IDX]
e_va = float(np.sqrt(np.mean((Thv @ cf - vth[m_va]) ** 2)))
path2 = C.threshold_path(lib2(VBg[m_tr], VDg[m_tr], phi_f), vth[m_tr], lib2(VBg[m_va], VDg[m_va], phi_f), vth[m_va],
                         np.logspace(-8, -1, 22))
sel2 = C.select_model_r2(path2, r2_tol=1e-6)
true_taus = [r["tau"] for r in path2 if set(np.flatnonzero(r["coef"][0])) == set(TRUE_IDX)]
part2 = dict(phi_fit=phi_f, phi_expected=EXP["phi_s"], fixed_form_coef=dict(zip([NAMES2[i] for i in TRUE_IDX],
                                                                                 map(float, cf))),
             expected=dict(c_sqrt=EXP["c_sqrt"], c_vbs=EXP["c_vbs"]), rmse_train_V=e_tr, rmse_val_V=e_va,
             sparse_selected=[NAMES2[i] for i in np.flatnonzero(sel2["coef"][0])],
             sparse_coef={NAMES2[i]: float(sel2["coef"][0][i]) for i in np.flatnonzero(sel2["coef"][0])},
             sparse_r2_val=sel2["r2_val"][0],
             tau_range_true_support=[min(true_taus), max(true_taus)] if true_taus else None,
             path=[dict(tau=r["tau"], support=[NAMES2[i] for i in np.flatnonzero(r["coef"][0])],
                        rmse_val=r["rmse_val"][0]) for r in path2])
print("(ii)", {k: v for k, v in part2.items() if k != "path"}, flush=True)
C.save_json("s1b2_bsim4.json", dict(impact_ionisation=part1, body_effect=part2, wall_s=time.time() - t_start))
print(f"done {time.time() - t_start:.0f}s")
