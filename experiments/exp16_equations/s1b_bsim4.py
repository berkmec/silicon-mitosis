"""Exp16 step 1b (ground truth): does sparse regression recover BSIM4's own laws from the ngspice tables?

Data: experiments/exp12_filament/spice/tab_vb{0.0..1.0}.dat (official SKY130 nfet_01v8 tt card, BSIM4 4.5, W = 1 um,
L = 0.15 um -> card bin model.62; columns Vds,id,Vds,isub,Vds,ibs,Vds,ibd,Vds,vth; rows V_gs-major 0.3..1.8 x V_ds 0..2.6).
(i) impact ionisation: BSIM4 uses I_sub/I_d = a x exp(-beta0/x), x = V_ds - V_dseff, a = alpha0/L_eff + alpha1.
    Per (V_gs, V_bs) curve: library {1, ln x, 1/x, x, x^2, 1/x^2, sqrt x} for y = ln(isub/id), offset V0 scanned;
    V0 and the sparsity threshold are chosen on held-out points (upper 30 % of V_ds: extrapolation test).
(ii) body effect: vth(V_bs, V_ds), V_bs <= 0.7 V (below BSIM4's V_bseff clamp 0.95 Phi_s), library
    {1, s, V_bs, V_bs^2, V_ds, V_ds^2, V_ds V_bs, V_ds s, s V_bs} with s = sqrt(phi - V_bs), phi scanned; trained on
    V_bs <= 0.5, validated on V_bs = 0.6, 0.7 (extrapolation).
Expected values are computed from the card's bin-62 parameters (binunit = 2: lengths in metres).
Output: results/exp16_equations/s1b_bsim4.json, fig_s1b_bsim4.png."""
import common as C  # noqa: I001

import time

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

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
    assert d.shape == (len(VG) * len(VD), 10) and np.allclose(d[:len(VD), 0], VD)
    ID[i], II[i], VTH[i] = (d[:, c].reshape(len(VG), len(VD)) for c in (1, 3, 9))

# ---------------------------------------------------------------- card (bin 62) expectations
LEFF, WEFF = 0.15e-6 - 2 * 1.1932e-8, 1e-6 - 2 * 2.1859e-8          # LINT, WINT
LW = LEFF * WEFF
beta0 = 9.447267864 + 5.553430207e-07 / LEFF + 5.354963848e-06 / WEFF - 6.754537199e-13 / LW
alpha0 = 1.907954345e-07 - 2.028209293e-14 / LEFF - 1.955725927e-13 / WEFF + 2.466874455e-20 / LW
alpha1 = 0.85
a_exp = alpha0 / LEFF + alpha1
k1, lpeb, k3b, toxe = 0.90707349, -7.082e-8, 0.54, 4.148e-9
k2 = -5.039519017e-02 - 8.580291123e-09 / LEFF - 9.232858671e-08 / WEFF + 1.043605365e-14 / LW
tnom = 303.15
vtm0 = 8.617333e-5 * tnom
eg0 = 1.16 - 7.02e-4 * tnom**2 / (tnom + 1108.0)
ni = 1.45e10 * (tnom / 300.15) ** 1.5 * np.exp(21.5565981 - eg0 / (2 * vtm0))
phis = 0.4 + vtm0 * np.log(1.7e17 / ni)
eps_si, q = 11.7 * 8.8541878e-12, 1.602176634e-19
xdep0 = np.sqrt(2 * eps_si * phis / (q * 1.7e23))
lt0 = np.sqrt(11.7 / 3.9 * toxe * xdep0)
th_dibl = 0.5 / (np.cosh(0.45862506 * LEFF / lt0) - 1)
EXPECT = dict(beta0=beta0, alpha0=alpha0, a=a_exp, ln_a=np.log(a_exp), phi_s=phis,
              c_sqrt=k1 * np.sqrt(1 + lpeb / LEFF), c_vbs=-k2 + k3b * toxe * phis / WEFF,
              c_vds=-0.00069413878 * th_dibl, c_vds_vbs=0.043998 * th_dibl, k1_bin=k1, k2_bin_eff=k2,
              theta_dibl=th_dibl, note="user-quoted k1=0.5409, k2=-0.0267 are the long-channel bins model.0/.1 "
                                       "(L >= 8 um); the simulated L = 0.15 um device uses bin model.62")
print("expected", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in EXPECT.items()}, flush=True)

# ---------------------------------------------------------------- (i) impact ionisation
LIB1 = ["1", "ln x", "1/x", "x", "x^2", "1/x^2", "sqrt x"]


def lib1(x):
    return np.column_stack([np.ones_like(x), np.log(x), 1 / x, x, x**2, 1 / x**2, np.sqrt(x)])


TAUS1 = np.logspace(-6, -1, 12)
curves = []
for ib in range(0, len(VB), 2):
    for ig in range(2, len(VG), 2):
        ratio = II[ib, ig] / np.maximum(ID[ib, ig], 1e-30)
        ok = (ratio > 1e-13) & (ID[ib, ig] > 1e-12) & (VD > 0)
        if ok.sum() < 14:
            continue
        vd, y = VD[ok], np.log(ratio[ok])
        best = None
        for V0 in np.arange(0.0, vd[0] - 0.049, 0.02):
            x = vd - V0
            n_tr = int(round(0.7 * len(x)))
            Th = lib1(x)
            path = C.threshold_path(Th[:n_tr], y[:n_tr], Th[n_tr:], y[n_tr:], TAUS1)
            sel = C.select_model(path)
            err = sel["rmse_val"][0] * np.sqrt(np.mean(y[n_tr:] ** 2))
            if best is None or err < best[0]:
                best = (err, V0, sel)
        if best is None:
            continue
        err, V0, sel = best
        coef = sel["coef"][0]
        curves.append(dict(vbs=float(VB[ib]), vgs=float(VG[ig]), n=int(ok.sum()), vd_min=float(vd[0]), V0=float(V0),
                           val_rmse_ln=float(err), support=[LIB1[i] for i in np.flatnonzero(coef)],
                           coef={LIB1[i]: float(coef[i]) for i in np.flatnonzero(coef)}))
print(f"(i) {len(curves)} curves fitted ({time.time() - t_start:.0f}s)", flush=True)
supports = {}
for c in curves:
    key = " + ".join(c["support"])
    supports[key] = supports.get(key, 0) + 1
chyn = [c for c in curves if set(c["support"]) == {"1", "ln x", "1/x"}]


def stats(vals):
    v = np.array(vals)
    return dict(median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)), n=len(v))


part1 = dict(n_curves=len(curves), support_counts=dict(sorted(supports.items(), key=lambda kv: -kv[1])),
             frac_chynoweth_support=len(chyn) / max(len(curves), 1),
             chynoweth_curves=dict(ln_x_coef=stats([c["coef"]["ln x"] for c in chyn]) if chyn else None,
                                   beta0=stats([-c["coef"]["1/x"] for c in chyn]) if chyn else None,
                                   ln_a=stats([c["coef"]["1"] for c in chyn]) if chyn else None,
                                   val_rmse_ln=stats([c["val_rmse_ln"] for c in chyn]) if chyn else None),
             curves=curves)
# fixed-form check (no sparsity): y - ln x = ln a - beta0/x with V0 from the sparse scan, per Chynoweth curve
print("(i) supports:", part1["support_counts"], flush=True)
print("(i) Chynoweth-support stats:", {k: v for k, v in part1["chynoweth_curves"].items()}, flush=True)

# ---------------------------------------------------------------- (ii) body effect
vg_spread = float(np.abs(VTH - VTH[:, :1, :]).max())
LIB2 = ["1", "s", "Vbs", "Vbs^2", "Vds", "Vds^2", "Vds·Vbs", "Vds·s", "s·Vbs"]
VBg, VDg = np.meshgrid(VB, VD, indexing="ij")
vth = VTH[:, 0, :]  # V_gs independent (checked)
mask_tr = (VBg <= 0.5 + 1e-9)
mask_va = (VBg > 0.55) & (VBg <= 0.7 + 1e-9)


def lib2(vbs, vds, phi):
    s = np.sqrt(phi - vbs)
    return np.column_stack([np.ones_like(vbs), s, vbs, vbs**2, vds, vds**2, vds * vbs, vds * s, s * vbs])


TAUS2 = np.logspace(-6, -1, 16)
scan = []
for phi in np.arange(0.71, 1.5001, 0.01):
    path = C.threshold_path(lib2(VBg[mask_tr], VDg[mask_tr], phi), vth[mask_tr],
                            lib2(VBg[mask_va], VDg[mask_va], phi), vth[mask_va], TAUS2)
    sel = C.select_model(path)
    rm = sel["rmse_val"][0] * np.sqrt(np.mean(vth[mask_va] ** 2))
    scan.append((rm, phi, sel))
scan.sort(key=lambda r: r[0])
rm, phi_b, sel = scan[0]
coef2 = sel["coef"][0]
# fixed BSIM4-form least squares (support forced) for comparison, same phi scan
best_fixed = None
for phi in np.arange(0.71, 1.5001, 0.005):
    Th = lib2(VBg[mask_tr], VDg[mask_tr], phi)[:, [0, 1, 2, 4, 6]]
    cf, *_ = np.linalg.lstsq(Th, vth[mask_tr], rcond=None)
    pv = lib2(VBg[mask_va], VDg[mask_va], phi)[:, [0, 1, 2, 4, 6]] @ cf
    e = float(np.sqrt(np.mean((pv - vth[mask_va]) ** 2)))
    if best_fixed is None or e < best_fixed[0]:
        best_fixed = (e, phi, cf)
# bootstrap of the sparse model at the chosen phi/tau
Th_all = lib2(VBg[mask_tr | mask_va], VDg[mask_tr | mask_va], phi_b)
bs = C.bootstrap(Th_all, vth[mask_tr | mask_va], sel["tau"], n_boot=200)
part2 = dict(vth_vgs_independence_max_abs_V=vg_spread,
             sparse=dict(phi=float(phi_b), tau=sel["tau"], val_rmse_V=float(rm),
                         support=[LIB2[i] for i in np.flatnonzero(coef2)],
                         coef={LIB2[i]: float(coef2[i]) for i in np.flatnonzero(coef2)},
                         bootstrap_p_incl={LIB2[i]: float(bs["p_incl"][0, i]) for i in range(len(LIB2))},
                         bootstrap_std={LIB2[i]: float(bs["std"][0, i]) for i in range(len(LIB2))}),
             fixed_bsim4_form=dict(phi=float(best_fixed[1]), val_rmse_V=best_fixed[0],
                                   coef=dict(zip(["1", "s", "Vbs", "Vds", "Vds·Vbs"], map(float, best_fixed[2])))),
             top5_phi=[dict(phi=float(p), val_rmse_V=float(r), support=[LIB2[i] for i in np.flatnonzero(s["coef"][0])])
                       for r, p, s in scan[:5]],
             vth_clamp_region=dict(vbs=VB[7:].tolist(), vth_at_vds0p05=vth[7:, 1].tolist()))
print("(ii) sparse:", part2["sparse"], flush=True)
print("(ii) fixed :", part2["fixed_bsim4_form"], flush=True)

C.save_json("s1b_bsim4.json", dict(expected=EXPECT, impact_ionisation=part1, body_effect=part2,
                                   wall_s=time.time() - t_start))

fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
for c in curves[::7]:
    ib, ig = list(VB).index(c["vbs"]), list(VG).index(c["vgs"])
    ratio = II[ib, ig] / np.maximum(ID[ib, ig], 1e-30)
    ok = (ratio > 1e-13) & (VD > c["V0"] + 1e-6)
    x = VD[ok] - c["V0"]
    cf = np.array([c["coef"].get(n, 0.0) for n in LIB1])
    axs[0].plot(1 / x, np.log(ratio[ok]) - np.log(x), "k.", ms=3)
    axs[0].plot(1 / x, lib1(x) @ cf - np.log(x), "r-", lw=0.8)
axs[0].set_xlabel("1/x = 1/(V_ds − V0)  (1/V)")
axs[0].set_ylabel("ln(I_sub/I_d) − ln x")
axs[0].set_title("Chynoweth/BSIM4 düz çizgi beklentisi: eğim −β0", fontsize=9)
if chyn:
    axs[1].hist([-c["coef"]["1/x"] for c in chyn], bins=20, color="tab:blue")
axs[1].axvline(beta0, color="r", ls="--", label=f"kart β0,eff = {beta0:.2f} V")
axs[1].set_xlabel("bulunan β0 (V), {1, ln x, 1/x} desteğini seçen eğriler")
axs[1].legend(fontsize=8)
for ib in range(len(VB)):
    axs[2].plot(VD, VTH[ib, 0], lw=1, label=f"V_bs={VB[ib]:.1f}" if ib % 2 == 0 else None)
    if VB[ib] <= 0.7 + 1e-9:
        axs[2].plot(VD, lib2(np.full_like(VD, VB[ib]), VD, phi_b) @ coef2, "k:", lw=0.8)
axs[2].set_xlabel("V_ds (V)")
axs[2].set_ylabel("vth (V)")
axs[2].set_title(f"gövde etkisi: renkli BSIM4, nokta seyrek model (φ={phi_b:.2f})", fontsize=9)
axs[2].legend(fontsize=6, ncol=2)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s1b_bsim4.png", dpi=110)
print(f"done {time.time() - t_start:.0f}s")
