"""Stage 2b: random-dopant potential landscape at the Si/SiO2 interface of a SKY130-like channel (126 nm x 420 nm),
first-order model (point ionised acceptors, grounded gate + neutral-bulk planes via method of images, no screening by carriers).
Outputs: sigma of the surface potential, sigma_VT (vs Stolk analytic formula and the model's Pelgrom coefficient),
and puddle (accidental dot) statistics near threshold."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import mos1d as M  # noqa: E402

OUT = ROOT / "results" / "exp6"
OUT.mkdir(parents=True, exist_ok=True)
Q, EPS0, EPS_SI, EPS_OX = M.Q, M.EPS0, M.EPS_SI, M.EPS_OX
LEFF, WCH = 126.0, 420.0  # nm (drawn 150 nm minus 2*lint, W=0.42 um)
TOX = 4.148
D_G = TOX * EPS_SI / EPS_OX  # gate plane distance in silicon-equivalent nm
KE = Q / (4 * np.pi * EPS_SI * EPS0 * 1e-9)  # V*nm : potential of charge q at 1 nm in Si
H = 2.0  # grid spacing nm


def surface_potential(xs, ys, dop, WD, nimg=4, corr=False):
    """Potential (V) at the interface plane from ionised acceptors (charge -q) at positions dop (x,y,z nm, z>0 into Si).
    corr=True multiplies each dopant's potential by a depth/distance factor fitted to the 3D FEM check (fem_check.py):
    f = 1 + (1.8/z0) * exp(-r^2 / (2 (z0+5)^2))  (FEM/image ratio at r~0: 1.18 at z0=10 nm, 1.07 at 30 nm, 1.035 at 60 nm)."""
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    phi = np.zeros_like(gx)
    D = WD + D_G
    for (x0, y0, z0) in dop:
        s0 = z0 + D_G  # distance of the charge from the gate plane
        r2 = (gx - x0) ** 2 + (gy - y0) ** 2
        loc = np.zeros_like(gx)
        for n in range(-nimg, nimg + 1):
            for sgn, sz in ((+1, s0 + 2 * n * D), (-1, -s0 + 2 * n * D)):
                dz = sz - D_G  # observation point at s = D_G
                loc += sgn * (-1.0) * KE / np.sqrt(r2 + dz**2)  # charge -q
        if corr:
            loc *= 1.0 + (1.8 / max(z0, 1.0)) * np.exp(-r2 / (2 * (z0 + 5.0) ** 2))
        phi += loc
    return phi


def one_sample(rng, na_cm3, WD, nimg=8, corr=False):
    vol = LEFF * WCH * WD * 1e-21  # m^3 (nm^3 = 1e-27 m^3 ... below) handled via cm3
    n_mean = na_cm3 * (LEFF * WCH * WD * 1e-21)  # nm^3 -> cm^3 (1 nm^3 = 1e-21 cm^3)
    n = rng.poisson(n_mean)
    dop = np.column_stack([rng.uniform(0, LEFF, n), rng.uniform(0, WCH, n), rng.uniform(0, WD, n)])
    xs = np.arange(H / 2, LEFF, H)
    ys = np.arange(H / 2, WCH, H)
    return xs, ys, surface_potential(xs, ys, dop, WD, nimg, corr), n


# ---- checks -----------------------------------------------------------------
# (1) single charge, no images (nimg=0 -> one pair of images only) vs analytic image formula at one point
xs = np.array([10.0]); ys = np.array([0.0])
dop = [(0.0, 0.0, 10.0)]
phi_num = surface_potential(xs, ys, dop, WD=1e9, nimg=0)[0, 0]  # bulk plane at infinity: only gate image
r = np.hypot(10.0, 0.0)
dz_dir = (10.0 + D_G) - D_G
dz_img = (-(10.0 + D_G)) - D_G
phi_ref = -KE / np.hypot(r, dz_dir) + KE / np.hypot(r, dz_img)
print(f"[check1] single charge + gate image: numeric {phi_num:.6f} V, analytic {phi_ref:.6f} V")
assert abs(phi_num - phi_ref) < 1e-9
# (2) convergence in number of image terms
rng0 = np.random.default_rng(1)
_, _, p4, _ = one_sample(rng0, 1.7e17, 81.0, nimg=4)
rng0 = np.random.default_rng(1)
_, _, p8, _ = one_sample(rng0, 1.7e17, 81.0, nimg=8)
print(f"[check2] image-sum convergence (nimg 4 vs 8): max |dphi| = {np.max(np.abs(p4 - p8))*1e3:.4f} mV, rms phi = {np.std(p8)*1e3:.2f} mV")

# ---- Monte Carlo ---------------------------------------------------------------
res = {"checks": dict(single_charge_diff_V=float(abs(phi_num - phi_ref)),
                      image_convergence_max_mV=float(np.max(np.abs(p4 - p8)) * 1e3))}
KT = 0.025852
NS = 80
for label, na, WD, corr in (("nominal_1.7e17", 1.7e17, 80.9, False), ("nominal_1.7e17_FEMcorr", 1.7e17, 80.9, True),
                            ("effective_3.9e17", 3.92e17, 54.6, False), ("effective_3.9e17_FEMcorr", 3.92e17, 54.6, True)):
    rng = np.random.default_rng(12345)
    sig_phi, vth_dev, dots = [], [], []
    first = None
    for s in range(NS):
        xs, ys, phi, n = one_sample(rng, na, WD, corr=corr)
        U = -phi  # electron potential energy (eV); repulsive near ionised acceptors
        U = U - U.mean()
        sig_phi.append(float(U.std()))
        vth_dev.append(float(U.mean() + 0.0))  # placeholder, replaced below
        if first is None:
            first = (xs, ys, U.copy())
        # puddles at E_F = U_min + e*dE for several fill levels; isolated = not touching source (x=0) or drain (x=Leff) edges
        U_ = U
        Umin = U_.min()
        for dE in (0.015, 0.030, 0.045):
            mask = U_ < Umin + dE
            lab, nl = ndimage.label(mask)
            touch = set(np.unique(lab[0, :])) | set(np.unique(lab[-1, :]))
            for k in range(1, nl + 1):
                if k in touch:
                    continue
                sel = lab == k
                area = float(sel.sum()) * H * H
                a_nm = float(np.sqrt(area))
                depth = dE - float(U_[sel].min() - Umin)  # eV, fill level above this puddle's own minimum
                e1_inf = 2 * 0.0381 / 0.19 * np.pi**2 / a_nm**2  # eV, square infinite well, m_par = 0.19 (overestimates E1)
                dots.append(dict(sample=s, dE=dE, area_nm2=area, a_nm=a_nm, depth_eV=depth, E1_inf_eV=float(e1_inf),
                                 bound_strict=bool(e1_inf < depth), bound_optimistic=bool(0.5 * e1_inf < depth)))
    sig = float(np.mean(sig_phi))
    # sigma of the channel-averaged potential across samples = sigma_VT (needs un-demeaned means): recompute
    rng = np.random.default_rng(777)
    means = []
    for s in range(NS):
        xs, ys, phi, n = one_sample(rng, na, WD, corr=corr)
        means.append(float(-phi.mean()))
    sigma_vt_mc = float(np.std(means))
    # Stolk et al. 1998 analytic (literature formula, not re-derived here)
    sigma_vt_stolk = float(Q * TOX * 1e-9 / (EPS_OX * EPS0) * np.sqrt(na * 1e6 * WD * 1e-9 / (3 * LEFF * 1e-9 * WCH * 1e-9)))
    summary = dict(na_cm3=na, W_dep_nm=WD, mean_dopants=float(na * LEFF * WCH * WD * 1e-21),
                   sigma_surface_potential_mV=sig * 1e3, sigma_VT_MC_mV=sigma_vt_mc * 1e3,
                   sigma_VT_Stolk_mV=sigma_vt_stolk * 1e3, n_classical_isolated_puddles=len(dots),
                   frac_samples_with_classical_puddle=float(len({d["sample"] for d in dots}) / NS))
    for crit in ("bound_strict", "bound_optimistic"):
        bd = [d for d in dots if d[crit]]
        A = np.array([d["a_nm"] for d in bd])
        entry = dict(n=len(bd), frac_samples=float(len({d["sample"] for d in bd}) / NS))
        if len(A):
            ec = 585.0 / A  # FEM-derived constant scaled to eps_eff=7.8 (rough, geometry mismatch)
            entry.update(a_nm_median=float(np.median(A)), a_nm_p10_p90=[float(np.percentile(A, 10)), float(np.percentile(A, 90))],
                         EC_over_kB_K_median=float(np.median(ec)),
                         EC_over_kB_K_p10_p90=[float(np.percentile(ec, 10)), float(np.percentile(ec, 90))],
                         E1_meV_median=float(np.median([d["E1_inf_eV"] for d in bd]) * 1e3))
        summary[crit] = entry
    res[label] = summary
    print(label, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in summary.items()}, flush=True)
    if label == "nominal_1.7e17":
        fx = first
(OUT / "rdf.json").write_text(json.dumps(res, indent=2))

xs, ys, U = fx
fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
im = ax[0].imshow(U.T * 1e3, origin="lower", extent=[0, LEFF, 0, WCH], aspect="equal", cmap="viridis")
ax[0].set(xlabel="x (nm), kaynak→savak", ylabel="y (nm)", title="Rastgele katkı potansiyeli (elektron U, meV), örnek")
fig.colorbar(im, ax=ax[0])
Umin = U.min()
for dE, c in ((0.015, "w"), (0.030, "r")):
    ax[0].contour(xs, ys, (U < Umin + dE).T.astype(float), levels=[0.5], colors=c, linewidths=0.8)
ax[1].hist([U.ravel() * 1e3], bins=40)
ax[1].set(xlabel="U (meV)", ylabel="piksel sayısı", title=f"σ = {U.std()*1e3:.1f} meV")
fig.tight_layout()
fig.savefig(OUT / "fig_rdf.png", dpi=140)
print("ok")
