"""Stage 2a: 1D Poisson-Schroedinger of the SKY130 nfet_01v8 gate stack.
Inputs from the official model (nfet_01v8__tt.pm3.spice): toxe=4.148 nm, ndep=1.7e17 cm^-3, epsrox=3.9, vth0=0.519 V.
Outputs: ns(Vg), extracted VT, comparison with vth0, calibrated effective doping, subbands and layer thickness at threshold."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import mos1d as M  # noqa: E402

OUT = ROOT / "results" / "exp5"
OUT.mkdir(parents=True, exist_ok=True)
VTH0_BSIM = 0.5190093


def sweep(na_cm3, tox=4.148):
    mos = M.Mos(na_cm3=na_cm3, tox_nm=tox, zmax_nm=140.0, dz_nm=0.2)
    phif_ = mos.phif
    pss = np.linspace(1.55 * phif_, 2 * phif_ + 0.42, 20)  # weak inversion -> strong inversion for this doping
    sols = [mos.solve(p) for p in pss]
    vg = np.array([s["Vg"] for s in sols])
    ns = np.array([max(s["ns_cm2"], 1e3) for s in sols])
    return mos, pss, sols, vg, ns


res = {}
# nominal SKY130 doping
mos, pss, sols, vg, ns = sweep(1.7e17)
ext = M.vt_from_curve(vg, ns)
phif = mos.phif
W_dep = np.sqrt(2 * M.EPS_SI * M.EPS0 * 2 * phif / (M.Q * mos.na))
vt_classical = mos.vfb + 2 * phif + M.Q * mos.na * W_dep / mos.cox
res["nominal"] = dict(na_cm3=1.7e17, phiF_V=phif, W_dep_nm=W_dep * 1e9, Cox_uF_cm2=mos.cox * 1e6 * 1e-4,
                      VT_classical_depletion=float(vt_classical), **ext, vth0_bsim4=VTH0_BSIM)
print("nominal:", {k: round(v, 4) if isinstance(v, float) else v for k, v in res["nominal"].items()}, flush=True)

# calibrate effective doping so that VT_linear = vth0 (NOT an independent validation)
lo, hi = 1.7e17, 3e18
for _ in range(18):
    mid = np.sqrt(lo * hi)
    _, _, _, vgm, nsm = sweep(mid)
    vt = M.vt_from_curve(vgm, nsm)["VT_linear_extrap"]
    if vt < VTH0_BSIM:
        lo = mid
    else:
        hi = mid
na_eff = np.sqrt(lo * hi)
mos2, pss2, sols2, vg2, ns2 = sweep(na_eff)
ext2 = M.vt_from_curve(vg2, ns2)
res["calibrated"] = dict(na_eff_cm3=float(na_eff), **ext2, phiF_V=mos2.phif,
                         W_dep_nm=float(np.sqrt(2 * M.EPS_SI * M.EPS0 * 2 * mos2.phif / (M.Q * mos2.na)) * 1e9),
                         spacing_nm=float((na_eff * 1e6) ** (-1 / 3) * 1e9))
print("calibrated N_eff:", res["calibrated"], flush=True)

# subbands and thickness at threshold (constant-charge ns=1e12 cm^-2) for the calibrated stack
def at_ns(mos_, pss_, sols_, target):
    nsv = np.array([max(s["ns_cm2"], 1e3) for s in sols_])
    ps = float(np.interp(np.log(target), np.log(nsv), pss_))
    return mos_.solve(ps)


info = {}
for label, (m_, p_, s_) in {"nominal": (mos, pss, sols), "calibrated": (mos2, pss2, sols2)}.items():
    s = at_ns(m_, p_, s_, 1e12)
    z = m_.z
    nz = s["n"]
    zc = np.trapezoid(z * nz, z) / np.trapezoid(nz, z)
    # lowest subband energies per valley relative to E_F
    sub = {("unprimed" if v["mz"] > 0.5 else "primed"): [float(x) for x in w[:3]] for w, _, v in s["bands"]}
    info[label] = dict(Vg_at_1e12=s["Vg"], psi_s=s["psi_s"], Es_MV_per_cm=s["Es"] * 1e-8, z_centroid_nm=float(zc), subbands_eV=sub)
    print(label, info[label], flush=True)
res["at_ns_1e12"] = info
(OUT / "mos1d.json").write_text(json.dumps(res, indent=2))

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
ax[0].semilogy(vg, ns, "o-", label="PS, N=1,7e17 (SKY130 ndep)")
ax[0].semilogy(vg2, ns2, "s-", label=f"PS, N_eff={na_eff:.2e} (vth0'a kalibre)")
ax[0].axvline(VTH0_BSIM, c="k", ls="--", lw=0.8, label="BSIM4 vth0 = 0,519 V")
ax[0].set(xlabel="V_g (V)", ylabel="n_s (cm⁻²)", title="Evirme yük yoğunluğu"); ax[0].legend(fontsize=8)
s = at_ns(mos2, pss2, sols2, 1e12)
ax[1].plot(mos2.z, s["n"] * 1e-6 / 1e18, label="n(z) (10¹⁸ cm⁻³)")
ax[1].set(xlim=(0, 10), xlabel="z (nm)", ylabel="n (10¹⁸ cm⁻³)", title="Evirme tabakası, n_s=10¹² cm⁻² (kalibre yığın)")
ax[2].plot(mos2.z, (M.EG / 2 + mos2.phif - s["psi"]), label="E_c(z) − E_F (eV)")
ax[2].set(xlim=(0, 20), xlabel="z (nm)", ylabel="eV", title="İletim bandı profili")
for w, _, v in s["bands"]:
    for E in w[:3]:
        ax[2].axhline(E, lw=0.7, c="r" if v["mz"] > 0.5 else "b")
fig.tight_layout()
fig.savefig(OUT / "fig_mos1d.png", dpi=140)
plt.close(fig)
print("ok")
