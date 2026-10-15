"""Review check (2026-10-11) of exp6: how much of sigma_psi is the deterministic 'dome' E[U](x,y) created by placing
dopants only inside the 126 x 420 nm box (no charge outside, no S/D electrodes)?
Campbell's theorem: E[phi](x,y) = N_A * sum_z dz * (K_z ** box_indicator)(x,y), with the SAME image kernel as run.py
(gate plane + neutral-bulk plane via images, nimg = 8). Random part: Var_rand = <sigma_psi^2> - Var_space(E[U]).
Uses rdf.json for the MC sigma_psi (image-model variants only; the FEM-corrected kernel is not handled here)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import mos1d as M  # noqa: E402

Q, EPS0, EPS_SI, EPS_OX = M.Q, M.EPS0, M.EPS_SI, M.EPS_OX
LEFF, WCH, TOX, H = 126.0, 420.0, 4.148, 2.0
D_G = TOX * EPS_SI / EPS_OX
KE = Q / (4 * np.pi * EPS_SI * EPS0 * 1e-9)


def kernel_z(rx, ry, z0, WD, nimg=8):
    """Potential (V) at the interface from one charge -q at depth z0, same image construction as run.py."""
    D = WD + D_G
    s0 = z0 + D_G
    r2 = rx**2 + ry**2
    out = np.zeros_like(r2)
    for n in range(-nimg, nimg + 1):
        for sgn, sz in ((+1, s0 + 2 * n * D), (-1, -s0 + 2 * n * D)):
            out += sgn * (-1.0) * KE / np.sqrt(r2 + (sz - D_G) ** 2)
    return out


def mean_potential(na_cm3, WD, pad=300.0, dz=2.0):
    nx = int(round((LEFF + 2 * pad) / H))
    ny = int(round((WCH + 2 * pad) / H))
    x = (np.arange(nx) - nx // 2) * H
    y = (np.arange(ny) - ny // 2) * H
    RX, RY = np.meshgrid(x, y, indexing="ij")
    box = np.zeros((nx, ny))
    ix0, iy0 = int(round(pad / H)), int(round(pad / H))
    box[ix0:ix0 + int(LEFF / H), iy0:iy0 + int(WCH / H)] = 1.0
    rho_area = na_cm3 * 1e-21 * H * H * dz  # dopants per (H x H x dz) cell
    phi = np.zeros((nx, ny))
    Fb = np.fft.rfft2(box)
    for z0 in np.arange(dz / 2, WD, dz):
        K = kernel_z(np.fft.ifftshift(RX), np.fft.ifftshift(RY), z0, WD)
        phi += np.fft.irfft2(Fb * np.fft.rfft2(K), s=box.shape) * rho_area
    U = -phi[ix0:ix0 + int(LEFF / H), iy0:iy0 + int(WCH / H)]  # electron potential energy (eV) inside the box
    return U


rdf = json.loads((ROOT / "results" / "exp6" / "rdf.json").read_text(encoding="utf-8"))
out = {}
for label, na, WD in (("nominal_1.7e17", 1.7e17, 80.9), ("effective_3.9e17", 3.92e17, 54.6)):
    U = mean_potential(na, WD)
    var_fixed = float(U.var())
    sig_mc = None
    for k, v in rdf.items():
        if isinstance(v, dict) and k == label:
            sig_mc = v.get("sigma_surface_potential_mV")
    cx, cy = U.shape[0] // 2, U.shape[1] // 2
    res = dict(sigma_fixed_mV=np.sqrt(var_fixed) * 1e3, center_minus_y_edge_meV=float(U[cx, cy] - U[cx, 0]) * 1e3,
               center_minus_x_edge_meV=float(U[cx, cy] - U[0, cy]) * 1e3,
               center_minus_corner_meV=float(U[cx, cy] - U[0, 0]) * 1e3, sigma_psi_MC_mV=sig_mc)
    if sig_mc:
        res["fixed_fraction_of_variance"] = var_fixed / (sig_mc * 1e-3) ** 2
        res["sigma_random_mV"] = float(np.sqrt(max((sig_mc * 1e-3) ** 2 - var_fixed, 0.0)) * 1e3)
    out[label] = res
    print(label, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res.items()}, flush=True)
(ROOT / "results" / "exp6" / "review_campbell.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
