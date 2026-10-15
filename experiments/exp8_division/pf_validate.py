"""Validation of the Ohta-Kawasaki phase-field solver against analytic results."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.phasefield import OKField  # noqa: E402

OUT = ROOT / "results" / "exp8"
OUT.mkdir(parents=True, exist_ok=True)
res = {}

# ---- V1: kink profile and wall energy (A = 0, plain Cahn-Hilliard) -----------------------------
kappa = 1.0
f = OKField(nx=256, ny=4, Lx=64.0, Ly=1.0, kappa=kappa, A=0.0)
x = np.arange(f.nx) * f.Lx / f.nx
phi0 = np.tanh((x - 16) / 2.0) - np.tanh((x - 48) / 2.0) - 1.0  # two walls, deliberately wrong width
phi = np.tile(phi0[:, None], (1, f.ny))
mass0 = phi.sum()
phi, _ = f.run(phi, t_end=3000.0, dt=0.1)
prof = phi[:, 0]
i0 = np.argmin(np.abs(prof[: f.nx // 2 + 16]))  # zero crossing near x~16
# slope at the centre (central difference) vs analytic 1/sqrt(2 kappa)
dx = f.Lx / f.nx
slope = (prof[i0 + 1] - prof[i0 - 1]) / (2 * dx)
# take the better of the two sign-crossings via interpolation
xz = np.interp(0.0, prof[i0 - 3:i0 + 4], x[i0 - 3:i0 + 4]) if prof[i0 - 3] < prof[i0 + 3] else x[i0]
dphi = np.gradient(prof, dx)
fdens = -prof**2 / 2 + prof**4 / 4 + 0.25
energy_two_walls = np.sum(kappa / 2 * dphi**2 + fdens) * dx
sigma_num = energy_two_walls / 2
sigma_an = 2 * np.sqrt(2) / 3 * np.sqrt(kappa)
res["kink"] = dict(slope_num=float(abs(slope)), slope_analytic=float(1 / np.sqrt(2 * kappa)),
                   sigma_num=float(sigma_num), sigma_analytic=float(sigma_an),
                   rel_err_sigma=float(abs(sigma_num - sigma_an) / sigma_an),
                   mass_drift=float(abs(phi.sum() - mass0)))
print("V1 kink:", {k: round(v, 5) for k, v in res["kink"].items()})

# ---- V2: linear growth rate sigma(q) --------------------------------------------------------
A, d = 0.1, np.inf
f = OKField(nx=128, ny=4, Lx=100.0, Ly=1.0, kappa=1.0, A=A, d=d)
x = np.arange(f.nx) * f.Lx / f.nx
ms = [3, 5, 7, 9, 12]
amp0 = 1e-7
phi = sum(amp0 * np.cos(2 * np.pi * m * x / f.Lx) for m in ms)
phi = np.tile(phi[:, None], (1, f.ny))
t_end, dt = 6.0, 0.002
phi_t, _ = f.run(phi.copy(), t_end=t_end, dt=dt)
rows = []
for m in ms:
    q = 2 * np.pi * m / f.Lx
    a0 = abs(np.fft.fft(phi[:, 0])[m]) * 2 / f.nx
    a1 = abs(np.fft.fft(phi_t[:, 0])[m]) * 2 / f.nx
    sig_num = np.log(a1 / a0) / t_end
    sig_an = float(f.growth_rate(q))
    rows.append(dict(q=float(q), sigma_num=float(sig_num), sigma_analytic=sig_an, rel_err=float(abs(sig_num - sig_an) / abs(sig_an))))
    print(f"V2 q={q:.3f}: sigma num={sig_num:.5f} analytic={sig_an:.5f}  rel err={rows[-1]['rel_err']:.2e}")
res["dispersion"] = rows

# ---- V3: bifurcation: critical A and amplitude exponent -----------------------------------
kappa = 1.0
Ac_an = 2 * kappa * (1 / (3 * kappa)) ** 1.5
f0 = OKField(nx=64, ny=4, Lx=64.0, Ly=1.0, kappa=kappa, A=0.1, d=np.inf)
Ac_num = f0.critical_A(0.0)
q0 = (Ac_an / (2 * kappa)) ** (1 / 3)
res["critical"] = dict(A_c_analytic=float(Ac_an), A_c_numeric=float(Ac_num), q_star_analytic=float(q0))
print(f"V3 A_c analytic={Ac_an:.5f} numeric={Ac_num:.5f}  q*={q0:.4f}")
m_per = 4
Lx = 2 * np.pi * m_per / q0
deltas = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32]
amps = []
rng = np.random.default_rng(3)
for dl in deltas:
    A = Ac_an * (1 - dl)
    f = OKField(nx=128, ny=4, Lx=Lx, Ly=1.0, kappa=kappa, A=A, d=np.inf)
    phi = 1e-3 * rng.standard_normal((f.nx, f.ny))
    phi, _ = f.run(phi, t_end=30000.0, dt=0.25)
    fh = np.abs(np.fft.fft(phi.mean(axis=1))) * 2 / f.nx
    amps.append(float(fh[m_per]))
    print(f"   delta={dl:5.2f}: A={A:.5f} fundamental amplitude = {amps[-1]:.5f}  (other modes max {np.delete(fh[1:20], m_per-1).max():.1e})", flush=True)
sel = np.array(deltas) <= 0.08
slope_fit = np.polyfit(np.log(np.array(deltas)[sel]), np.log(np.array(amps)[sel]), 1)[0]
res["amplitude"] = dict(deltas=deltas, amplitudes=amps, loglog_slope_small_delta=float(slope_fit))
print(f"V3 amplitude ~ delta^{slope_fit:.3f}  (supercritical pitchfork expects 0.5)")
(OUT / "pf_validation.json").write_text(json.dumps(res, indent=2))
