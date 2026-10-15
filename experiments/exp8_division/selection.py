"""Domain-size selection: noise -> pattern in the gate-screened Coulomb-frustrated Cahn-Hilliard model.
Compare the measured structure-factor wavenumber with the linear-theory fastest-growing wavenumber and report A_c(d)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.phasefield import OKField  # noqa: E402

OUT = ROOT / "results" / "exp8"
OUT.mkdir(parents=True, exist_ok=True)
N, L = 128, 100.0
DS = (np.inf, 12.0, 6.0, 3.0, 1.5)
rows = []


def q_star(f, phi0):
    qs = np.linspace(1e-3, 3.0, 30000)
    s = f.growth_rate(qs, phi0)
    return float(qs[np.argmax(s)]), float(s.max())


# ---- A_c(d) (uniform state at phi0 = 0 loses stability) ----------------------------------------
crit = {}
for d in DS:
    f = OKField(nx=64, ny=64, Lx=L, Ly=L, A=0.1, d=d)
    crit["inf" if np.isinf(d) else d] = f.critical_A(0.0)
print("A_c(d), phi0=0:", {k: round(v, 4) for k, v in crit.items()}, flush=True)
res = {"A_c_phi0_0": crit, "runs": []}

grid = [(0.0, A) for A in (0.1, 0.2, 0.3)] + [(-0.4, A) for A in (0.05, 0.15)]
for phi0, A in grid:
    for d in DS:
        f = OKField(nx=N, ny=N, Lx=L, Ly=L, kappa=1.0, A=A, d=d)
        qs, smax = q_star(f, phi0)
        rng = np.random.default_rng(7)
        phi = phi0 + 0.05 * rng.standard_normal((N, N))
        phi, _ = f.run(phi, t_end=3000.0, dt=0.25)
        qpk = f.structure_peak(phi)
        amp = float((phi.max() - phi.min()) / 2)
        patterned = bool(phi.std() > 0.2)
        row = dict(phi0=phi0, A=A, d=("inf" if np.isinf(d) else d), q_lin_star=qs, sigma_max=smax, q_peak=qpk,
                   wavelength_peak=float(2 * np.pi / qpk), wavelength_lin=float(2 * np.pi / qs) if qs > 0 else None,
                   phi_std=float(phi.std()), patterned=patterned, amplitude=amp, mass_drift=float(abs(phi.mean() - phi0)))
        res["runs"].append(row)
        print(f"phi0={phi0:>5} A={A:<5} d={str(row['d']):>4}: sigma_max={smax:.4f} q_lin*={qs:.3f} (lam {2*np.pi/qs:5.1f}) | "
              f"q_peak={qpk:.3f} (lam {2*np.pi/qpk:5.1f}) std={phi.std():.2f} patterned={patterned}", flush=True)
        np.save(OUT / f"sel_phi0{phi0}_A{A}_d{row['d']}.npy", phi.astype(np.float32))
        (OUT / "selection.json").write_text(json.dumps(res, indent=2))
print("done")
