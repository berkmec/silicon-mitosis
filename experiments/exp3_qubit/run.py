"""Steps B+C: double quantum dot (two cells) in Si/SiO2, FEM Schroedinger + single-phonon golden rule
+ qubit figures of merit at 4/77/300 K + Qiskit Aer gate fidelity.

IMPORTANT: the phonon part gives only the single-phonon acoustic limit (an UPPER bound on T2);
RT dephasing from optical/two-phonon/charge-noise channels is NOT computed and is handled as explicit scenarios."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import model as m  # noqa: E402
from qcell import phonon, schro  # noqa: E402

OUT = ROOT / "results" / "exp3"
OUT.mkdir(parents=True, exist_ok=True)
M_IN, M_OUT, V0 = 0.26, 0.5, 3.1
MARGIN = 0.6
TEMPS = (4.0, 77.0, 300.0)
TAU_PHI = (10e-15, 100e-15, 1e-12)  # RT dephasing scenarios (s): NOT predictions
rows = []
configs = [(L, d) for L in (1.0, 1.5, 2.0) for d in (0.5, 0.8, 1.2)]
for L, d in configs:
    h_int = 0.12
    gm = schro.build_dqd(L, d, MARGIN, h_int, 0.3)
    nt = gm.mesh.t.shape[1]
    inwell = np.isin(gm.region, gm.well_vols)
    ev, nfree, basis, psis = schro.solve_eigs(gm, np.where(inwell, 0.0, V0), np.where(inwell, 1 / M_IN, 1 / M_OUT),
                                              k=4, return_vecs=True)
    delta_eV = ev[1] - ev[0]
    delta_J = delta_eV * m.E
    ov = phonon.Overlap(basis, psis[0], psis[1])
    q0 = delta_J / (m.HBAR * phonon.V_LA)
    f2 = ov.mean_F2(q0)
    x12 = ov.dipole(0)
    f2_dip = (q0 * x12 * 1e-9) ** 2 / 3.0
    row = dict(L_nm=L, d_nm=d, dofs=nfree, E1_eV=float(ev[0]), E2_eV=float(ev[1]), delta_eV=float(delta_eV),
               gap_to_next_eV=float(ev[2] - ev[1]), x12_nm=x12, q0_per_nm=q0 * 1e-9, meanF2=f2, meanF2_dipole=f2_dip, temps={})
    for T in TEMPS:
        c = phonon.coherence(delta_J, f2, T)
        c["P_thermal"] = c.pop("P")
        c["Q_acoustic"] = (delta_J / m.HBAR) * c["T2_max"]
        sc = {}
        for tp in TAU_PHI:
            T2 = 1.0 / (1.0 / c["T2_max"] + 1.0 / tp)
            sc[f"{tp*1e15:.0f}fs"] = dict(T2_s=T2, Q=float((delta_J / m.HBAR) * T2))
        c["scenarios"] = sc
        row["temps"][str(int(T))] = {k: (float(v) if not isinstance(v, dict) else v) for k, v in c.items()}
    rows.append(row)
    t300 = row["temps"]["300"]
    print(f"L={L} d={d}: Delta={delta_eV*1e3:7.1f} meV  x12={x12:.3f} nm  <F2>={f2:.2e} (dip {f2_dip:.2e})  "
          f"P(300K)={t300['P_thermal']:.3f}  T1_ac(300K)={1/t300['Gamma1']:.2e} s", flush=True)
(OUT / "dqd.json").write_text(json.dumps(rows, indent=2))
