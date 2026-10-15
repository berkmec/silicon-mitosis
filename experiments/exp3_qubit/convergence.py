"""Convergence of Delta and <|F(q0)|^2> for the double dot (L=1 nm, d=0.5 nm): mesh size and quadrature order."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import model as m  # noqa: E402
from qcell import phonon, schro  # noqa: E402

OUT = ROOT / "results" / "exp3"
OUT.mkdir(parents=True, exist_ok=True)
M_IN, M_OUT, V0 = 0.26, 0.5, 3.1
L, d = 1.0, 0.5
res = []
for h_int, h_max in ((0.22, 0.6), (0.16, 0.4), (0.12, 0.3)):
    gm = schro.build_dqd(L, d, 0.6, h_int, h_max)
    inwell = np.isin(gm.region, gm.well_vols)
    ev, nfree, basis, psis = schro.solve_eigs(gm, np.where(inwell, 0.0, V0), np.where(inwell, 1 / M_IN, 1 / M_OUT),
                                              k=4, return_vecs=True)
    delta_J = (ev[1] - ev[0]) * m.E
    q0 = delta_J / (m.HBAR * phonon.V_LA)
    for order, levels in (((4, 0), (9, 0), (6, 1)) if h_int == 0.16 else ((6, 1),)):
        f2 = phonon.Overlap(basis, psis[0], psis[1], order=order, levels=levels).mean_F2(q0)
        res.append(dict(h_int=h_int, h_max=h_max, intorder=f"order{order}/sub{levels}", dofs=nfree,
                        delta_meV=float((ev[1] - ev[0]) * 1e3), q0_per_nm=q0 * 1e-9, meanF2=f2,
                        tets=int(gm.mesh.t.shape[1])))
        print(res[-1], flush=True)
(OUT / "convergence.json").write_text(json.dumps(res, indent=2))
