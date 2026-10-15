"""Validation 1: concentric spheres, analytic C = 4*pi*eps*a*b/(b-a). Convergence in h.
Also exports a mesh for the independent dolfinx cross-check."""
import json
import sys
from pathlib import Path

import gmsh
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import fem  # noqa: E402

OUT = ROOT / "results" / "exp1"
OUT.mkdir(parents=True, exist_ok=True)
A, B = 0.5, 1.0
C_EXACT = 4 * np.pi * A * B / (B - A)

rows = []
for h in (0.30, 0.20, 0.14, 0.10, 0.07, 0.05):
    fem.gmsh_begin("sph")
    big = gmsh.model.occ.addSphere(0, 0, 0, B)
    small = gmsh.model.occ.addSphere(0, 0, 0, A)
    gmsh.model.occ.cut([(3, big)], [(3, small)])
    gmsh.model.occ.synchronize()
    fem.set_mesh_size(h)
    gmsh.model.mesh.generate(3)
    gm = fem.GmshMesh()
    inner = fem.surfaces_in_box((-A, A, -A, A, -A, A), tol=1e-6)
    outer = [s for _, s in gmsh.model.getEntities(2) if s not in inner]
    ns = {"inner": gm.surface_nodes(inner), "outer": gm.surface_nodes(outer)}
    names, C, nfree = fem.capacitance_matrix(gm.mesh, np.ones(gm.mesh.t.shape[1]), ns)
    c_in = C[0, 0]
    rows.append(dict(h=h, nodes=int(gm.mesh.p.shape[1]), tets=int(gm.mesh.t.shape[1]),
                     C_fem=float(c_in), rel_err=float(abs(c_in - C_EXACT) / C_EXACT),
                     row_sum=float(C[0].sum())))
    if h == 0.14:
        np.savez(OUT / "sphere_h014_mesh.npz", p=gm.p, t=gm.mesh.t.T.copy(), C_fem=c_in)
    fem.gmsh_end()
    print(f"h={h:<5} nodes={rows[-1]['nodes']:<7} C={c_in:.6f} err={rows[-1]['rel_err']:.3e}", flush=True)

hs = np.array([r["h"] for r in rows])
errs = np.array([r["rel_err"] for r in rows])
slope = np.polyfit(np.log(hs), np.log(errs), 1)[0]
res = dict(C_exact=C_EXACT, rows=rows, loglog_slope=float(slope))
(OUT / "validate_sphere.json").write_text(json.dumps(res, indent=2))
print("C_exact", C_EXACT, "| yakınsama eğimi (beklenen ~2):", round(float(slope), 2))
