"""3D FEM check of the image-charge building block used in run.py:
single ionised acceptor (Gaussian-smeared, sigma=1.5 nm) at depth z0 in Si, oxide slab 4.148 nm (eps 3.9) on top with grounded gate,
grounded neutral plane at z = W_d, grounded lateral walls. Compare phi(z=0, r) with the image-method model."""
import json
import sys
from pathlib import Path

import gmsh
import numpy as np
from skfem import Basis, BilinearForm, ElementTetP1, LinearForm
from skfem.helpers import dot, grad
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import fem  # noqa: E402
from qcell import mos1d as M  # noqa: E402

OUT = ROOT / "results" / "exp6"
OUT.mkdir(parents=True, exist_ok=True)
TOX, WD, HALF = 4.148, 81.0, 250.0
EPS_SI, EPS_OX = 11.7, 3.9
SIGMA = 1.5
E_OVER_EPS0 = 18.0953  # V*nm  (e / eps0 in V nm)
KE = E_OVER_EPS0 / (4 * np.pi * EPS_SI)  # V nm : potential of e in Si at 1 nm
D_G = TOX * EPS_SI / EPS_OX


def image_phi(z0, r, nimg=12):
    D = WD + D_G
    s0 = z0 + D_G
    phi = 0.0
    for n in range(-nimg, nimg + 1):
        for sgn, sz in ((+1, s0 + 2 * n * D), (-1, -s0 + 2 * n * D)):
            phi += sgn * (-1.0) * KE / np.sqrt(r**2 + (sz - D_G) ** 2)
    return phi


def fem_phi(z0, h_max, h_fine, rs):
    fem.gmsh_begin("cs")
    si = gmsh.model.occ.addBox(-HALF, -HALF, 0, 2 * HALF, 2 * HALF, WD)
    ox = gmsh.model.occ.addBox(-HALF, -HALF, -TOX, 2 * HALF, 2 * HALF, TOX)
    gmsh.model.occ.fragment([(3, si)], [(3, ox)])
    gmsh.model.occ.synchronize()
    f1 = gmsh.model.mesh.field.add("Ball")
    gmsh.model.mesh.field.setNumber(f1, "XCenter", 0); gmsh.model.mesh.field.setNumber(f1, "YCenter", 0)
    gmsh.model.mesh.field.setNumber(f1, "ZCenter", z0); gmsh.model.mesh.field.setNumber(f1, "Radius", 6.0)
    gmsh.model.mesh.field.setNumber(f1, "Thickness", 25.0)
    gmsh.model.mesh.field.setNumber(f1, "VIn", h_fine); gmsh.model.mesh.field.setNumber(f1, "VOut", h_max)
    f2 = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(f2, "XMin", -70); gmsh.model.mesh.field.setNumber(f2, "XMax", 70)
    gmsh.model.mesh.field.setNumber(f2, "YMin", -70); gmsh.model.mesh.field.setNumber(f2, "YMax", 70)
    gmsh.model.mesh.field.setNumber(f2, "ZMin", -TOX - 0.1); gmsh.model.mesh.field.setNumber(f2, "ZMax", 4.0)
    gmsh.model.mesh.field.setNumber(f2, "VIn", 2.0); gmsh.model.mesh.field.setNumber(f2, "VOut", h_max)
    gmsh.model.mesh.field.setNumber(f2, "Thickness", 30.0)
    f3 = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(f3, "FieldsList", [f1, f2])
    gmsh.model.mesh.field.setAsBackgroundMesh(f3)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(3)
    gm = fem.GmshMesh()
    vols = {v: gmsh.model.getBoundingBox(3, v) for _, v in gmsh.model.getEntities(3)}
    ox_vol = [v for v, b in vols.items() if b[5] <= 1e-6][0]
    fem.gmsh_end()
    eps = np.where(gm.region == ox_vol, EPS_OX, EPS_SI)
    mesh = gm.mesh
    basis = Basis(mesh, ElementTetP1(), intorder=4)
    p0 = basis.with_element(__import__("skfem").ElementTetP0())
    epsf = p0.interpolate(eps)

    @BilinearForm
    def a(u, v, w):
        return w["eps"] * dot(grad(u), grad(v))

    @LinearForm
    def f(v, w):
        x = w.x
        g = np.exp(-((x[0]) ** 2 + (x[1]) ** 2 + (x[2] - z0) ** 2) / (2 * SIGMA**2)) / ((2 * np.pi) ** 1.5 * SIGMA**3)
        return -g * v  # negative unit charge (acceptor)

    K = a.assemble(basis, eps=epsf).tocsr()
    b = f.assemble(basis)
    bd = basis.get_dofs().all()
    free = np.setdiff1d(np.arange(K.shape[0]), bd)
    u = np.zeros(K.shape[0])
    u[free] = splu(K[free][:, free].tocsc()).solve(b[free])
    phi = E_OVER_EPS0 * u  # volts (charge -e is in the source sign)
    pts = np.array([[r, 0.0, 0.0] for r in rs]).T
    vals = basis.probes(pts) @ phi
    return vals, int(mesh.p.shape[1])


rs = np.array([0.0, 5.0, 10.0, 20.0, 40.0, 80.0])
res = {}
for z0 in (10.0, 30.0, 60.0):
    for label, (hm, hf) in {"coarse": (14.0, 2.0), "fine": (9.0, 1.2)}.items():
        v, n = fem_phi(z0, hm, hf, rs)
        ref = np.array([image_phi(z0, r) for r in rs])
        rel = (v - ref) / np.abs(ref)
        res[f"z0={z0}_{label}"] = dict(nodes=n, r_nm=rs.tolist(), phi_fem_mV=(v * 1e3).tolist(), phi_image_mV=(ref * 1e3).tolist(),
                                       rel_diff=rel.tolist())
        print(f"z0={z0:>4} {label:<6} nodes={n:>7}  phi_fem(r=0,10,40) = {v[0]*1e3:7.2f} {v[2]*1e3:7.2f} {v[4]*1e3:7.2f} mV | "
              f"image = {ref[0]*1e3:7.2f} {ref[2]*1e3:7.2f} {ref[4]*1e3:7.2f} mV | rel diff r=0,10,40: "
              f"{rel[0]*100:+.1f}% {rel[2]*100:+.1f}% {rel[4]*100:+.1f}%", flush=True)
(OUT / "fem_check.json").write_text(json.dumps(res, indent=2))
