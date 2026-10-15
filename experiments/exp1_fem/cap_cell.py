"""3D FEM capacitance of the nested cell.

Study A: island + source + drain inside a gate cavity of edge b = a/r. Sweep r, mesh level, scale.
Study B: island inside a floating shell with a window (width w), inside a grounded box.
Results in units eps0*a (study A) / eps0 (study B, fixed geometry).
"""
import csv
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


def run_cell(r, a=1.0, t_frac=0.1, s_frac=0.3, level=1):
    t, s, b = t_frac * a, s_frac * a, a / r
    xs = a / 2 + t + s / 2
    I = fem.box_bounds(0, 0, 0, a)
    S = fem.box_bounds(-xs, 0, 0, s)
    D = fem.box_bounds(+xs, 0, 0, s)
    cav = fem.box_bounds(0, 0, 0, b)
    fem.gmsh_begin("cell")
    vc = fem.add_box(cav)
    tools = [(3, fem.add_box(I)), (3, fem.add_box(S)), (3, fem.add_box(D))]
    gmsh.model.occ.cut([(3, vc)], tools)
    gmsh.model.occ.synchronize()
    sI, sS, sD = (fem.surfaces_in_box(x) for x in (I, S, D))
    sG = [x for _, x in gmsh.model.getEntities(2) if x not in sI + sS + sD]
    div, frac = {1: (2.0, 0.15), 2: (3.0, 0.10), 3: (4.5, 0.07)}[level]
    h_fine = t / div
    h_max = b * frac
    fem.set_mesh_size(h_max, h_fine, sI + sS + sD, dist_min=t, dist_max=0.5 * b)
    gmsh.model.mesh.generate(3)
    gm = fem.GmshMesh()
    ns = {"I": gm.surface_nodes(sI), "S": gm.surface_nodes(sS), "D": gm.surface_nodes(sD), "G": gm.surface_nodes(sG)}
    names, C, nfree = fem.capacitance_matrix(gm.mesh, np.ones(gm.mesh.t.shape[1]), ns)
    nn, nt = int(gm.mesh.p.shape[1]), int(gm.mesh.t.shape[1])
    fem.gmsh_end()
    C = C / a  # -> units eps0*a
    return dict(r=r, a=a, level=level, nodes=nn, tets=nt,
                c_sigma=float(C[0, 0]), c_s=float(-C[0, 1]), c_d=float(-C[0, 2]), c_g=float(-C[0, 3]),
                row_sum_I=float(C[0].sum()), sym_err=float(np.abs(C - C.T).max()))


SHELL_LEVELS = {1: (0.08, 0.8), 2: (0.05, 0.6), 3: (0.035, 0.45)}


def run_shell(w, level=2, B=2.0, wall=0.2, a=0.5, box=6.0, with_shell=True):
    h_fine, h_max = SHELL_LEVELS[level]
    b = B - 2 * wall
    Gb = fem.box_bounds(0, 0, 0, box)
    Ib = fem.box_bounds(0, 0, 0, a)
    fem.gmsh_begin("shell")
    g = fem.add_box(Gb)
    isl = fem.add_box(Ib)
    tools = [(3, isl)]
    if with_shell:
        outer = fem.add_box(fem.box_bounds(0, 0, 0, B))
        inner = fem.add_box(fem.box_bounds(0, 0, 0, b))
        shell, _ = gmsh.model.occ.cut([(3, outer)], [(3, inner)])
        if w > 0:
            win = gmsh.model.occ.addBox(-w / 2, -w / 2, b / 2 - 1e-9 - wall * 0.0, w, w, wall + 0.5)
            shell, _ = gmsh.model.occ.cut(shell, [(3, win)])
        tools += shell
    gmsh.model.occ.cut([(3, g)], tools)
    gmsh.model.occ.synchronize()
    sI = fem.surfaces_in_box(Ib)
    sP = [x for x in fem.surfaces_in_box(fem.box_bounds(0, 0, 0, B)) if x not in sI] if with_shell else []
    sG = [x for _, x in gmsh.model.getEntities(2) if x not in sI + sP]
    fem.set_mesh_size(h_max, h_fine, sI + sP, dist_min=0.1, dist_max=1.2)
    gmsh.model.mesh.generate(3)
    gm = fem.GmshMesh()
    ns = {"I": gm.surface_nodes(sI), "G": gm.surface_nodes(sG)}
    if with_shell:
        ns["P"] = gm.surface_nodes(sP)
    names, C, _ = fem.capacitance_matrix(gm.mesh, np.ones(gm.mesh.t.shape[1]), ns)
    nn = int(gm.mesh.p.shape[1])
    fem.gmsh_end()
    ix = {n: i for i, n in enumerate(names)}
    if not with_shell:
        c_eff = C[ix["I"], ix["I"]]
        return dict(w=None, level=level, nodes=nn, C_II=float(c_eff), C_charging=float(c_eff))
    i, p_ = ix["I"], ix["P"]
    # P floating (Q_P = 0): eliminate P -> charging capacitance of the island
    c_ch = C[i, i] - C[i, p_] * C[p_, i] / C[p_, p_]
    return dict(w=w, level=level, nodes=nn, C_II=float(C[i, i]), C_IP=float(-C[i, p_]),
                C_charging=float(c_ch), row_sum_I=float(C[i].sum()))


if __name__ == "__main__":
    A = []
    for r in (0.5, 0.4, 0.3, 0.2):
        for level in (1, 2, 3):
            row = run_cell(r, level=level)
            A.append(row)
            print("A", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    with open(OUT / "cap_cell_A.csv", "w", newline="") as f:
        wri = csv.DictWriter(f, fieldnames=list(A[0]))
        wri.writeheader()
        wri.writerows(A)

    Bres = []
    for level in (1, 2, 3):
        Bres.append(run_shell(0.0, level=level, with_shell=False))
        print("B baseline (kabuk yok)", Bres[-1], flush=True)
    for w in (0.0, 0.2, 0.4, 0.8, 1.2):
        for level in ((1, 2, 3) if w in (0.0, 1.2) else (2,)):
            row = run_shell(w, level=level)
            Bres.append(row)
            print("B", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    (OUT / "cap_shell_B.json").write_text(json.dumps(Bres, indent=2))
