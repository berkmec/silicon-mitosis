"""3D effective-mass Schroedinger eigenproblem (BenDaniel-Duke) with scikit-fem P2 elements.

Units: nm, eV. HB = hbar^2 / (2 m0) in eV nm^2.
-HB * div( (1/m) grad psi ) + V psi = E psi,  psi = 0 on the outer box boundary.
"""
import gmsh
import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import eigsh
from skfem import Basis, BilinearForm, ElementTetP0, ElementTetP2
from skfem.helpers import dot, grad

from . import fem

HB = 0.0380998212  # eV nm^2


def build_well(L, margin, h_int, h_max, barrier=True):
    """Return GmshMesh. barrier=True: well cube + barrier shell; else only the well cube."""
    fem.gmsh_begin("well")
    if barrier:
        outer = gmsh.model.occ.addBox(-(L / 2 + margin), -(L / 2 + margin), -(L / 2 + margin),
                                      L + 2 * margin, L + 2 * margin, L + 2 * margin)
        well = gmsh.model.occ.addBox(-L / 2, -L / 2, -L / 2, L, L, L)
        gmsh.model.occ.fragment([(3, outer)], [(3, well)])
    else:
        gmsh.model.occ.addBox(-L / 2, -L / 2, -L / 2, L, L, L)
    gmsh.model.occ.synchronize()
    iface = fem.surfaces_in_box((-L / 2, L / 2, -L / 2, L / 2, -L / 2, L / 2))
    fem.set_mesh_size(h_max, h_int, iface, dist_min=0.3, dist_max=max(L / 3, 1.0))
    gmsh.model.mesh.generate(3)
    gm = fem.GmshMesh()
    # which gmsh volume is the well?
    well_vol = None
    for _, v in gmsh.model.getEntities(3):
        bx = gmsh.model.getBoundingBox(3, v)
        if abs((bx[3] - bx[0]) - L) < 1e-6:
            well_vol = v
    fem.gmsh_end()
    gm.well_vol = well_vol
    return gm


def build_dqd(L, d, margin, h_int, h_max):
    """Two cubic wells (edge L) separated by a barrier of thickness d along x, inside a barrier box.
    Returns GmshMesh with .well_vols (the two well volume tags)."""
    cx = (L + d) / 2
    fem.gmsh_begin("dqd")
    outer = gmsh.model.occ.addBox(-(cx + L / 2 + margin), -(L / 2 + margin), -(L / 2 + margin),
                                  2 * (cx + L / 2 + margin), L + 2 * margin, L + 2 * margin)
    w1 = gmsh.model.occ.addBox(-cx - L / 2, -L / 2, -L / 2, L, L, L)
    w2 = gmsh.model.occ.addBox(cx - L / 2, -L / 2, -L / 2, L, L, L)
    gmsh.model.occ.fragment([(3, outer)], [(3, w1), (3, w2)])
    gmsh.model.occ.synchronize()
    iface = fem.surfaces_in_box((-cx - L / 2, -cx + L / 2, -L / 2, L / 2, -L / 2, L / 2)) + \
        fem.surfaces_in_box((cx - L / 2, cx + L / 2, -L / 2, L / 2, -L / 2, L / 2))
    fem.set_mesh_size(h_max, h_int, iface, dist_min=0.3, dist_max=1.2)
    gmsh.model.mesh.generate(3)
    gm = fem.GmshMesh()
    wells = []
    for _, v in gmsh.model.getEntities(3):
        bx = gmsh.model.getBoundingBox(3, v)
        if abs((bx[3] - bx[0]) - L) < 1e-6 and abs((bx[4] - bx[1]) - L) < 1e-6:
            wells.append(v)
    fem.gmsh_end()
    gm.well_vols = wells
    return gm


def solve_eigs(gm, V_elem, minv_elem, k=5, potential_fn=None, return_vecs=False):
    mesh = gm.mesh
    basis = Basis(mesh, ElementTetP2())
    p0 = basis.with_element(ElementTetP0())
    invm = p0.interpolate(np.asarray(minv_elem, float))
    vel = p0.interpolate(np.asarray(V_elem, float))

    @BilinearForm
    def kform(u, v, w):
        pot = w["V"] if potential_fn is None else potential_fn(w.x)
        return HB * w["invm"] * dot(grad(u), grad(v)) + pot * u * v

    @BilinearForm
    def mform(u, v, w):
        return u * v

    K = kform.assemble(basis, V=vel, invm=invm).tocsr()
    M = mform.assemble(basis).tocsr()
    bd = basis.get_dofs().all()
    free = np.setdiff1d(np.arange(K.shape[0]), bd)
    Kf, Mf = K[free][:, free].tocsc(), M[free][:, free].tocsc()
    vals, vecs = eigsh(Kf, k=k, M=Mf, sigma=0.0, which="LM")
    order = np.argsort(vals)
    if not return_vecs:
        return vals[order], len(free)
    psis = []
    for j in order:
        psi = np.zeros(K.shape[0])
        psi[free] = vecs[:, j]
        psi /= np.sqrt(psi @ (M @ psi))  # M-normalised
        psis.append(psi)
    return vals[order], len(free), basis, psis


def e1d_finite(L, V0, m):
    """Ground state of a 1D finite square well, equal masses (eV)."""
    def f(E):
        kk = np.sqrt(E * m / HB)
        kap = np.sqrt((V0 - E) * m / HB)
        return kk * np.tan(kk * L / 2) - kap
    emax = min(V0, HB * (np.pi / L) ** 2 / m) * (1 - 1e-9)
    return brentq(f, 1e-9, emax)


def hardwall(L, m, n=(1, 1, 1)):
    return HB * np.pi ** 2 / (m * L ** 2) * sum(i * i for i in n)
