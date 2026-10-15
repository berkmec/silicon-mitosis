"""FEM helpers: gmsh geometry -> scikit-fem mesh, Maxwell capacitance matrix, node sets.

Units of the capacitance solver: epsilon0 * [length of the geometry].
"""
import gmsh
import numpy as np
from scipy.sparse.linalg import splu
from skfem import Basis, BilinearForm, ElementTetP0, ElementTetP1, MeshTet
from skfem.helpers import dot, grad


def gmsh_begin(name="m"):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.clear()
    gmsh.model.add(name)


def gmsh_end():
    if gmsh.isInitialized():
        gmsh.finalize()


def box_bounds(cx, cy, cz, ex, ey=None, ez=None):
    ey = ex if ey is None else ey
    ez = ex if ez is None else ez
    return (cx - ex / 2, cx + ex / 2, cy - ey / 2, cy + ey / 2, cz - ez / 2, cz + ez / 2)


def add_box(bnd):
    x0, x1, y0, y1, z0, z1 = bnd
    return gmsh.model.occ.addBox(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0)


def surfaces_in_box(bnd, tol=1e-6):
    """Surface tags whose bounding box lies inside bnd (+tol)."""
    x0, x1, y0, y1, z0, z1 = bnd
    out = []
    for _, s in gmsh.model.getEntities(2):
        bx = gmsh.model.getBoundingBox(2, s)
        if (bx[0] >= x0 - tol and bx[3] <= x1 + tol and bx[1] >= y0 - tol and bx[4] <= y1 + tol
                and bx[2] >= z0 - tol and bx[5] <= z1 + tol):
            out.append(s)
    return out


def set_mesh_size(h_max, h_fine=None, near_surfaces=None, dist_min=None, dist_max=None):
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    if h_fine is not None and near_surfaces:
        f1 = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f1, "SurfacesList", list(near_surfaces))
        gmsh.model.mesh.field.setNumber(f1, "Sampling", 40)
        f2 = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f2, "InField", f1)
        gmsh.model.mesh.field.setNumber(f2, "SizeMin", h_fine)
        gmsh.model.mesh.field.setNumber(f2, "SizeMax", h_max)
        gmsh.model.mesh.field.setNumber(f2, "DistMin", dist_min if dist_min is not None else h_fine)
        gmsh.model.mesh.field.setNumber(f2, "DistMax", dist_max if dist_max is not None else 4 * h_max)
        gmsh.model.mesh.field.setAsBackgroundMesh(f2)


class GmshMesh:
    """Current (meshed) gmsh model converted to a compact scikit-fem MeshTet."""

    def __init__(self):
        ntags, coords, _ = gmsh.model.mesh.getNodes()
        coords = coords.reshape(-1, 3)
        order = np.argsort(ntags)
        self._ntags = ntags[order]
        coords = coords[order]
        tets, regions = [], []
        for _, vt in gmsh.model.getEntities(3):
            etypes, _, enodes = gmsh.model.mesh.getElements(3, vt)
            for et, nn in zip(etypes, enodes):
                if et == 4:
                    tets.append(np.searchsorted(self._ntags, np.asarray(nn)).reshape(-1, 4))
                    regions.append(np.full(len(tets[-1]), vt))
        t_all = np.vstack(tets)
        self.region = np.concatenate(regions)
        used = np.unique(t_all)
        self._remap = -np.ones(len(self._ntags), dtype=int)
        self._remap[used] = np.arange(len(used))
        self.p = coords[used]
        self.mesh = MeshTet(self.p.T.copy(), self._remap[t_all].T.copy())

    def surface_nodes(self, surface_tags):
        out = []
        for s in surface_tags:
            nt, _, _ = gmsh.model.mesh.getNodes(2, s, includeBoundary=True)
            idx = self._remap[np.searchsorted(self._ntags, nt)]
            out.append(idx[idx >= 0])
        return np.unique(np.concatenate(out))


def capacitance_matrix(mesh, eps_elem, node_sets):
    """Maxwell capacitance matrix C_ij = phi_i^T K phi_j (units eps0*length).

    node_sets: dict name -> node index array (Dirichlet conductors). Every conductor is
    explicit, so rows sum to zero up to solver error (checked by the caller).
    """
    basis = Basis(mesh, ElementTetP1())
    p0 = basis.with_element(ElementTetP0())
    eps = p0.interpolate(np.asarray(eps_elem, dtype=float))

    @BilinearForm
    def a(u, v, w):
        return w["eps"] * dot(grad(u), grad(v))

    K = a.assemble(basis, eps=eps).tocsr()
    n = K.shape[0]
    names = list(node_sets)
    fixed = np.concatenate([node_sets[k] for k in names])
    free = np.setdiff1d(np.arange(n), fixed)
    lu = splu(K[free][:, free].tocsc())
    Kff_fix = K[free][:, fixed]
    phis = []
    for nm in names:
        phi = np.zeros(n)
        phi[node_sets[nm]] = 1.0
        phi[free] = lu.solve(-(Kff_fix @ phi[fixed]))
        phis.append(phi)
    P = np.array(phis)
    return names, P @ (K @ P.T), len(free)
