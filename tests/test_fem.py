import gmsh
import numpy as np

from qcell import fem, schro


def _sphere_C(h):
    fem.gmsh_begin("t")
    big = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
    small = gmsh.model.occ.addSphere(0, 0, 0, 0.5)
    gmsh.model.occ.cut([(3, big)], [(3, small)])
    gmsh.model.occ.synchronize()
    fem.set_mesh_size(h)
    gmsh.model.mesh.generate(3)
    gm = fem.GmshMesh()
    inner = fem.surfaces_in_box((-0.5, 0.5, -0.5, 0.5, -0.5, 0.5))
    outer = [s for _, s in gmsh.model.getEntities(2) if s not in inner]
    ns = {"in": gm.surface_nodes(inner), "out": gm.surface_nodes(outer)}
    _, C, _ = fem.capacitance_matrix(gm.mesh, np.ones(gm.mesh.t.shape[1]), ns)
    fem.gmsh_end()
    return C


def test_sphere_capacitance_matches_analytic_and_rows_sum_to_zero():
    C = _sphere_C(0.12)
    assert abs(C[0, 0] - 4 * np.pi) / (4 * np.pi) < 0.03
    assert np.allclose(C.sum(axis=1), 0, atol=1e-9)
    assert np.allclose(C, C.T, atol=1e-10)


def test_sphere_error_decreases_with_mesh():
    e = [abs(_sphere_C(h)[0, 0] - 4 * np.pi) for h in (0.25, 0.12)]
    assert e[1] < e[0]


def test_hardwall_cube_eigenvalues():
    L, m = 3.0, 0.26
    gm = schro.build_well(L, 0.0, 0.4, 0.4, barrier=False)
    n = gm.mesh.t.shape[1]
    ev, _ = schro.solve_eigs(gm, np.zeros(n), np.ones(n) / m, k=4)
    unit = schro.HB * np.pi**2 / (m * L**2)
    assert np.allclose(ev, unit * np.array([3, 6, 6, 6]), rtol=0.01)
