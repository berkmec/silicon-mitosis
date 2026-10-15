import numpy as np

from qcell import chain, model, phonon, schro


def test_chain_single_island_matches_single_dot_theory():
    Cs, Cd, Cg = 1e-18, 1e-18, 1.5e-18
    C = Cs + Cd + Cg
    ch = chain.Chain([[C]], [Cs], [Cd], [Cg], R=1e6, nmax=6)
    for Vg, Vsd, T in [(0.05, 0.03, 4.0), (0.12, 0.05, 20.0), (0.2, 0.08, 4.0)]:
        _, i_chain = ch.stationary(Vg, Vsd, T)
        i_ref = model.single_dot_current(Cs, Cd, Cg, [Vg], [Vsd], T)[0, 0]
        assert np.isclose(i_chain, i_ref, rtol=1e-6, atol=1e-18), (Vg, Vsd, T, i_chain, i_ref)


def test_chain_zero_bias_zero_current():
    C = np.array([[3e-18, -0.5e-18], [-0.5e-18, 2e-18]])
    ch = chain.Chain(C, [1e-18, 0], [0, 1e-18], [0.8e-18, 0.6e-18], nmax=3)
    _, i = ch.stationary(0.1, 0.0, 10.0)
    assert abs(i) < 1e-20


def test_overlap_normalisation_and_dipole_limit():
    L, m = 3.0, 0.26
    gm = schro.build_well(L, 0.0, 0.4, 0.4, barrier=False)
    n = gm.mesh.t.shape[1]
    ev, _, basis, psis = schro.solve_eigs(gm, np.zeros(n), np.ones(n) / m, k=2, return_vecs=True)
    ov11 = phonon.Overlap(basis, psis[0], psis[0])
    assert np.isclose(ov11.F(np.zeros((3,))).real, 1.0, atol=1e-6)
    ov12 = phonon.Overlap(basis, psis[0], psis[1])
    assert abs(ov12.F(np.zeros((3,)))) < 1e-8  # orthogonality
    q = 0.02 / (L * 1e-9)  # q*L = 0.02
    for ax in range(3):
        d = ov12.dipole(ax) * 1e-9
        qv = np.zeros(3)
        qv[ax] = q
        assert np.isclose(abs(ov12.F(qv)) ** 2, (q * d) ** 2, rtol=1e-2, atol=1e-30)


def test_refined_quadrature_exact_for_polynomials():
    for lv in (0, 1):
        X, W = phonon.refined_tet_quadrature(6, lv)
        assert np.isclose(W.sum(), 1 / 6)
        assert np.isclose((W * X[0] ** 2).sum(), 1 / 60)
        assert np.isclose((W * X[0] * X[1] * X[2]).sum(), 1 / 720)
        assert np.isclose((W * X[0] ** 4).sum(), 24 / 5040)
