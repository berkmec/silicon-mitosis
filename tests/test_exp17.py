"""Exp17: filament split law (experiments/exp17_silicon_split/law.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

EXP = Path(__file__).resolve().parents[1] / "experiments" / "exp17_silicon_split"
pytestmark = pytest.mark.skipif(not (Path(__file__).resolve().parents[1] / "experiments" / "exp12_filament" / "spice"
                                     / "tab_vb0.0.dat").exists(), reason="ngspice tables not generated")
sys.path.insert(0, str(EXP))


@pytest.fixture(scope="module")
def c():
    import law
    return law.constants()


def test_constants_match_exp12_and_exp16(c):
    assert c["V_M"] == pytest.approx(0.8018, abs=5e-4)          # exp12 analytic Maxwell voltage
    assert c["j_on"] * 1e6 == pytest.approx(2.55, abs=0.05)       # exp16 data-derived filament current density
    assert 0.78 < c["V_down"] < c["V_M"]


def test_split_law_limits(c):
    import law
    lc = law.lambda_c(10.0, c)
    assert np.isinf(law.w_split(0.99 * lc, 10.0, c))              # no split below the critical inhibition length
    assert law.w_split(100 * lc, 10.0, c) == pytest.approx(2 * lc, rel=0.02)   # minimum split width 2 lambda_c
    w = [law.w_split(x * lc, 10.0, c) for x in (1.1, 1.5, 2, 4)]
    assert all(a > b for a, b in zip(w, w[1:]))                   # monotonically decreasing towards 2 lambda_c


def test_centre_to_edge_drop_reaches_margin_at_split(c):
    """Direct check of the derivation with the discrete distributed-ballast rail solver of exp13."""
    import law
    sys.path.insert(0, str(EXP.parent / "exp12_filament"))
    import spatial as sp
    lam, r_D, dy = 60.0, 10.0, 0.25
    w = law.w_split(lam, r_D, c)
    n = int(round(2000 / dy))
    y = (np.arange(n) - n / 2 + 0.5) * dy
    jd = np.where(np.abs(y) < w / 2, c["j_on"], 0.0)
    V = sp.rail_voltage(jd, 1.0, r_D, dy, "distributed", r_D * lam**2)
    edge = np.interp(w / 2, y, V)
    assert edge - V[n // 2] == pytest.approx(c["dV_c"], rel=0.01)
