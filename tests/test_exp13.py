"""Validation of the exp13 resistive drain (local inhibitor) in experiments/exp12_filament/spatial.py."""
import sys
from pathlib import Path

import numpy as np
import pytest

EXP = Path(__file__).resolve().parents[1] / "experiments" / "exp12_filament"
pytestmark = pytest.mark.skipif(not (EXP / "spice" / "tab_vb0.0.dat").exists(), reason="ngspice tables not generated")
sys.path.insert(0, str(EXP))


def test_end_fed_rail_centre_drop_is_analytic():
    import spatial as sp

    jd, dy, W, r_D = np.full(200, 1e-6), 0.5, 100.0, 100.0
    V = sp.rail_voltage(jd, 1.0, r_D, dy, "both")
    assert 1.0 - V.min() == pytest.approx(r_D * jd.sum() * dy * W / 8, rel=1e-6)


def test_distributed_ballast_decay_length():
    import spatial as sp

    n, dy, r_D, lam = 8000, 0.05, 10.0, 20.0
    jd = np.zeros(n)
    jd[n // 2] = 1e-3 / dy
    V = sp.rail_voltage(jd, 1.0, r_D, dy, "distributed", r_D * lam**2)
    y = (np.arange(n) - n // 2) * dy
    sel = (y > 5) & (y < 60)
    assert -1.0 / np.polyfit(y[sel], np.log(1.0 - V[sel]), 1)[0] == pytest.approx(lam, rel=1e-3)


def test_zero_rail_resistance_is_global_coupling():
    import spatial as sp

    assert np.allclose(sp.rail_voltage(np.full(10, 1e-6), 0.9, 0.0, 0.5), 0.9)
