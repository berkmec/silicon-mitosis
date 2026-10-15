"""Regression / validation tests for exp12 (needs the ngspice BSIM4 tables in experiments/exp12_filament/spice)."""
import sys
from pathlib import Path

import numpy as np
import pytest

EXP = Path(__file__).resolve().parents[1] / "experiments" / "exp12_filament"
pytestmark = pytest.mark.skipif(not (EXP / "spice" / "tab_vb0.0.dat").exists(), reason="ngspice tables not generated")
sys.path.insert(0, str(EXP))

KW = dict(vov0=-0.1, ss=0.02, fii=10, von=1.05, vs=0.005, g_sub=1e-13)


def test_single_device_latch_hysteresis():
    import reduced_model as rm

    V_up, V_down, _, _ = rm.sweep(np.round(np.arange(0.6, 3.0001, 0.02), 3), **KW)
    assert V_up == pytest.approx(1.56, abs=0.021) and V_down == pytest.approx(0.78, abs=0.021)


def test_uniform_spatial_reproduces_single_device():
    import spatial as sp

    r = sp.simulate(W=5.0, dy=0.5, R_L=0.0, mismatch=False, vdd_up=(0.6, 2.0), dvdd=0.02, R_sq=1e8)
    up = [h for h in r["hist"] if h["dir"] == "up"]
    dn = [h for h in r["hist"] if h["dir"] == "down"]
    assert next(h["vdd"] for h in up if h["frac_high"] > 0.5) == pytest.approx(1.56, abs=0.021)
    assert next(h["vdd"] for h in dn if h["frac_high"] < 0.5) == pytest.approx(0.78, abs=0.021)


def test_filament_clamps_at_maxwell_voltage():
    import spatial as sp
    from scipy.integrate import quad
    from scipy.optimize import brentq
    import reduced_model as rm

    def area(vds):
        r = rm.roots_at(vds, **KW)
        lo, hi = min(r), max(r)
        f = lambda v: float((rm.j_in(np.array([v]), vds, KW["vov0"], KW["ss"], KW["fii"])  # noqa: E731
                             - rm.j_out(np.array([v]), KW["von"], KW["vs"], KW["g_sub"]))[0])
        return quad(f, lo, hi, limit=400)[0]

    v_maxwell = brentq(area, 0.80, 0.82)
    r = sp.simulate(R_sq=1e12, dy=0.5, mismatch=True, vdd_up=(1.3, 1.7), dvdd=0.1)
    mid = [h for h in r["hist"] if h["dir"] == "up" and 0.1 < h["frac_high"] < 0.9]
    assert mid and all(abs(h["vds"] - v_maxwell) < 0.01 for h in mid)
