"""Sensitivity of the corrected dopant-potential Monte Carlo (exp6) to the depth clamp of the FEM correction factor.

exp6/run.py applies f = 1 + (1.8/max(z0, 1)) exp(-r^2/(2(z0+5)^2)). The FEM check only covers z0 = 10, 30, 60 nm, so the
factor below 10 nm is an extrapolation. This script re-executes exp6's sampling (same seeds, NS=80, grid) and changes
ONLY max(z0, 1.0) -> max(z0, ZCAP). ZCAP=1 must reproduce results/exp6/rdf.json (self-check).
Output: results/exp6/cap_sensitivity.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "experiments" / "exp6_rdf" / "run.py"
src = RUN.read_text(encoding="utf-8").split("# ---- checks")[0]
assert src.count("max(z0, 1.0)") == 1, "clamp expression not found exactly once"
src = src.replace("max(z0, 1.0)", "max(z0, ZCAP)")
G = {"__file__": str(RUN), "__name__": "exp6_sens", "ZCAP": 1.0, "NS": 80}
exec(compile(src, str(RUN), "exec"), G)
H = G["H"]
NS = G["NS"]


def puddles(U):
    """Same criteria as exp6/run.py: isolated classical puddles at three fill levels, bound if E1(inf well) < depth."""
    Umin = U.min()
    out = []
    for dE in (0.015, 0.030, 0.045):
        lab, nl = ndimage.label(U < Umin + dE)
        touch = set(np.unique(lab[0, :])) | set(np.unique(lab[-1, :]))
        for k in range(1, nl + 1):
            if k in touch:
                continue
            sel = lab == k
            a = float(np.sqrt(float(sel.sum()) * H * H))
            depth = dE - float(U[sel].min() - Umin)
            e1 = 2 * 0.0381 / 0.19 * np.pi ** 2 / a ** 2
            out.append(dict(strict=bool(e1 < depth), a=a))
    return out


def run(na, WD, zcap):
    G["ZCAP"] = zcap
    rng = np.random.default_rng(12345)
    sig, strict_samples = [], set()
    for s in range(NS):
        _, _, phi, _ = G["one_sample"](rng, na, WD, corr=True)
        U = -phi
        U = U - U.mean()
        sig.append(float(U.std()))
        if any(d["strict"] for d in puddles(U)):
            strict_samples.add(s)
    rng = np.random.default_rng(777)
    means = []
    for s in range(NS):
        _, _, phi, _ = G["one_sample"](rng, na, WD, corr=True)
        means.append(float(-phi.mean()))
    return dict(zcap_nm=zcap, sigma_psi_mV=float(np.mean(sig) * 1e3), sigma_VT_MC_mV=float(np.std(means) * 1e3),
                strict_frac=len(strict_samples) / NS)


if __name__ == "__main__":
    out = {}
    for label, na, WD in (("nominal_1.7e17", 1.7e17, 80.9), ("effective_3.9e17", 3.92e17, 54.6)):
        for zcap in (1.0, 10.0):
            t0 = time.time()
            r = run(na, WD, zcap)
            r["seconds"] = round(time.time() - t0, 1)
            out[f"{label}_ZCAP{zcap:g}"] = r
            print(label, r, flush=True)
    (ROOT / "results" / "exp6" / "cap_sensitivity.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("done", flush=True)
