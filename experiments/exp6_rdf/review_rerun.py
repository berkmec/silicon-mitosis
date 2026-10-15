"""Review rerun (2026-10-11) of the exp6 accidental-dot statistics with two fixes found by the audit:
 (1) the deterministic 'dome' E[U](x,y) (dopants only inside the box, no S/D electrodes) is subtracted, so puddles are
     defined by the RANDOM dopant potential only (review_campbell.py: dome = 55-59 % of sigma_psi^2);
 (2) a dot must not touch ANY box edge (the original run excluded only the source/drain x-edges, not the STI y-edges).
Same kernel (image model, nimg = 8), same seeds and samples as run.py; FEM-corrected variants are dropped (their fitted
factor exceeds the physical r=0 limit 2 eps_Si/(eps_Si+eps_ox) = 1.5). Also reports the original rule on the same samples
and Clopper-Pearson 95 % intervals."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from review_campbell import mean_potential, kernel_z, LEFF, WCH, H  # noqa: E402,F401


def surface_potential(xs, ys, dop, WD, nimg=8):
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    phi = np.zeros_like(gx)
    for (x0, y0, z0) in dop:
        phi += kernel_z(gx - x0, gy - y0, z0, WD, nimg)
    return phi


def cp(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return [float(lo), float(hi)]


def dots_in(U, edges_all):
    found = []
    Umin = U.min()
    for dE in (0.015, 0.030, 0.045):
        lab, nl = ndimage.label(U < Umin + dE)
        touch = set(np.unique(lab[0, :])) | set(np.unique(lab[-1, :]))
        if edges_all:
            touch |= set(np.unique(lab[:, 0])) | set(np.unique(lab[:, -1]))
        for k in range(1, nl + 1):
            if k in touch:
                continue
            sel = lab == k
            area = float(sel.sum()) * H * H
            a_nm = float(np.sqrt(area))
            depth = dE - float(U[sel].min() - Umin)
            e1 = 2 * 0.0381 / 0.19 * np.pi**2 / a_nm**2
            pts = np.argwhere(sel).astype(float)
            asp = float(np.sqrt(np.max(np.linalg.eigvalsh(np.cov(pts.T))) / max(np.min(np.linalg.eigvalsh(np.cov(pts.T))), 1e-9))) if len(pts) > 2 else 1.0
            found.append(dict(dE=dE, a_nm=a_nm, depth_eV=depth, strict=bool(e1 < depth), opt=bool(0.5 * e1 < depth),
                              aspect=asp))
    return found


out = {}
NS = 80
for label, na, WD in (("nominal_1.7e17", 1.7e17, 80.9), ("effective_3.9e17", 3.92e17, 54.6)):
    dome = mean_potential(na, WD)  # E[U] on the same 2 nm grid inside the box
    rng = np.random.default_rng(12345)
    xs, ys = np.arange(H / 2, LEFF, H), np.arange(H / 2, WCH, H)
    tally = {"orijinal kural (x-kenarları, kubbe dahil)": [0, 0], "tüm kenarlar, kubbe dahil": [0, 0],
             "tüm kenarlar, kubbe çıkarılmış": [0, 0]}
    sizes, aspects = [], []
    for s in range(NS):
        n = rng.poisson(na * (LEFF * WCH * WD * 1e-21))
        dop = np.column_stack([rng.uniform(0, LEFF, n), rng.uniform(0, WCH, n), rng.uniform(0, WD, n)])
        U = -surface_potential(xs, ys, dop, WD)
        variants = ((list(tally)[0], U - U.mean(), False), (list(tally)[1], U - U.mean(), True),
                    (list(tally)[2], (U - dome) - (U - dome).mean(), True))
        for name, UU, alledges in variants:
            ds = dots_in(UU, alledges)
            tally[name][0] += any(d["strict"] for d in ds)
            tally[name][1] += any(d["opt"] for d in ds)
            if name == list(tally)[2]:
                sizes += [d["a_nm"] for d in ds if d["strict"]]
                aspects += [d["aspect"] for d in ds if d["strict"]]
        if s % 20 == 19:
            print(label, s + 1, {k: v for k, v in tally.items()}, flush=True)
    res = {k: dict(strict=f"{v[0]}/{NS}", strict_pct=100 * v[0] / NS, strict_CI95_pct=[100 * x for x in cp(v[0], NS)],
                   optimistic=f"{v[1]}/{NS}", optimistic_pct=100 * v[1] / NS,
                   optimistic_CI95_pct=[100 * x for x in cp(v[1], NS)]) for k, v in tally.items()}
    res["kubbe çıkarılmış, katı noktalar"] = dict(n=len(sizes), a_median_nm=float(np.median(sizes)) if sizes else None,
                                                   aspect_median=float(np.median(aspects)) if aspects else None)
    res["sigma_random_rms_mV"] = None
    out[label] = res
    print(label, json.dumps(res, ensure_ascii=False, indent=1), flush=True)
(ROOT / "results" / "exp6" / "review_rerun.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                             encoding="utf-8")
