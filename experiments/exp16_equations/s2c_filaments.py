"""Exp16 step 2c: filament laws in exp13 (SKY130 4 K reduced model, current filaments along the channel width).

For every run (results/exp13/A*,B*,C*,D*.json) and sweep direction, points with >= 1 filament and total filament width
below the channel width W are used: V_DD, I, V_ds, w_total.
S-NDC + load line expectation: V_ds clamps at a Maxwell-type voltage V_M and I = j_f w + j_0 (W - w) ->
w_total = (I - j_0 W)/(j_f - j_0), i.e. linear in I (and in V_DD through I = (V_DD - V_M)/R_L).
Sparse regression (PySINDy STLSQ) of w_total on {1, I, I^2, sqrt I}: trained on the up sweep, validated on the down
sweep (hysteresis branch, independent states) when it has >= 5 points; V_ds(I) on {1, I, I^2}.
Division (C runs, distributed ballast): width at the split vs inhibition length lambda.
Output: results/exp16_equations/s2c_filaments.json, fig_s2c_filaments.png."""
import common as C  # noqa: I001

import glob
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C.low_priority()
NAMES = ["1", "I", "I^2", "sqrt I"]


def lib(I):
    return np.column_stack([np.ones_like(I), I, I**2, np.sqrt(I)])


TAUS = np.logspace(-5, -0.5, 20)
runs = {}
fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
for f in sorted(glob.glob(str(C.RESULTS / "exp13" / "[ABCD]*.json"))):
    name = Path(f).stem
    d = json.loads(Path(f).read_text(encoding="utf-8"))
    if "result" not in d:  # *_fine.json files store the result dict at top level
        d = {"result": d}
    W = d["result"]["params"]["W"]
    pts = {"up": [], "down": []}
    split = None
    prev = None
    for h in d["result"]["hist"]:
        wt = float(np.sum(h["widths_um"])) if h["n_fil"] else 0.0
        if h["dir"] == "up" and h["splits"] > 0 and prev is not None and split is None:
            split = dict(vdd=h["vdd"], w_before=prev, daughters=h["widths_um"])
        if h["dir"] == "up":
            prev = wt
        if h["n_fil"] >= 1 and 2.0 < wt < W - 1.0:
            pts[h["dir"]].append((h["vdd"], h["I_uA"], h["vds"], wt, h["n_fil"]))
    up, dn = np.array(pts["up"]), np.array(pts["down"])
    res = dict(W=W, lambda_um=d["result"]["params"].get("lambda_um"), R_L=d["result"]["params"]["R_L"],
               r_D=d["result"]["params"].get("r_D"), contacts=d["result"]["params"].get("contacts"),
               n_up=len(up), n_down=len(dn), split=split,
               max_single_width_up=float(max([p[3] for p in pts["up"] if p[4] == 1], default=0.0)))
    if len(up) >= 6:
        if len(dn) >= 5:
            tr, va = up, dn
            res["validation"] = "aşağı tarama (histerezis dalı)"
        else:
            tr, va = up[::2], up[1::2]
            res["validation"] = "yukarı taramada araya serpiştirilmiş yarı"
        path = C.threshold_path(lib(tr[:, 1]), tr[:, 3], lib(va[:, 1]), va[:, 3], TAUS)
        sel = C.select_model_r2(path, r2_tol=1e-3)
        cf = sel["coef"][0]
        lin = np.polyfit(up[:, 1], up[:, 3], 1)
        vs = np.polyfit(up[:, 1], up[:, 2], 1)
        res.update(support=[NAMES[i] for i in np.flatnonzero(cf)], coef={NAMES[i]: float(cf[i]) for i in np.flatnonzero(cf)},
                   r2_val=sel["r2_val"][0], tau=sel["tau"],
                   linear_fit_up=dict(slope_um_per_uA=float(lin[0]), intercept_um=float(lin[1]),
                                      j_f_uA_per_um=float(1 / lin[0]), I_at_w0_uA=float(-lin[1] / lin[0])),
                   vds_up=dict(median=float(np.median(up[:, 2])), p10=float(np.percentile(up[:, 2], 10)),
                               p90=float(np.percentile(up[:, 2], 90)), slope_V_per_uA=float(vs[0]),
                               change_over_range_V=float(vs[0] * np.ptp(up[:, 1]))),
                   vdd_linear=dict(zip(["slope_um_per_V", "intercept_um"],
                                       map(float, np.polyfit(up[:, 0], up[:, 3], 1)))))
        if name.startswith(("A", "C_lam40_RL5000", "D")) and not name.endswith("fine"):
            ax = axs[0] if name.startswith("A") else axs[1]
            ax.plot(up[:, 1], up[:, 3], ".", ms=3, label=name)
            if len(dn):
                ax.plot(dn[:, 1], dn[:, 3], "x", ms=3)
            axs[2].plot(up[:, 1], up[:, 2], "-", lw=0.8, label=name)
    runs[name] = res
    print(name, {k: res.get(k) for k in ("n_up", "n_down", "support", "r2_val", "split", "max_single_width_up")},
          res.get("linear_fit_up"), res.get("vds_up"), flush=True)
# division vs lambda
div = {n: dict(lam=r["lambda_um"], R_L=r["R_L"], split=r["split"], max_single_width=r["max_single_width_up"])
       for n, r in runs.items() if n.startswith("C")}
rob = json.loads((C.RESULTS / "exp13" / "split_robustness.json").read_text(encoding="utf-8"))
out = dict(runs=runs, division_table=div, split_robustness_events={k: v["event_steps"] for k, v in rob.items()},
           exp12_maxwell_V=json.loads((C.RESULTS / "exp12" / "maxwell.json").read_text(encoding="utf-8"))["V_ds_maxwell"])
C.save_json("s2c_filaments.json", out)
axs[0].set_title("global bağlaşım (A): w_toplam – I (• yukarı, x aşağı)", fontsize=9)
axs[1].set_title("dağıtık balast λ=40 µm ve D (R□=1e10)", fontsize=9)
for ax in axs[:2]:
    ax.set_xlabel("I (µA)")
    ax.set_ylabel("toplam filaman genişliği (µm)")
    ax.legend(fontsize=6)
axs[2].set_xlabel("I (µA)")
axs[2].set_ylabel("V_ds (V)")
axs[2].set_title("kenetlenme: V_ds neredeyse sabit", fontsize=9)
axs[2].legend(fontsize=5)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s2c_filaments.png", dpi=110)
