"""Exp16: Pareto plot of the Gray-Scott PDE-FIND threshold paths on a log scale (from s1a_grayscott.json; no refit).
Output: results/exp16_equations/fig_s1a_pareto.png"""
import common as C  # noqa: I001

import json

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

d = json.loads((C.OUT / "s1a_grayscott.json").read_text(encoding="utf-8"))
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
for j, ax in enumerate(axs):
    for v in d["variants"]:
        nt = [p["n_terms"][j] for p in v["path"]]
        e = [max(1 - p["r2_val"][j], 1e-16) for p in v["path"]]
        ax.loglog(nt, e, "o-", ms=3, lw=1, label=v["name"])
        ex = [(p["n_terms"][j], max(1 - p["r2_val"][j], 1e-16)) for p in v["path"] if p["exact"][j]]
        if ex:
            ax.plot(*zip(*ex), "k*", ms=9, mfc="none")
    ax.axvline([5, 4][j], color="k", ls=":", lw=0.8)
    ax.set_xlabel(["u_t", "v_t"][j] + " denklemindeki terim sayısı")
    ax.set_ylabel("1 − R² (ayrı doğrulama koşusu)")
    ax.set_title(f"{['u', 'v'][j]} denklemi: eşik yolu (★ = tam doğru destek, dikey çizgi = gerçek terim sayısı)",
                 fontsize=8)
axs[0].legend(fontsize=6)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s1a_pareto.png", dpi=110)
print("ok")
