"""Exp17 2D figures from the saved snapshots in results/exp17/runs2d/*.npz (keys V_<V_DD>: body potential, VD_<V_DD>: local
drain voltage).
- fig_2d_spot.png: V64_small (walled 420 um field, lambda = 64 um), spot born in a corner; top row V_b, bottom row V_D.
- fig_2d_periodic.png: P64 (periodic 540 um field, lambda = 64 um, designed nucleation site at the centre); V_b only."""
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "exp17"


def snapshots(name, vdds, with_vd, title, fname):
    z = np.load(OUT / "runs2d" / f"{name}.npz")
    rows = 2 if with_vd else 1
    fig, axs = plt.subplots(rows, len(vdds), figsize=(1.8 * len(vdds), 1.9 * rows + 0.3), squeeze=False)
    for i, v in enumerate(vdds):
        axs[0, i].imshow(z[f"V_{v:.2f}"].astype(float).T, origin="lower", cmap="magma", vmin=0, vmax=1)
        axs[0, i].set_title(f"V_DD={v:.2f}", fontsize=7)
        if with_vd:
            axs[1, i].imshow(z[f"VD_{v:.2f}"].astype(float).T, origin="lower", cmap="viridis")
    for a in axs.ravel():
        a.set_xticks([])
        a.set_yticks([])
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / fname, dpi=90)
    plt.close(fig)
    print("ok", OUT / fname)


if __name__ == "__main__":
    snapshots("V64_small", (1.8, 3.3, 3.8, 5.8, 8.8, 9.3, 10.8, 11.8), True,
              "2B dizi (λ = 64 µm, 420 µm alan): leke büyüyor → uzayıp şeride dönüşüyor → bükülüyor (bölünme yok)",
              "fig_2d_spot.png")
    snapshots("P64_ckpt", (3.2, 6.2, 8.2, 9.2, 9.7, 10.7, 11.7), False,
              "2B, periyodik (duvarsız), merkezde tasarlanmış çekirdek, λ = 64 µm: yuvarlak leke R_c'yi geçiyor → uzayıp şerit "
              "oluyor (bölünme yok)", "fig_2d_periodic.png")
