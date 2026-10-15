"""Mechanism figure for the README (docs/figures/mechanism.png).
(a) device schematic, (b) drain-voltage profile of a filament and the centre-quench criterion,
(c) a lineage from a simulation (results/exp17/runs/H_lam60_s1.json), (d) universal split law with the clean first splits."""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))
import law  # noqa: E402
import spatial as sp  # noqa: E402

c = law.constants()
fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])


def zigzag(ax, x0, y0, x1, y1, n=6, amp=0.12, **kw):
    t = np.linspace(0, 1, 2 * n + 1)
    x, y = x0 + (x1 - x0) * t, y0 + (y1 - y0) * t
    dx, dy = x1 - x0, y1 - y0
    L = np.hypot(dx, dy)
    nx_, ny_ = -dy / L, dx / L
    off = amp * np.array([0] + [(-1) ** i for i in range(2 * n - 1)] + [0])
    ax.plot(x + nx_ * off, y + ny_ * off, **kw)


# (a) schematic
ax = fig.add_subplot(gs[0, 0])
W, nseg = 10.0, 10
ax.add_patch(Rectangle((0, 0), W, 1.0, fc="#d9d9d9", ec="k"))                    # source (n+)
ax.add_patch(Rectangle((0, 1.0), W, 0.45, fc="#4a6fa5", ec="k"))                 # gate over the channel
ax.add_patch(Rectangle((0, 1.45), W, 0.9, fc="#d9d9d9", ec="k"))                 # drain (n+)
ax.add_patch(Rectangle((3.6, 1.0), 2.8, 0.45, fc="#ff8c1a", ec="none", alpha=0.9))  # current filament
ax.text(5.0, 1.22, "filament (on)", ha="center", va="center", fontsize=9, color="k")
ax.text(W + 0.15, 0.5, "source (0 V)", va="center", fontsize=9)
ax.text(W + 0.15, 1.22, "gate: V_GS ≈ V_T − 0.1 V", va="center", fontsize=9)
ax.text(W + 0.15, 1.9, "drain", va="center", fontsize=9)
xs = (np.arange(nseg) + 0.5) * W / nseg
yr, yb = 2.9, 4.4
for i, x in enumerate(xs):
    ax.plot([x, x], [2.35, yr], "k-", lw=1)                                      # drain contact of segment i
    ax.plot(x, yr, "ko", ms=4)
    zigzag(ax, x, yr, x, yb, n=4, amp=0.12, c="C3", lw=1.3)                       # ballast r_v
    if i < nseg - 1:
        zigzag(ax, x, yr, xs[i + 1], yr, n=3, amp=0.08, c="C0", lw=1.3)          # rail r_D
ax.plot([xs[0], xs[-1]], [yb, yb], "k-", lw=2)
zigzag(ax, W / 2, yb, W / 2, 5.6, n=4, amp=0.15, c="k", lw=1.5)
ax.text(W / 2 + 0.3, 5.0, "R_L", fontsize=10)
ax.text(W / 2 - 0.4, 5.85, "V_DD", fontsize=10, ha="center")
ax.text(xs[-1] + 0.2, yb, "common node V_c", va="center", fontsize=9)
ax.text(xs[-1] + 0.2, yr, "drain rail  r_D (Ω/µm)", va="center", fontsize=9, color="C0")
ax.text(xs[-1] + 0.2, (yr + yb) / 2, "ballast  r_v (Ω·µm)", va="center", fontsize=9, color="C3")
ax.add_patch(FancyArrowPatch((0, -0.45), (W, -0.45), arrowstyle="<->", mutation_scale=12))
ax.text(W / 2, -0.8, "channel width y  (W ≈ 300–600 µm, L = 0.15 µm)", ha="center", fontsize=9)
ax.add_patch(FancyArrowPatch((3.6, 0.85), (6.4, 0.85), arrowstyle="<->", mutation_scale=10, color="C1"))
ax.text(5.0, 0.62, "w", ha="center", color="C1", fontsize=10)
ax.set_xlim(-0.3, W + 3.3)
ax.set_ylim(-1.1, 6.3)
ax.axis("off")
ax.set_title("(a) Device: wide SKY130 nFET with a segmented, ballasted drain\n"
             "inhibition length λ = √(r_v / r_D)", fontsize=11)

# (b) drain-voltage profile and centre quench
ax = fig.add_subplot(gs[0, 1])
lam, r_D, dy = 60.0, 10.0, 0.25
ws = law.w_split(lam, r_D, c)
y = np.arange(-300, 300, dy) + dy / 2
for w, col, lab in ((0.6 * ws, "C0", "w = 0.6 w_s (stable)"), (ws, "C3", f"w = w_s = {ws:.0f} µm (split)")):
    jd = np.where(np.abs(y) < w / 2, c["j_on"], 0.0)
    V = sp.rail_voltage(jd, 1.0, r_D, dy, "distributed", r_D * lam**2)
    V = V - np.interp(w / 2, y, V) + c["V_M"]                                      # fronts sit at V_M
    ax.plot(y, V * 1e3, c=col, lw=2, label=lab)
    ax.axvspan(-w / 2, w / 2, color=col, alpha=0.07)
ax.axhline(c["V_M"] * 1e3, c="k", ls="--", lw=1)
ax.axhline(c["V_down"] * 1e3, c="C3", ls=":", lw=1.5)
ax.text(-290, c["V_M"] * 1e3 + 1.2, "V_M (front / Maxwell voltage)", fontsize=9)
ax.text(-290, c["V_down"] * 1e3 - 3.0, "V_down (high state ceases to exist)", fontsize=9, color="C3")
ax.annotate("", xy=(0, c["V_down"] * 1e3), xytext=(0, c["V_M"] * 1e3),
            arrowprops=dict(arrowstyle="<->", color="k"))
ax.text(60, (c["V_M"] + c["V_down"]) * 500 - 2, f"Δu = V_M − V_down = {c['dV_c']*1e3:.2f} mV", fontsize=9)
ax.set_xlabel("y (µm)")
ax.set_ylabel("local drain voltage V_D (mV)")
ax.set_ylim(c["V_down"] * 1e3 - 8, c["V_M"] * 1e3 + 40)
ax.legend(fontsize=9, loc="upper right")
ax.set_title("(b) The centre is the most inhibited point: the filament splits when its\n"
             "centre-to-edge drain drop reaches V_M − V_down  (λ = 60 µm, r_D = 10 Ω/µm)", fontsize=11)

# (c) a lineage from a simulation (clean protocol, reference device)
ax = fig.add_subplot(gs[1, 0])
RUN_C = "H_lam60_s1"
r = json.loads((ROOT / "results" / "exp17" / "runs" / f"{RUN_C}.json").read_text(encoding="utf-8"))
for h in r["hist"]:
    for cc, w in zip(h["centres_um"], h["widths_um"]):
        ax.plot([h["I_uA"], h["I_uA"]], [cc - w / 2, cc + w / 2], c="#ff8c1a", lw=1.3)
mech = json.loads((ROOT / "results" / "exp17" / "mechanism.json").read_text(encoding="utf-8"))
for e in r["splits"]:
    q = next(m["ratio"] for m in mech["events"] if m["run"] == RUN_C and abs(m["vdd"] - e["vdd"]) < 1e-6)
    ax.plot(e["I_uA"], e["parent_c"], "k*", ms=13)
    ax.annotate(f"split at w = {e['parent_w']:.0f} µm\nΔu / (V_M − V_down) = {q:.2f}", (e["I_uA"], e["parent_c"]),
                textcoords="offset points", xytext=(-150, 12) if e["I_uA"] > 800 else (8, 10), fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
ax.set_xlabel("total drain current I (µA)")
ax.set_ylabel("position along the width y (µm)")
ax.set_ylim(0, 600)
ax.set_title(f"(c) A lineage (simulation, λ = 60 µm, W = 600 µm): one filament → 2 → 4\n"
             f"(isolated-filament prediction w_s = {ws:.0f} µm; total width = I / j_on, j_on = {c['j_on']*1e6:.2f} µA/µm)",
             fontsize=11)

# (d) universal curve
ax = fig.add_subplot(gs[1, 1])
an = json.loads((ROOT / "results" / "exp17" / "analysis.json").read_text(encoding="utf-8"))
x = np.linspace(1.02, 5, 400)
ax.plot(x, 2 * x * np.log(x / (x - 1)), "k-", lw=2, label="theory: w_s/λ_c = 2x ln[x/(x−1)],  x = λ/λ_c")
mk = {"base": ("o", "C0", "reference device set (r_D = 5, 10, 20 Ω/µm)"), "U1": ("s", "C1", "device set U1 (r_D = 10 Ω/µm)"),
      "U2": ("D", "C2", "device set U2 (r_D = 10 Ω/µm)")}
seen = set()
for t in an["table"]:
    if t["w_obs"] is None or t["lam_c"] is None or t["set"] == "N" or t["R_L"] < 4e3 or t["mm"] != 1.0:
        continue
    if t["n_before"] is None or t["n_before"] > 2:
        continue
    m, col, lab = mk[t["set"]]
    ax.plot(t["lam"] / t["lam_c"], t["w_obs"] / t["lam_c"], m, c=col, ms=8, mec="k",
            label=None if t["set"] in seen else lab)
    seen.add(t["set"])
ax.axhline(2, c="0.6", ls=":", lw=1)
ax.text(3.55, 2.08, "minimum: w_s → 2λ_c", fontsize=8, color="0.4")
ax.axvline(1, c="0.6", ls=":", lw=1)
ax.text(1.12, 6.3, "no division\nfor λ ≤ λ_c", fontsize=8, color="0.4")
ax.set_xlim(0.9, 4.6)
ax.set_ylim(0, 7)
ax.set_xlabel("λ / λ_c")
ax.set_ylabel("split width  w_s / λ_c")
ax.set_title("(d) The 14 clean first splits (three device sets) against one parameter-free curve\n"
             "λ_c = √(2(V_M − V_down) / (r_D j_on)); constants from the device model", fontsize=11)
ax.legend(fontsize=9)
fig.suptitle("Silicon mitosis: how a current filament in a ballasted SKY130 nFET divides", fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = ROOT / "docs" / "figures"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "mechanism.png", dpi=110)
print("ok", out / "mechanism.png")
