"""Exp17: mechanism check at every observed split. For the filament configuration one step BEFORE each split (all
filaments, device ends included), compute the exact centre-to-edge drain drop of the parent with the discrete rail
(law.du_config, uniform j_on). If the split is a centre quench at V_M - V_down, all values cluster at dV_c.
Also compares the periodic-array closed form (law.du_array) for the low-mismatch triplets. -> results/exp17/mechanism.json"""
import glob
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import law  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp17"


def main():
    c0 = law.constants()
    cache = {}
    rows = []
    for f in sorted(glob.glob(str(OUT / "runs" / "*.json"))):
        r = json.loads(Path(f).read_text(encoding="utf-8"))
        j, h = r["job"], r["hist"]
        key = json.dumps(j.get("p"), sort_keys=True)
        if key not in cache:
            cache[key] = law.constants(j["p"]) if j.get("p") else c0
        c = cache[key]
        for k in range(1, len(h)):
            prev, cur = h[k - 1], h[k]
            P = [(cc - w / 2, cc + w / 2) for cc, w in zip(prev["centres_um"], prev["widths_um"])]
            Cn = [(cc - w / 2, cc + w / 2) for cc, w in zip(cur["centres_um"], cur["widths_um"])]
            for p in P:
                kids = [q for q in Cn if min(p[1], q[1]) - max(p[0], q[0]) > 0]
                if len(kids) >= 2:
                    du = law.du_config(P, p, j["W"], j["lam"], j["r_D"], c)
                    du_iso = law.du_array(p[1] - p[0], 1e9, j["lam"], j["r_D"], c)
                    if j.get("p") and j["p"]["g_sub"] < 1e-14:
                        prot = "düşük sızıntı (g 1e-15)"
                    elif j.get("p"):
                        prot = "evrensellik (U)"
                    else:
                        prot = "hızlı (R_L 2k)" if j["R_L"] < 4e3 else "temiz (R_L 5k)"
                    rows.append(dict(run=Path(f).stem, protocol=prot,
                                     mm=j["sigma_scale"], lam=j["lam"], r_D=j["r_D"], vdd=cur["vdd"], w=p[1] - p[0],
                                     n_before=len(P), du_mV=du * 1e3, du_isolated_mV=du_iso * 1e3,
                                     ratio=du / c["dV_c"]))
    out = dict(dV_c_mV=c0["dV_c"] * 1e3, events=rows)
    for prot in ("temiz (R_L 5k)", "hızlı (R_L 2k)", "evrensellik (U)", "düşük sızıntı (g 1e-15)"):
        for mm in (1.0, 0.1):
            x = [e["ratio"] for e in rows if e["protocol"] == prot and e["mm"] == mm]
            if x:
                out[f"{prot}, mm={mm}"] = dict(n=len(x), mean=float(np.mean(x)), sd=float(np.std(x)),
                                               min=float(np.min(x)), max=float(np.max(x)))
    (OUT / "mechanism.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for e in rows:
        print(f"{e['run']:16s} {e['protocol']:15s} mm={e['mm']:.1f} lam={e['lam']:5.0f} rD={e['r_D']:4.0f} V_DD={e['vdd']:.2f} "
              f"w={e['w']:6.1f} n={e['n_before']} du={e['du_mV']:5.2f} mV (yalıtık formül {e['du_isolated_mV']:5.2f}) "
              f"oran={e['ratio']:.3f}")
    for k, v in out.items():
        if isinstance(v, dict):
            print(k, {kk: round(vv, 3) for kk, vv in v.items()})


if __name__ == "__main__":
    main()
