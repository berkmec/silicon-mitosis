"""Exp17: unbiased symbolic regression (gplearn genetic programming) of the split width w(lam, r_D) from the clean-protocol
data saved by symreg.py (results/exp17/symreg_xy.npy). Run with env/venv_sindy (gplearn is not in the main venv).
Several seeds; reports the expressions, their training error and whether they contain the log structure."""
import json
import sys
from pathlib import Path

import numpy as np
from gplearn.genetic import SymbolicRegressor

OUT = Path(__file__).resolve().parents[2] / "results" / "exp17"
d = np.load(OUT / "symreg_xy.npy")
X, y = d[:, :2], d[:, 2]
res = []
for seed in range(int(sys.argv[1]) if len(sys.argv) > 1 else 5):
    gp = SymbolicRegressor(population_size=3000, generations=40, function_set=("add", "sub", "mul", "div", "log", "sqrt"),
                           metric="mean absolute error", parsimony_coefficient=0.02, p_crossover=0.7, p_subtree_mutation=0.1,
                           p_hoist_mutation=0.05, p_point_mutation=0.1, max_samples=1.0, random_state=seed, n_jobs=1,
                           const_range=(-100.0, 100.0), init_depth=(2, 5), feature_names=["lam", "rD"])
    gp.fit(X, y)
    pred = gp.predict(X)
    res.append(dict(seed=seed, expr=str(gp._program), rel_rmse=float(np.sqrt(np.mean(((pred - y) / y) ** 2))),
                    length=int(gp._program.length_)))
    print(res[-1], flush=True)
(OUT / "symreg_gp.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
