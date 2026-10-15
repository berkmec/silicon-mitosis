#!/bin/sh
# Exp16: run the remaining analyses ONE AFTER ANOTHER (single Python process at a time, each sets below-normal priority).
cd /e/quantum_cell/experiments/exp16_equations || exit 1
PY=/e/quantum_cell/env/venv_sindy/Scripts/python.exe
LOG=/e/quantum_cell/results/exp16_equations
for s in s1c_exp10_hf s1d_exp1 s1b_bsim4 s1a2_weak_reduced s2e_charging s2c_filaments s2d_scaling s2b_gs_count s2a_droplets; do
  echo "=== $s start $(date +%H:%M:%S)"
  PYTHONPATH= PYTHONIOENCODING=utf-8 "$PY" "$s.py" > "$LOG/log_$s.txt" 2>&1
  echo "=== $s exit $? $(date +%H:%M:%S)"
done
echo "ALL DONE"
