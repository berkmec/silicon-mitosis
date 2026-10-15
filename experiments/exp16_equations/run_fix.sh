#!/bin/sh
# Exp16 follow-ups, sequential (one Python process at a time).
cd /e/quantum_cell/experiments/exp16_equations || exit 1
PY=/e/quantum_cell/env/venv_sindy/Scripts/python.exe
LOG=/e/quantum_cell/results/exp16_equations
for s in s2c_filaments s2e_charging fig_pareto_gs s1b2_bsim4_followup; do
  echo "=== $s start $(date +%H:%M:%S)"
  PYTHONPATH= PYTHONIOENCODING=utf-8 "$PY" "$s.py" > "$LOG/log_$s.txt" 2>&1
  echo "=== $s exit $? $(date +%H:%M:%S)"
done
echo "ALL DONE"
