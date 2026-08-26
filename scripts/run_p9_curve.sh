#!/bin/zsh
# Sequential certified LKR runs for the P9 curve (each ~1.5-2 h with 8 workers on the x4 grid).
# Waits for any running p9.run_lkr_certified process to finish first.
cd "$(dirname "$0")/.."
export PYTHONPATH=problems/P9/src
while pgrep -f "p9.run_lkr_certified" >/dev/null; do sleep 60; done
for L in 1 2 3 5 10; do
  uv run python -u -m p9.run_lkr_certified --L $L --Delta 4 --refine 2 --workers 8 --passes 8 > problems/P9/results/lkr_cert_L${L}_D4_r2.log 2>&1
done
for D in 1 9; do
  uv run python -u -m p9.run_lkr_certified --L 1.5 --Delta $D --refine 2 --workers 8 --passes 8 > problems/P9/results/lkr_cert_L1.5_D${D}_r2.log 2>&1
done
echo "P9 curve queue finished" > problems/P9/results/curve_queue_done.txt
