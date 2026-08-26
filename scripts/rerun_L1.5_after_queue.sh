#!/bin/zsh
cd "$(dirname "$0")/.."
export PYTHONPATH=problems/P9/src
while [ ! -f problems/P9/results/curve_queue_done.txt ]; do sleep 120; done
uv run python -u -m p9.run_lkr_certified --L 1.5 --Delta 4 --refine 2 --workers 8 --passes 8 > problems/P9/results/lkr_cert_L1.5_D4_r2_v1.log 2>&1
echo done > problems/P9/results/rerun_L1.5_done.txt
