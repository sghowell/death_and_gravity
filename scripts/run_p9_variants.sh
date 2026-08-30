#!/bin/zsh
cd "$(dirname "$0")/.."
export PYTHONPATH=problems/P9/src
run() { uv run python -u -m p9.run_lkr_certified --L 1.5 --Delta 4 --refine 2 --workers 8 --passes 8 "$@" > "problems/P9/results/lkr_cert_variant_$(echo "$@" | tr -d ' -').log" 2>&1; }
run --sn union3
run --sn union3 --dv
run --dv
run --sn dessn5yr
run --sn dessn5yr --dv
uv run python -u -m p9.report_certified --certify-feasible > problems/P9/results/variants_report.log 2>&1
echo done > problems/P9/results/variants_queue_done.txt
