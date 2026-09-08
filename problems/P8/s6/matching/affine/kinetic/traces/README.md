# S6.39: a rank-complete constant two-trace test

This child covers every constant real symmetric kinetic matrix for the
two specified projectively invariant distortion-trace curls, with the
exact original CD/M1 affine action otherwise unchanged. Every nonzero
matrix fails the literal parent-health test on the actual rolling
background. The zero matrix remains the original auxiliary theory.

The result includes off-diagonal mixing, both kinetic signs, and the
two exceptional null-Schur directions. Those directions require a fresh
full-connection constraint calculation; a singular Proca limit would
miss their unconstrained higher-derivative mode. Free chi is retained
throughout the metric/matter constraints.

See the [family and limits](FORMULATION.md),
[all-component geometry](notes/geometry.md),
[rank-one dynamics](notes/dynamics.md),
[rank-two proof](notes/rank-two.md), and
[independent interfaces](notes/interfaces.md).

This closes this constant two-trace family, not general metric-affine
parents or P8. No unhealthy low-frequency EFT band is asserted.

## Read-only replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/affine/kinetic/traces/tests", "-q", "-p", "no:faulthandler"]))'
```

The following standalone CLI also discovers the local P8 source
directories and sets both random seeds explicitly:

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'import sys, runpy; from pathlib import Path; from sympy.core.random import seed; seed(0); sys.path[:0]=[str(p.resolve()) for p in sorted(Path("problems/P8").rglob("src"))]; sys.argv=["p8_affine_traces.verify", "--check"]; runpy.run_module("p8_affine_traces.verify", run_name="__main__")'
```

Without `--check` it prints the report without writing files. The
ordinary tests and this CLI use no broad-regression GCD adapter.
The diagnosed host faulthandler timer plugin is disabled in pytest.
