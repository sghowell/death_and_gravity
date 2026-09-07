# S6.35: own-f mixing survives on a fixed physical-time window

For the unchanged exact own-f branch, a uniform bound on the complete
canonical tensor potential proves a finite-window homogeneous transfer
error below 1/80 and a surviving mixing coefficient above 1/80. Unlike
the sibling causal-pulse result, the physical duration does not shrink
with delta. Incoming hidden data are supplied; no zero-data source or
physical EFT conclusion is inferred.

See [the exact statement](FORMULATION.md), [the potential proof](notes/potential.md),
[the connection and propagation proof](notes/connection.md),
[the action/domain interfaces](notes/interfaces.md), and
[the primary-source audit](notes/sources.md).

The primary Fraction partial-jet calculation and an independently written
Arb total-jet calculation cover the full parameter box. Acb special-function
checks corroborate the analytic connection; numerical residual containment
does not replace the written proof. The report recursively rebuilds S6.33
and its frozen ancestry without changing them.

## Ordinary replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/variable/reduction/exact/scattering/tests", "-q", "-p", "no:faulthandler"]))'
```

With local P8 src directories on the Python path, the separate CLI is
`python -m p8_own_scattering.verify --check`. Without --check it prints the
computed report and writes no files. The ordinary checks use no broad-
regression exact-GCD adapter. Both seeds are zero and the diagnosed host
timer plugin is disabled.

The free exterior is a mathematical endpoint convention, not a physical
vacuum. The physical metric and clock remain the original off-shell probe.
Original P8 remains open.
