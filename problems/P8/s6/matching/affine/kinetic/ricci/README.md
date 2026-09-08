# S6.40: a curved auxiliary deformation, not a new heavy sector

This new Ricci-difference square is outside the closed two-trace curl
family. Its flat connection Schur symbol is null, but the curved
commutator survives. Retaining it and every original source component
gives an exact quadratic addition q*A*n² after all connection equations.

For positive coupling, the full scalar kinetic and principal gradient
matrices are positive on the stated regular charts. Negative coupling
has a regular-tail high-momentum kinetic obstruction. Zero coupling
returns the original auxiliary theory. The original rolling background
is preserved, but a nonzero regular deformation changes its physical
quadratic action. No new heavy spectrum or UV completion follows.

See the [literal scope](FORMULATION.md),
[all-64 curved connection proof](notes/geometry.md), and
[source, constraints, principal signs and kinetic bound](notes/dynamics.md).

## Read-only replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/affine/kinetic/ricci/tests", "-q", "-p", "no:faulthandler"]))'
```

The standalone CLI discovers all local P8 sources and seeds both RNGs:

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'import sys, runpy; from pathlib import Path; from sympy.core.random import seed; seed(0); sys.path[:0]=[str(p.resolve()) for p in sorted(Path("problems/P8").rglob("src"))]; sys.argv=["p8_affine_ricci.verify", "--check"]; runpy.run_module("p8_affine_ricci.verify", run_name="__main__")'
```

Without `--check` the CLI prints the report without writing files.
Neither ordinary command uses the broad-regression GCD adapter.
The diagnosed host faulthandler timer plugin is disabled in pytest.
