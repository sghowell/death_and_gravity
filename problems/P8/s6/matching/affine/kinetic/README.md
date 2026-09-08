# S6.38: kinetic terms need their own constraint test

This child tests three explicit additions to the exact auxiliary
CD/M1 affine lift. A homothetic-curvature square only adds a Maxwell
spectator. The simple trace-free curvature matrix-trace square has
opposite signs in two independent physical connection sectors.

A more selective projectively invariant quotient-vector curl preserves
the actual rolling background, but its full coupled scalar constraints
give one negative physical scalar kinetic direction at sufficiently
large momentum on both sides of the bounce. Its isolated Proca block
would miss this failure. The zero-coupling limit changes constraint rank.
See the [fixed scope](FORMULATION.md), [curvature screens](notes/curvature.md),
[vector Schur reduction](notes/vector.md), and
[coupled scalar test](notes/scalar.md).

No general affine no-go or low-frequency EFT cutoff follows from these
specified-action tests. The original P8 UV problem remains open.

## Read-only replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/affine/kinetic/tests", "-q", "-p", "no:faulthandler"]))'
```

With the local P8 source directories on the Python path, run
`python -m p8_affine_kinetic.verify --check`. Without `--check` it prints
the report without writing files. Ordinary tests and this CLI use no
broad-regression GCD adapter. Both random seeds are zero and the
diagnosed host faulthandler timer plugin is disabled.
