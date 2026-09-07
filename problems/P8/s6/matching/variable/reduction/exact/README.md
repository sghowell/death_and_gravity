# S6.33: exact background, still a separate inverse problem

The unchanged parent's own-f equations have a unique positive isotropic
branch through b(0)=2 on the specified off-shell physical probe. A whole-
box contraction proof gives a positive lapse uniformly as c approaches 2
from above. The literal action remains undefined at c=2.

The exact tensor action distinguishes the prescribed-metric own-f inner
coefficient 16 from the coupled relative-mode coefficient 80. Constructing
the background does not select tensor homogeneous data or a Green inverse.

See [the statement](FORMULATION.md), [the Euler/uniqueness proof](notes/stationary.md),
[the rational whole-box proof](notes/intervals.md), and
[the independent-source boundary](notes/sources.md).
The report recursively rebuilds S6.32 and its unchanged ancestors.

## Ordinary replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/variable/reduction/exact/tests", "-q", "-p", "no:faulthandler"]))'
```

With every local P8 src directory on the Python path, the independent CLI
is `python -m p8_exact_stationary.verify --check`. It rebuilds and compares
the report read-only. Without --check it prints JSON and writes no files.

Ordinary report and test replays do not use the opt-in broad-regression
exact-GCD adapter. Both seeds are zero and the diagnosed host timer plugin
is disabled. Symbolic bridges clear nonzero denominators before expansion;
this is an exact polynomial check, not a numeric surrogate.

Original P8 remains open. No physical metric, matter coupling, full parent
background, momentum band or cutoff is reassigned by this checkpoint.
