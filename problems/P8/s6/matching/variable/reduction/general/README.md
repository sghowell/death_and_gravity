# S6.31: regular roots and next-order freedom

The general ordinary two-metric interaction does not repair the retained-
order CD mismatch by choosing beta2 or beta3. On a positive simple
proportional root, the complete formal action through four derivatives
has the same specified scalar-tensor relation A1=2F2_X, with the separate
curvature-square operator retained.

This is not an all-order no-go. An explicit smooth coefficient deformation
leaves every retained-order datum unchanged but changes the formal
degree-six action. A compact conformal variation proves that the change
is genuinely non-boundary. An exact potential-only operator table and
an independent Euler boundary control show why a center-first projection
does not determine the full higher-order CD coefficient.

Read [the frozen statement](FORMULATION.md),
[the general stationary proof](notes/general-stationary.md),
[the representative control](notes/basis.md), and
[the source boundary](notes/sources.md).

## Replay

Discover local P8 source packages through the test conftest or explicitly
add their src directories to the Python path. With the repository's
environment, run the ordinary tests:

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/variable/reduction/general/tests", "-q", "-p", "no:faulthandler"]))'
```

The separate ordinary replay entry point is
`python -m p8_general_reduction.verify --check` with the local src paths
available. It rebuilds S6.30 and its ancestry before exact report
comparison. Certificate generation prints JSON to standard output and
does not overwrite frozen files. No floating-point fit or cached parent
verdict substitutes for the exact algebra and written proof.

Both random seeds are fixed for reproducibility. The faulthandler timer
plugin is disabled because of the previously diagnosed host behavior;
ordinary certificate checks do not use the broad-regression GCD adapter.

Original P8, a controlled higher-order match, states/sources, positivity
and global healthy completion remain open. The deformed parent is not
asserted to preserve the frozen example's actual rolling background.
