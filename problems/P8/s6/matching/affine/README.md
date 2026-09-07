# S6.37: an exact auxiliary dictionary, not yet a UV parent

An independent affine connection can reproduce the complete original
CD/M1 action without changing its physical metric or free scalar matter.
The full quotient connection Hessian, including unsourced components,
has an explicit uniform inverse bound on the target clock tube.
A regular coefficient ODE matches the lower-order action as well.

The extra degeneracy p^2=1/8 lies outside this tube and is kept as a
negative control. The published lower-order formula is checked directly
from the action, including its missing-X transcription problem.

See [the exact scope](FORMULATION.md), [the dictionary proof](notes/dictionary.md),
[the literal action and rank proof](notes/connection.md),
[the independent audit](notes/audit.md), and [source references](notes/sources.md).

The connection is auxiliary here. Neither its indefinite algebraic
Hessian nor its inverse is a propagating-health or heavy-gap certificate.
The adopted vacuum, finite-gravity and full matching-control obligations
remain separate research. Original P8 remains open.

## Read-only replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/affine/tests", "-q", "-p", "no:faulthandler"]))'
```

With local P8 src directories on the Python path, run
`python -m p8_affine.verify --check`. Without --check the verifier prints
its report and writes no files. The ordinary tests and CLI do not use
the broad-regression exact-GCD adapter. Both random seeds are zero and
the diagnosed host faulthandler timer plugin is disabled.
