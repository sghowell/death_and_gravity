# S6.34: a specified causal response remembers the earlier metric

With zero initial hidden tensor data, a smooth prescribed physical-metric
pulse leaves a quantitatively nonzero hidden response after the pulse
has ended. This is proved for the exact finite-parameter own-f equation,
using whole-box coefficient bounds and a positive-flux Green argument.

The result distinguishes the response from any map of the final metric
germ alone. Its pulses have shrinking physical duration, so it does not
establish a fixed low-frequency EFT exclusion. The prescribed metric is
off shell and may require driving; it is not a free-chi source response.

See [the statement](FORMULATION.md), [the coefficient proof](notes/bounds.md),
[the Green/pulse proof](notes/green.md), and
[the independent-check boundary](notes/sources.md). S6.33 and its frozen
ancestors are recursively rebuilt, not edited.

## Ordinary replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/variable/reduction/exact/retarded/tests", "-q", "-p", "no:faulthandler"]))'
```

With the local P8 src directories on the Python path, the separate CLI is
`python -m p8_own_retarded.verify --check`. It rebuilds the report read-only.
Without --check it prints the report and writes no files. The ordinary
checks use no broad-regression exact-GCD adapter; both seeds are zero and
the diagnosed host timer plugin is disabled.

Original P8 remains open. A declared causal operator is progress beyond
an unselected inverse, but its physical domain and source contract remain
part of the result, not optional caveats.
