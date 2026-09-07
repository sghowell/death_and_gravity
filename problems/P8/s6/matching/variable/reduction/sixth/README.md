# S6.32: differentiated coefficients and their cancellation

The first omitted formal action term of the unchanged parent has a
nonzero limiting contribution to a specified off-shell tensor equation,
despite its center action coefficients tending to zero. The previous
formal order contributes an opposite finite term, and their combined
center coefficient is O(c-2). Both results are verified together.

This is neither a full CD propagation equation nor a proof for or
against a controlled all-order reduction. It identifies concrete
effects that a proposed reduction must retain.

See [the statement](FORMULATION.md), [the full proof](notes/tensor.md),
and [the independent-source boundary](notes/sources.md). The report
pins and recursively rebuilds S6.30; no frozen parent is edited.

## Ordinary replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/variable/reduction/sixth/tests", "-q", "-p", "no:faulthandler"]))'
```

With all local P8 src directories on the Python path, the separate CLI
is `python -m p8_sixth_reduction.verify --check`. It regenerates the exact
report in memory and compares it read-only. Without --check it prints
the report and does not overwrite a certificate.

The ordinary checks do not use the broad-regression exact-GCD adapter.
Both random seeds are zero; the diagnosed host faulthandler timer
plugin is disabled. The theorem uses continuous symbolic identities
and a written polynomial inequality, not a fitted parameter grid.

Original P8 remains open. The source/matter dictionary, controlled
full inverse and higher-operator remainder are not discharged by this
single off-shell tensor probe.
