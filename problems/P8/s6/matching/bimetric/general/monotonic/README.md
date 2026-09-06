# Normalized-Hubble monotonicity

The regular common-flat HR screen now covers degenerate as well as
nondegenerate contraction-to-expansion transitions. Its sharp CD endpoint
bound needs no derivative-error hypothesis. General parent matching is
still open.

See [FORMULATION](FORMULATION.md), [proof](notes/monotonicity.md),
[source audit](notes/sources.md), and
[certificate](certificates/regular-hr-monotonic.json).

From the repository root:

```sh
PYTHONPATH=problems/P8/s6/matching/bimetric/general/monotonic/src:problems/P8/s6/matching/bimetric/general/src:problems/P8/s6/matching/bimetric/src:problems/P8/s6/matching/src:problems/P8/s6/src:problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/control/src:problems/P8/s5/scattering/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src:problems/P8/s5/matter/physical/control/tree/src:problems/P8/s5/matter/physical/control/tree/loops/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/scalars/src .venv/bin/python -m p8_bimetric_monotonic.verify --check
.venv/bin/pytest -q problems/P8/s6/matching/bimetric/general/monotonic/tests
.venv/bin/ruff check problems/P8/s6/matching/bimetric/general/monotonic
```

The verifier is read-only and prints JSON when `--check` is absent.
