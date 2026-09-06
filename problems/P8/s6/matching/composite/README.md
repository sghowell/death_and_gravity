# Composite-metric background gate

The new shared-coupling physical metric admits a regular free-scalar bounce.
A reconstructed canonical potential also realizes CD exactly on
`|T|<=tau/64`, with explicit positive-lapse/kinetic bounds. Global CD,
old-action matching and EFT stability remain open.

See [FORMULATION](FORMULATION.md), [proof](notes/background.md),
[primary-source audit](notes/sources.md), and
[certificate](certificates/composite-background.json).

From the repository root:

```sh
PYTHONPATH=problems/P8/s6/matching/composite/src:problems/P8/s6/matching/bimetric/general/monotonic/src:problems/P8/s6/matching/bimetric/general/src:problems/P8/s6/matching/bimetric/src:problems/P8/s6/matching/src:problems/P8/s6/src:problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/control/src:problems/P8/s5/scattering/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src:problems/P8/s5/matter/physical/control/tree/src:problems/P8/s5/matter/physical/control/tree/loops/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/scalars/src .venv/bin/python -m p8_composite.verify --check
.venv/bin/pytest -q problems/P8/s6/matching/composite/tests
.venv/bin/ruff check problems/P8/s6/matching/composite
```

The verifier is read-only and prints JSON when `--check` is absent.
