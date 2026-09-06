# Exact beta1 parent screen

The specified positive-Einstein, canonical-matter bimetric family has a
healthy quadratic Minkowski spectrum and a source-preserving flat quadratic
tree-level Weyl coefficient. It nevertheless fails the exact regular CD
bounce equations. General S6 matching and P8 remain open.

See [FORMULATION](FORMULATION.md), [proof](notes/parent-screen.md),
[primary-source audit](notes/sources.md), and
[certificate](certificates/beta1-parent-screen.json).

From the repository root:

```sh
PYTHONPATH=problems/P8/s6/matching/bimetric/src:problems/P8/s6/matching/src:problems/P8/s6/src:problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/control/src:problems/P8/s5/scattering/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src:problems/P8/s5/matter/physical/control/tree/src:problems/P8/s5/matter/physical/control/tree/loops/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/scalars/src .venv/bin/python -m p8_bimetric.verify --check
.venv/bin/pytest -q problems/P8/s6/matching/bimetric/tests
.venv/bin/ruff check problems/P8/s6/matching/bimetric
```

The verifier is read-only. Without `--check` it prints a freshly replayed
JSON report; it does not overwrite the checked-in certificate.
