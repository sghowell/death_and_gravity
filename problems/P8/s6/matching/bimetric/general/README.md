# Arbitrary-beta regular-flat HR screen

All constant HR interaction coefficients are covered, including Bianchi
polynomial roots and branch-switching points. Positive Einstein terms and
separate classical NEC matter rule out a regular nondegenerate physical
bounce in the common flat FLRW sector. This rejects a parent class, not
general bimetric models or UV completions, and does not close P8.

See [FORMULATION](FORMULATION.md), [proof](notes/no-bounce.md),
[source audit](notes/sources.md), and
[certificate](certificates/regular-flat-hr.json).

From the repository root:

```sh
PYTHONPATH=problems/P8/s6/matching/bimetric/general/src:problems/P8/s6/matching/bimetric/src:problems/P8/s6/matching/src:problems/P8/s6/src:problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/control/src:problems/P8/s5/scattering/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src:problems/P8/s5/matter/physical/control/tree/src:problems/P8/s5/matter/physical/control/tree/loops/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/scalars/src .venv/bin/python -m p8_bimetric_general.verify --check
.venv/bin/pytest -q problems/P8/s6/matching/bimetric/general/tests
.venv/bin/ruff check problems/P8/s6/matching/bimetric/general
```

The verifier is read-only. Without `--check` it prints JSON; it never writes
or replaces certificates.
