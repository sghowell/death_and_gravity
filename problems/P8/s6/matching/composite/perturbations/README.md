# S6.7: the frozen composite candidates fail a local vector-principal screen

The exact tensor/vector calculation finds positive, regular vector inertia
but a negative formal physical vector-gradient coefficient in a punctured
neighborhood of both frozen S6.6 bounces. The vector shift constraint is
regular at the zero-stiffness bounce itself. Changing the independent
interaction-to-duration scale can reverse the leading sign, so this is not
a general composite-model exclusion.

- [Scope and claim](FORMULATION.md)
- [Written derivation and analytic sign proof](notes/modes.md)
- [Primary-source and normalization audit](notes/sources.md)
- [Replay certificate](certificates/composite-modes.json)

Run the complete child tests, including the read-only lineage replay:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/composite/perturbations/tests
.venv/bin/ruff check problems/P8/s6/matching/composite/perturbations
```

The replay entry point is `python -m p8_composite_modes.verify --check`,
with this `src` and the pinned S6/S5 ancestor package `src` directories on
PYTHONPATH, as enumerated in `tests/conftest.py`. Without `--check` it prints
the deterministic JSON report to stdout; neither mode writes files.

Public APIs are `tensor.derive/checks`, `shift.derive/checks`,
`vector.derive/checks/canonical_frequency`, `jets.derive/checks`, the separate
Fraction engine `independent.checks`, and `verify.build_report`.

Only the formal local vector-principal screen is closed for the frozen
examples. Scalar gradients, physical rolling gaps, nonlinear interactions,
cutoff hierarchy, controlled matching and P8 completion remain open. Mu=0
does not establish that a canonically normalized physical gap vanishes.
