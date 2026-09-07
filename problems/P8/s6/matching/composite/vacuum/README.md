# Composite physical vacuum-branch matching audit

`P8-S6.12.COMPOSITE` uses the physical composite metric as a field
coordinate. Exchange symmetry makes the zero-relative branch exact and
source preserving: its classical light action is Einstein gravity plus
the unchanged input matter. This directly obstructs a small-error match
of that branch to the original CD/M1 action and physical bounce.

Read [FORMULATION.md](FORMULATION.md) and [notes/proof.md](notes/proof.md)
for the distinct matter classes, literal physical source normalization,
conditional functional-uniqueness theorem, vacuum shifts and error
budgets. [notes/sources.md](notes/sources.md) records the primary-source
boundary. No loop, general parent, cutoff or P8-closure verdict follows.

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/composite/vacuum/tests
```

With this child and the pinned ancestor `src` directories on `PYTHONPATH`,
`python -m p8_composite_vacuum.verify` prints the source-hashed report;
`--check` compares it with `certificates/composite-vacuum.json`.
The verifier performs no writes and pins the old CD/M1 classification
as well as S6.11. Ancestor certificates remain immutable.
