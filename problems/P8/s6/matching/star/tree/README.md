# S6.17 finite HR-tree no-bounce screen

A finite tree of constant pairwise HR interactions cannot produce an
actual regular spatially flat contraction-to-expansion transition at a
positive-Einstein physical vertex when all Einstein coefficients are
nonnegative and matter sectors are separately conserved and NEC-respecting.
Zero-Einstein auxiliary vertices and arbitrary algebraic branch changes
are included; cycles, shared sources and a zero-kinetic physical target
without another guarantee are not.

Read [FORMULATION.md](FORMULATION.md), the [proof](notes/proof.md) and
[source audit](notes/sources.md). Actual expanding/mixed-algebraic scalar
solutions and a genuine zero-target exception check the scope.

```sh
.venv/bin/pytest -q problems/P8/s6/matching/star/tree/tests
.venv/bin/ruff check problems/P8/s6/matching/star/tree
```

With this package and pinned ancestor `src` directories on PYTHONPATH,
`python -m p8_hr_tree.verify --check` recursively replays the read-only
`certificates/hr-tree-no-bounce.json`. Without `--check` the verifier emits
a candidate JSON report to stdout; it never writes a certificate. Only
this gate's immediate source/test/document files enter its manifest,
excluding future children.

This is a background and conditional controlled-background matching
screen, not a universal multigravity, perturbative-health, cutoff or UV
verdict. The completed scoped photon objective and frozen linear
classification are unchanged; original P8 remains open.
