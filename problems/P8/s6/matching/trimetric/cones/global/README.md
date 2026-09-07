# S6.15 actual auxiliary-background no-bounce screen

For this constant-link auxiliary parent with positive Einstein terms and
an actual isotropic physical NEC source,
`(H_u/sqrt(K))' = -n_h/(2 K^(3/2)) <= 0`.
It cannot contract and then expand, even at a degenerate bounce. This main
result requires neither subluminal tensor cones nor a flat vacuum and
survives a singular auxiliary tensor denominator.

Read [FORMULATION.md](FORMULATION.md), the [proof](notes/proof.md), and
[source audit](notes/sources.md). Additional, separate results cover
finite CD matching errors, a cone/three-slice alternative and calibrated
mixed-sign tensor denominators. Genuine disconnected and rolling controls
keep their hypotheses distinct.

Exact replay:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/trimetric/cones/global/tests
.venv/bin/ruff check problems/P8/s6/matching/trimetric/cones/global
```

With this package and its pinned ancestor `src` directories on PYTHONPATH,
`python -m p8_trimetric_global.verify --check` performs read-only recursive
replay of `certificates/auxiliary-no-bounce.json`. Without `--check` it
prints a proposed JSON report; it never writes a file. Source hashes cover
only this gate, not future children. The manifest includes independent
exact arithmetic and the separately authored covariant audit.

This excludes the specified parent action as an exact actual CD-bounce
route and gives a conditional controlled-background mismatch threshold.
It does not exclude changed actions or all parents, establish a cutoff,
or change the completed scoped photon objective or frozen linear
classification. Original P8 remains open.
