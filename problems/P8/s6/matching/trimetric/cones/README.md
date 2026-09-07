# S6.14 actual auxiliary-metric tensor cones

The positive-link auxiliary background equations force at least one full
tensor cone outside the actual matter cone when canonical matter rolls.
On the exchange-symmetric branch, the physical common tensor itself has
that excess; its exact auxiliary response does not hide it in a heavy mode.

Read [FORMULATION.md](FORMULATION.md), the [proof](notes/proof.md), and
[source/convention audit](notes/sources.md). The proof includes an actual
local free-scalar rolling solution and distinct mixed-sign controls.

Exact replay:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/trimetric/cones/tests
.venv/bin/ruff check problems/P8/s6/matching/trimetric/cones
```

With this package and all pinned ancestor `src` directories on PYTHONPATH,
`python -m p8_trimetric_cones.verify --check` performs read-only recursive
replay of `certificates/actual-tensor-cones.json`. Without `--check` it emits
a proposed JSON report to stdout; it never writes a certificate. Source
hashes deliberately do not recurse into future children. Independent exact
Fraction/polynomial checks and the separately authored covariant audit are
part of the manifest.

The asymmetric statement concerns full-parent principal cones. The stronger
common-channel statement requires exact exchange symmetry. Neither is a
universal obstruction to all light-only EFTs or all auxiliary-parent actions.
The result is a scoped matching screen; original P8 remains open.
