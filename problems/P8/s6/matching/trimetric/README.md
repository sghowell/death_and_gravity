# Auxiliary trimetric parent prerequisite audit

`P8-S6.13.AUXILIARY` tests the full field equations of one specified
auxiliary parent. With no separate g/f cosmological terms, either
nonzero link obstructs every regular constant flat triple, irrespective
of a constant matter vacuum energy. A separately named beta4 extension
restores a flat vacuum and has positive relative quadratic mass.

Read [FORMULATION.md](FORMULATION.md), [notes/proof.md](notes/proof.md)
and [notes/sources.md](notes/sources.md). The physical matter metric,
source-induced Lorentz constraint, and restricted effective potential
are retained; this is not imported as the old beta2 composite parent.

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/trimetric/tests
```

With this and the pinned ancestor source directories on PYTHONPATH,
`python -m p8_trimetric.verify` prints the report; `--check` compares it
with `certificates/auxiliary-parent.json`. The verifier performs no
writes. No ancestor certificate, P8 row, or matching contract is changed.
