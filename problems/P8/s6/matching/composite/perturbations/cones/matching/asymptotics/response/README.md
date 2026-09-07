# Composite finite-band physical source response

`P8-S6.11.COMPOSITE` adds an explicit finite-positive-parameter result to
the frozen S6.10 limit. A conserved external TT pulse produces an
order-one source-normalized response missed by the specified locked
two-tensor action, even though its individual metric amplitudes can be
made uniformly small. Prepared source-free light data provide the
positive control. General matching and P8 remain open.

Read [FORMULATION.md](FORMULATION.md) and [notes/proof.md](notes/proof.md)
for the exact experiment, source units, very conservative parameter
range, amplitude scaling, and exclusion boundaries. Primary-source
bookkeeping is in [notes/sources.md](notes/sources.md).

Tests from the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/composite/perturbations/cones/matching/asymptotics/response/tests
```

The module `p8_composite_response.verify` prints a replayed report; its
`--check` option compares against
`certificates/composite-response.json`. Set `PYTHONPATH` to this child
and the pinned ancestor `src` directories, as for the prior gates.
The CLI never writes files. Reports and proofs are source-hashed.

The sufficient gates `epsilon<=2^(-187116)` and
`sigma<=2^(-67162)` are factorized exact rationals. They are not numerical
estimates of an optimal background window or a physical cutoff.
