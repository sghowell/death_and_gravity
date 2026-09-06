# Composite large-asymmetry response

Claim: **P8-S6.10.COMPOSITE**. The exact scaled limit has a growing
fixed-charge homogeneous relative tensor, and its composite-metric
projection does not disappear. This advances the earlier centre-only
screen to an actual time-dependent limiting equation.

The result is deliberately not a general EFT exclusion. Finite-parameter
convergence is proved only on each fixed compact forward interval; the
limiting individual \(g\) metric is degenerate, and no finite-\(k\),
cutoff, nonlinear or full-duration conclusion is made.

Read [FORMULATION.md](FORMULATION.md), [the proof](notes/proof.md), and
[the source audit](notes/sources.md).

## Implementation

`scaled` contains the full parameter-dependent CD ODE and regularized
canonical coefficient templates. `limiting` contains the exact
integrable flow, full time-dependent relative coefficient, canonical
limits and physical projection checks. `independent` provides separate
coefficientwise Fraction arithmetic.

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/composite/perturbations/cones/matching/asymptotics/tests
.venv/bin/ruff check problems/P8/s6/matching/composite/perturbations/cones/matching/asymptotics
```

The read-only verifier is `p8_composite_asymptotics.verify`, with the
source-hashed report at
[certificates/composite-asymptotics.json](certificates/composite-asymptotics.json).
With this package and its pinned ancestors' `src` directories on
`PYTHONPATH`, run:

```sh
.venv/bin/python -m p8_composite_asymptotics.verify --check
```

Without `--check` the verifier prints recomputed JSON to stdout and does
not write a report.
