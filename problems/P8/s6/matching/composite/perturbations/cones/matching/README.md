# Composite light-tensor matching screen

Claim: **P8-S6.9.COMPOSITE**. Status: exact two-TT canonical/retarded algebra
and two specified mass-led hierarchy screens; general matching and P8 open.

The new result keeps the rolling normalization, rotating weighted field,
canonical boundary, and heavy initial data. The symmetric bounce has zero
algebraic relative mass but can have a positive fixed-charge \(k=0\)
frequency. The asymmetric bounce has a subluminal locked speed and a
positive algebraic mass, but changing \(m\tau\) cannot make that mass
dominate both the physical curvature and mixing rate. Neither statement
excludes every finite-band nonadiabatic reduction.

Read [FORMULATION.md](FORMULATION.md), the
[written proof](notes/proof.md), and the [source audit](notes/sources.md).
The read-only report is
[certificates/composite-light.json](certificates/composite-light.json).
It pins and replays the immutable S6.8 certificate and lineage.

## Replay

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/composite/perturbations/cones/matching/tests
.venv/bin/ruff check problems/P8/s6/matching/composite/perturbations/cones/matching
```

For the standalone read-only CLI, put this package and each prior P8
package's `src` directory on `PYTHONPATH`, then run:

```sh
.venv/bin/python -m p8_composite_light.verify --check
```

Without `--check`, the CLI prints the recomputed JSON to stdout. It does
not write a report. The exact tests replay all positive margins,
coefficient identities, omission controls, and inherited hashes.

## Public modules

`canonical` exposes the full physical-clock action, canonical map
coefficients, formal adjoint, inverse-mass representative, exact residual,
and fixed-charge identities. `background.initial_jets(1 or 2)` evaluates
the full reconstruction ODE, including second derivatives. `bounds.build`
certifies the rational hierarchy bounds; `independent.checks` replays
separate coefficientwise Fraction arithmetic. `verify.build_report` and
`verify.validate_report` provide source-hashed read-only replay.
