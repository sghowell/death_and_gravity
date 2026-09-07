# P8 S6.20: variable-beta local bounce and exact inner tensor audit

This isolated gate supplies an actual local canonical-clock/free-`chi`
solution with the CD scale factor, for a specified scalar-dependent HR
interaction. It retains the sourced clock in the undivided Bianchi equation,
and computes the literal two-TT and transverse-shift action.

The prospective large mass at `c=2+epsilon^2` has a shrinking variation layer:
the full canonical inner equation is
`Q_xx+80 Q/(1+8x^2)=0`, not an increasingly adiabatic algebraic constraint.
The associated adiabatic ratios remain nonzero. This identifies a required
nonadiabatic calculation, not a light-only matching exclusion.

Read [FORMULATION.md](FORMULATION.md), the [proof](notes/proof.md), and the
[source audit](notes/sources.md). The report is
`certificates/variable-beta-local.json`. It pins and replays S6.17 and the
adopted S6 matching contract; it does not depend on S6.19.

From the repository root, the ordinary tests discover the local source paths:

```sh
.venv/bin/pytest problems/P8/s6/matching/variable/tests -q
.venv/bin/ruff check problems/P8/s6/matching/variable
```

The CLI is read-only. With the new package and its frozen ancestor packages
on `PYTHONPATH`, run `python -m p8_variable_beta.verify --check`. Without
`--check` it emits the candidate report on stdout; it never writes artifacts.
Source hashes cover immediate modules/tests/docs only, allowing future
children without changing this gate.

No complete-CD solution, full stability, physical retarded matching error,
cutoff, loop/UV completion, or original C/D operator map is supplied. This
does not change the completed scoped photon objective or frozen linear
classification. Original P8 remains open.
