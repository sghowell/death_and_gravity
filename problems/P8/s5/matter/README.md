# P8(b) S5.5: CD/M1 nonlinear preparation and causal robustness

This checkpoint advances the canonical-matter branch independently of
[P8(a)](../../a/README.md). It leaves all earlier certificates unchanged.

The exact nonlinear spatial chart now includes the boundary term required by
the CD model's time-dependent F2. A seven-velocity Legendre transform recovers
the pinned coupled linear coefficients, and a stationary lapse expansion gives
the Hamiltonian through fourth **invariant** order. All 85 recorded lapse and
Hamiltonian coefficients have explicit uniform bounds in local background units.
Spatial momentum constraints, physical cubic/quartic vertices and an interacting
frequency window remain to be established for M1.

Separately, a sharp principal-cone normal form shows that a gravitational diagonal
buffer cannot cure a kinetic/gradient mixing mismatch with frozen canonical
matter. An explicit F deformation preserves the exceptional mixing relation and
gives speeds squared `{1,1/(1+epsilon/d^2)}`; it is healthy and causal for fixed
`epsilon>=0`. Neither this family nor the mismatch theorem computes radiative
corrections or proves technical naturalness.

## Verification

From the repository root:

```sh
PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/matter/src .venv/bin/python -m p8_m1.verify --check
.venv/bin/python -m pytest problems/P8/s5/matter/tests -q
.venv/bin/ruff check problems/P8/s5/matter
```

The replay checks pinned prior sources and reports. Without `--check` it emits
JSON to stdout; it never rewrites a certificate. The earlier full P8 linear
reconstruction remains part of the integration verification, not silently
replaced by loading its output.

- [Frozen scope](FORMULATION.md)
- [Nonlinear and uniform-bound proof](notes/nonlinear.md)
- [Principal-cone theorem, deformation and caveats](notes/robustness.md)
- [Certificate](certificates/cd-matter.json)

Next: include the canonical matter source in the spatial momentum constraints,
derive fully reduced physical vertices with all time-dependent canonical
boundaries, then establish finite-time M1 interaction bounds. The normal form
specifies conditions to test in a later radiative-correction calculation. Neither
P8(b) nor full P8 is complete, and no user choice is currently blocking this work.
