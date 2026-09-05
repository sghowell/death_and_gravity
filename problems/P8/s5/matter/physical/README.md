# CD/M1 physical cubic and quartic reduction

The rolling-matter branch now has fully spatially reduced cubic/quartic phase
kernels and coupled **unnormalized** velocity kernels. All three matter-sourced
momentum constraints, the mixed canonical boundary, both scalars and both tensor
polarizations are retained.

The finite-momentum kinetic matrix differs from the prior principal matrix.
The full two-scalar inverse Hessian is used in quartic Legendre contacts;
diagonal-only or D-only substitution is insufficient.

This is P8-S5.6.CD, not M1 interaction control or full P8 completion.

- [Scope and velocity domains](FORMULATION.md)
- [Reduction and fixed-momentum regularity proof](notes/reduction.md)
- [Independent quadratic derivation and audit](notes/audit.md)
- [Pinned read-only certificate](certificates/cd-interactions.json)

From the repository root:

```sh
PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src .venv/bin/python -m p8_m1_physical.verify --check
.venv/bin/python -m pytest problems/P8/s5/matter/physical/tests -q
.venv/bin/ruff check problems/P8/s5/matter/physical
```

Without `--check`, the verifier prints a candidate report and never writes it.
It pins both reused nonlinear-invariant and metric/York packages and their
previous evidence chain. Replay checks the complete symbolic-time quadratic
bridge, general identities, exact rational-time cubic/quartic examples and
omission controls. Tests additionally check parity/permutation symmetry and
internal TT basis invariance. Symbolic-time quartic expansion is supported but
can be expensive; its expansion is not required by the report.

Next: normalize the coupled finite-momentum modes with time-dependent terms,
control their free propagation, and derive a specified M1 interaction window.
UV matching, loops, inclusive infrared limits and nonlinear stability remain
separate problems.
