# S6.230: complete original flat-vacuum tensor response

The actual retained Minkowski Einstein-plus-Proca TT response has a conjugate pair of growing first-sheet poles. Its complete original causal inverse nevertheless exists on an explicit all-momentum flat TT graph for every finite time window.

This is a mathematical statement about the specified conditional Gaussian mean equation. It does not establish a physical UV no-go or a curved-bounce instability, and it does not close original P8.

## What is proved

- The exact flat three-polarization covariance and original adiabatic variational matching fix the lower polynomial, including the original second vertex and finite m^2 R term.
- With p=lambda^2+P^2, the full physical force operator is p[C+pA2(p)]/(16pi^2 kappa), C=16pi^2 kappa+5m^2/6. The original massive A2 and all fixed finite terms are retained.
- There is exactly one simple zero in each open p half-plane and none on the real axis. Both have negative real p. The complete closing arc is included in the analytic count.
- At mass1000 and kappa10^800, 10^796<|z|<kappa and Re sqrt(z)>10^393. The original physical low-disk propagator ratio has error below10^-796 for |p|<=m^2.
- The full massless/pair/cut dispersion and two finite moment cancellations give an all-spatial-momentum inverse bound4T^2 exp(sqrt(kappa)T), with both graph identities and original causal atoms/germ.
- Real smooth compact-time TT wave-packet forcing can excite the growing modes. Their scale is not a certified physical cutoff or a regime of established full-parent validity.

See [FORMULATION.md](FORMULATION.md) and the seven proof notes in [notes](notes). Exact algebra, scope controls and the source-pinned report are replayable; the continuum proofs are not FORMALIZED.

## Read-only replay from the repository root

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_230
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_full_flat_tensor_response.verify
PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full
```

Native, direct science, ordinary and CLI use original SymPy. Only full regression uses the audited exact-GCD adapter. The report and its17 source files are immutable after their first native build.

The unrestricted curved scalar/clock/matter graph, controlled nonlinear/finite-coupling parent matching and original V/G/B/P8 obligations remain open.
