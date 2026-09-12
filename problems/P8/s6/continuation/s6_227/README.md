# S6.227 — full original finite local reference inverse

Original P8 remains OPEN. This checkpoint retains the COMPLETE original finite local Hessian inside a controlled curvature-coordinate reference inverse. The exact factorization and source-time adjoint estimate prove an all-momentum weighted bound for the previously unbounded-looking local remainder.

The actual nonlocal curved mass/state/contact/classical/matter correction and full S222 coupled graph remain open. The time weight is not physical damping; unweighting can be extremely costly, and neither stability nor a controlled bounce is claimed.

See [FORMULATION.md](FORMULATION.md) and the seven proof notes. Replay:

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_227
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_full_finite_local_reference.verify
PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full
```

All scientific/proof/test/report bytes are frozen after preflight. Fresh outcomes are recorded externally. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter.
