# P8 S6.213: complete matched tracefree spatial Gaussian current

The same prepared Gaussian sector now has a complete original-prescription tracefree spatial weak current and a uniform original-regulator error, anchored to the exact homogeneous current. The full dimensional limit is justified on an explicit disk, all contact/finite terms remain, and the derivative loss is stated explicitly. Original P8 remains OPEN.

See [FORMULATION.md](FORMULATION.md) for the claim, [notes/scope.md](notes/scope.md) for the remaining obligations and historical metadata erratum, and [notes/validation.md](notes/validation.md) for the verification boundary.

Read-only replay from the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_213
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_matched_spatial_current.verify
    .venv/bin/python -u scripts/p8_replay.py full

Native, direct science, ordinary and CLI retain original SymPy. Only full regression uses the audited exact-GCD adapter. Do not edit frozen evidence to regenerate a passing report.
