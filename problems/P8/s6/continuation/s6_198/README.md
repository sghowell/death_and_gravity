# S6.198: retarded reference boundaries and finite bulk

Six exact integrations isolate a convergent reference bulk
and endpoint, with all spatial momenta retained. The full
metric contact and five earlier equal-time terms remain
for covariant spatial matching. Original P8 remains OPEN.

Read FORMULATION.md and the seven notes.

## Read-only replay

From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_198
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_retarded_reference_boundary.verify
    .venv/bin/python -u scripts/p8_replay.py full

Eighteen scientific sources and the native report are
hash-pinned. Native, direct science, ordinary and CLI
use original SymPy; only full regression uses the audited
exact GCD adapter. Analytic arguments are not FORMALIZED.
