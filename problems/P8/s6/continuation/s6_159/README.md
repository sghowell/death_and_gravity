# S6.159: complete second heavy-source reference

This follows S6.158's complete vacuum reference and
regulated first-source-square cancellation. It completes
the named fixed-order Phi/vacuum/source reference assembly,
not the physical truncation or original P8 obligations.

From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_159
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_full_heavy_source.verify
    .venv/bin/python -u scripts/p8_replay.py full

Native, ordinary, CLI and direct science use unmodified
SymPy with interpreter-only recursion and trusted exact
integer-formatting allowances. Only full regression uses
the separately audited exact GCD adapter. Native report
generation and replay are read-only; never reserialize
frozen certificate bytes.

See FORMULATION.md and the six notes. No user-intervention
blocker is asserted; original P8 remains open.
