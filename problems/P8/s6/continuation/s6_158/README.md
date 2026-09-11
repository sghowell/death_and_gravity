# S6.158: complete two-loop vacuum reference

This follows the complete matched Phi pole and forward
coefficient of S6.157. It advances the local vacuum reference
and first-source-square ownership, not full source or P8
closure.

Read FORMULATION.md and the six notes. From the repository:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_158
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_full_vacuum_reference.verify
    .venv/bin/python -u scripts/p8_replay.py full

Only full regression uses the separately audited exact
GCD adapter. Native, ordinary, CLI and direct science use
unmodified SymPy with interpreter-only recursion and trusted
integer-formatting allowances. Native generation and replay
are read-only. Never replace a frozen report by a different
serialization.

Original P8 remains open; no user-intervention blocker is
asserted by this checkpoint.
