# Scalar quadratic finite MS slope and converted OS remainder

Read [FORMULATION.md](FORMULATION.md) and the six proof notes.
A finite MS slope must be recorded before the outer on-shell
operation removes it. The existing 32-refinement forest and
first-sheet integral estimates are retained.

The report pins eighteen source/proof/test files and S6.153.
From the repository:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_154
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_scalar_ms_slopes.verify

Native, ordinary, CLI and direct science use unmodified SymPy.
Only full regression uses the audited exact GCD adapter.
Original P8 remains OPEN.
