# Complete double-bubble MS interaction forest

Read [FORMULATION.md](FORMULATION.md) and the six proof notes.
The report pins eighteen source/proof/test files and the frozen
S6.151 parent. From the repository:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_152
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_double_bubble_ms.verify

Native, ordinary, CLI and direct science use unmodified SymPy.
Only full regression uses the audited exact GCD adapter.
The assigned graph family is not the complete matched amplitude.
Original P8 is OPEN.
