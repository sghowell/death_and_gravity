# Scalar OS insertion with its finite MS outer reference

Read [FORMULATION.md](FORMULATION.md) and the six proof notes.
The finite outer reference retains the epsilon coefficient of
the fixed inner on-shell slope; its leading value alone is
insufficient.

The report pins eighteen source/proof/test files and the frozen
S6.152 parent. From the repository:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_153
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_scalar_insertion_ms.verify

Native, ordinary, CLI and direct science use unmodified SymPy.
Only full regression uses the audited exact GCD adapter.
Original P8 remains OPEN.
