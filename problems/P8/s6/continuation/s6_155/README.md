# First finite Phi normalization with its MS projection commutator

Read [FORMULATION.md](FORMULATION.md) and the six proof notes.
A finite field factor may multiply poles before taking the
regulator limit. The same factor must be used on vertices,
propagators and assigned counterterms.

The report pins eighteen source/proof/test files and S6.154.
From the repository:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_155
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_finite_field_covariance.verify

Native, ordinary, CLI and direct science use unmodified SymPy.
Only full regression uses the audited exact GCD adapter.
Original P8 remains OPEN.
