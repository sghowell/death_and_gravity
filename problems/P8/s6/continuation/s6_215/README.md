# P8 S6.215: prepared Ward reconstruction

The conditional Gaussian ADM metric response is reduced to synchronous spatial response and known one-point Ward/chart terms. Different initial-source and final-detector primitives avoid an unaccounted boundary term. Beyond the matched tracefree block, three ordered scalar kernels still require construction and fixed-prescription matching.

See [FORMULATION.md](FORMULATION.md), [notes/scope.md](notes/scope.md), and [notes/validation.md](notes/validation.md). Original P8 remains OPEN.

Read-only repository-root replays:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_215
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_prepared_ward_reconstruction.verify
    .venv/bin/python -u scripts/p8_replay.py full

Native/direct/ordinary/CLI use original SymPy; the exact-GCD adapter is full-regression only. Frozen scientific bytes must not be changed.
