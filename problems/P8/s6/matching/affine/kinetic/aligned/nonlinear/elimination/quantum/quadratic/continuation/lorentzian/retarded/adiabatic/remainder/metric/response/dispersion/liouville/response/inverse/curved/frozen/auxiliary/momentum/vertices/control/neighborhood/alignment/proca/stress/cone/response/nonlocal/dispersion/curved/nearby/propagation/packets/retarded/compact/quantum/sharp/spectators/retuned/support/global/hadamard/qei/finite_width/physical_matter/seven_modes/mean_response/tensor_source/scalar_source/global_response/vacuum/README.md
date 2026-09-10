# Exceptional vacuum and tree-matching certificate

The [formulation](FORMULATION.md) distinguishes a smooth
exact-tube extension from a real-analytic finite-jet family.
Read the [classical family](notes/family.md),
[analytic extension](notes/analytic.md),
[uniform bounds](notes/bounds.md),
[canonical limit](notes/decoupling.md),
[amplitude and remainder](notes/amplitude.md) and
[scope](notes/scope.md).

The independent native replay, with all repository P8 src
directories on sys.path, is

    python -m p8_exceptional_vacuum.verify --check

The separate complete regression is

    .venv/bin/python -u scripts/p8_snapshot_regression.py

Scientific code, report generation, ordinary replay and CLI
do not monkeypatch SymPy. Only the separately audited complete
regression uses its exact GCD adapter. Source, proof and test
bytes are pinned in the [report](certificates/exceptional-vacuum-tree-matching.json).

The positive heavy tree exchange is a vacuum amplitude match,
not an all-order UV completion or a bounce-parent certificate.
Original P8 remains open.
