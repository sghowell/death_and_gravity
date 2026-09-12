# P8 S6.214: full constrained ADM Gaussian vertices

Complete lapse, shift and spatial-trace variations now have explicit constrained first/second vertices and all-momentum energy-relative form bounds. The full temporal constraint, longitudinal shift and nonlinear four-metric chart contacts remain. This is not full new-sector continuum matching or a reduced inverse; original P8 remains OPEN.

See [FORMULATION.md](FORMULATION.md), [notes/scope.md](notes/scope.md) and [notes/validation.md](notes/validation.md).

Read-only repository-root replays:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_214
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_full_adm_vertices.verify
    .venv/bin/python -u scripts/p8_replay.py full

Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. Never edit frozen evidence to force a passing report.
