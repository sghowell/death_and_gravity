# S6.212: dimensional spatial ultraviolet finite part

The actual dimension-dependent Proca modes and full field-strength vertex determine the original spatial UV finite difference after the fixed pole subtraction. Evanescent odd endpoints, fixed-source tensor reconstruction, curvature and volume variations are retained.

This is a local finite UV coefficient, not the completed dimension-limit interchange, assembled spatial response, reduced inverse, background or original P8 closure.

See [formulation](FORMULATION.md), [geometry](notes/geometry.md), [modes](notes/modes.md), [dimension jets](notes/dimension.md), [curvature](notes/curvature.md), [finite part](notes/finite.md), [scope](notes/scope.md), and [verification](notes/validation.md).

Read-only frozen replay commands:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_212
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_dimensional_spatial_symbol.verify
    .venv/bin/python -u scripts/p8_replay.py full

Only full regression uses the audited exact GCD adapter. Native, direct science, ordinary and CLI retain original SymPy.
