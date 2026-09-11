# Actual matrix adiabatic reference and covariance tail

This continuation replaces the momentum-growing
covariance comparison by a uniform finite reference
estimate on the specified homogeneous C12 shear class.
It retains the actual all-order state and all three
constrained modes, with no commuting-polarization
assumption or new physical cutoff.

Read FORMULATION.md, then notes/frame.md through
notes/error.md. The scope and validation notes distinguish
this covariance remainder from the remaining renormalized
response and full P8 obligations.

From the repository root, fresh read-only replay is

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_190
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_matrix_adiabatic.verify
    .venv/bin/python -u scripts/p8_replay.py full

Only the final exit status and unchanged source/report
hashes establish a pass. Original P8 remains OPEN.
