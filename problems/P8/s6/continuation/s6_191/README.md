# C2 covariance-tail response and finite-amplitude remainder

The same actual-state matrix covariance tail now has
first and second shear-amplitude derivative bounds.
A common infinite energy-weighted majorant justifies
C2 differentiation and a Taylor remainder on the stated
compact-time homogeneous C12 family.

This is not the full covariantly subtracted stress
response or gravitational feedback problem. Read the
formulation, mixed.md through tail.md, then scope.md
and validation.md. Original P8 remains OPEN.

From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_191
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_matrix_response_tail.verify
    .venv/bin/python -u scripts/p8_replay.py full

Only final exit code0 and unchanged source/report
hashes establish a pass. No frozen source may change.
