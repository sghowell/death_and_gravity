# P8 S6.199: flat Proca full tensor cut

This continuation derives the complete canonical flat-reference Proca
stress cut, including the longitudinal mode, and a controlled
three-subtracted dispersive part. It independently compares the full
tensor ultraviolet weights with the four-dimensional curvature Hessians.

The actual curved CD state and its fixed prescription are unchanged.
This is a matching benchmark, not the outstanding curved contact/endpoint
matching, full parent scattering amplitude, response inverse or P8 closure.

See [FORMULATION.md](FORMULATION.md) and the seven proof notes.
The read-only report is
[polynomial-vacuum-affine-flat-tensor-cut.json](certificates/polynomial-vacuum-affine-flat-tensor-cut.json).

From the repository root:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_199
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_flat_tensor_cut.verify
    .venv/bin/python scripts/p8_replay.py full

Only the full regression uses the separately audited exact GCD adapter.
Original P8(b), V/G/B and P8 remain OPEN.
