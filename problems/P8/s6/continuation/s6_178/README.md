# Finite-time complete affine-source response

This continuation controls the retarded sourced mean of the same
conditional Gaussian vector on the fixed CD metric. Its full-source
Sobolev and cubic light-force estimates retain the temporal constraint,
source contact, homogeneous cancellation and nonclosed source cases.

Read FORMULATION.md and the seven written proofs. Verification:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_178
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_retarded_energy.verify
    .venv/bin/python scripts/p8_replay.py full

Original P8(b), V/G/B and P8 remain OPEN.
