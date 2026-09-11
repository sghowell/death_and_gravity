# Current-parent retained homogeneous coupled response

The actual classical constraint and uncancelled conditional vector
response give a unique prepared causal LINEAR response, with an
explicit local physical-lapse recovery bound. This is not a
quantum-corrected background or a small full inverse.

Read FORMULATION.md and all seven written proofs. Verification:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_180
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_coupled_response.verify
    .venv/bin/python scripts/p8_replay.py full

Original P8(b), V/G/B and P8 remain OPEN.
