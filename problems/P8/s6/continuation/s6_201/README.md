# P8 S6.201: curved reference initial-state prefactor

The complete initial-alpha prefactor can be removed from the W8
comparison response with a controlled finite correction. All retarded
memory pairs, the full local contact, occupation products and distinct
regulator tails are retained. The actual all-order CD state is unchanged.

The resulting unit-W8 readout is Wronskian-normalized but is not an
exact bisolution or a new physical state. Its joint spatial estimates
and fixed covariant contact/endpoint matching remain open.

See [FORMULATION.md](FORMULATION.md), the seven proof notes and the
[read-only report](certificates/polynomial-vacuum-affine-reference-state-prefactor.json).

From the repository root:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_201
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_reference_state_prefactor.verify
    .venv/bin/python scripts/p8_replay.py full

Native, direct science, ordinary and CLI use original SymPy.
Only full regression uses the independently audited exact GCD adapter.
Original P8(b), V/G/B and P8 remain OPEN.
