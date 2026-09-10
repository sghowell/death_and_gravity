# Paired scalar and gauge fermion vacuum primitives

Direct dimensional vacuum tensors and fixed-scale mass derivatives
agree. The actual scalar vacuum uses its complete physical
mass/residue forest, three beta anchors and a finite remainder.
The gauge vacuum includes all fourteen Dirac flavors, not just
the two active Yukawa flavors.

Run from the repository root:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_149
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_fermion_vacuum.verify

All native, ordinary, CLI and direct science checks use unmodified
SymPy. Only complete regression uses its separately audited adapter.
Exact algebra and written proofs are not formalization or
independent peer review. Original P8 remains OPEN.
