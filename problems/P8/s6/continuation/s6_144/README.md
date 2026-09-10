# Complete paired mixed quartic primitive

See [the formulation](FORMULATION.md) and six written proofs in notes.
Run the read-only checkpoint with the repository replay helper:

    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_fermion_mixed_quartic.verify

The report records exact algebra, fixed source hashes, proof gates,
adversarial input controls and the actual conservative enclosure.
Only full regression uses the previously audited exact GCD adapter.
Native, direct science, ordinary and CLI use unmodified SymPy.

This checkpoint bounds the fourth and last primitive quartic row.
It is not full matched two-loop matching or original P8 closure.
