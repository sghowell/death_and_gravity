# Complete GY14 two-loop Phi normalization and unit-disc pole

Read [FORMULATION.md](FORMULATION.md) and the six proof notes.
The full dimensional coordinate shift matters: its epsilon
coefficient cancels a field-pole product and induces a finite
second-order heavy-mass parameter shift.

The report pins eighteen source/proof/test files and S6.155.
From the repository:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_156
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_full_phi_normalization.verify

Native, ordinary, CLI and direct science use unmodified SymPy.
Only full regression uses the audited exact GCD adapter.
Original P8 remains OPEN.
