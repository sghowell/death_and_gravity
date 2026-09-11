# S6.170 — Actual curved free Dirac in/out state

This checkpoint extends the free quadratic state calculation to the
actual CD bounce geometry. All geometric fermion copies are retained.
An all-order exact-frame argument identifies the same in/out states as
Hadamard, and twenty frames bound their energy and pressure differences.

Read FORMULATION.md and notes/ before interpreting the report.
All quoted tiny bounds are differences between two states of one
operator. No absolute renormalized curved stress is supplied.

From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_170
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_curved_dirac_state.verify
    .venv/bin/python -u scripts/p8_replay.py full

Only full regression uses the separately audited exact GCD adapter.
The native report, direct science, ordinary tests and CLI use original
SymPy. Every scientific/proof/test/report byte is immutable after freeze.
Original V/G/B and P8 remain OPEN.
