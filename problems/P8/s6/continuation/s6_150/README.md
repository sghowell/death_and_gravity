# Isolated fixed finite-contact conversion

Read [FORMULATION.md](FORMULATION.md) and the six notes before
using the report. This checkpoint removes one assigned nonlocal
order-two matching contribution by an exact regulated identity,
not by treating a loop-inserted contact as momentum independent.

The report is generated read-only and pins all eighteen source,
proof and test files. Its ancestry is S6.149. From the repository:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_150
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_finite_contact_conversion.verify

Only complete regression uses the audited exact GCD adapter.
Native, ordinary, CLI and direct science use unmodified SymPy.
The isolated affine inverse is not an all-loop physics result.
Original P8 is OPEN.
