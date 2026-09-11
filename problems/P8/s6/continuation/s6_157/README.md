# S6.157: complete matched two-loop Phi b2

The next checkpoint after S6.156 completes the formal
canonical GY14 Phi forward coefficient through two loops.
It keeps the already established unit-disc light pole and
all original P8 scope boundaries.

From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_157
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_full_two_loop_amplitude.verify
    .venv/bin/python -u scripts/p8_replay.py full

Only full regression uses the separately audited exact GCD
adapter. Native, ordinary, CLI and direct science use
unmodified SymPy with interpreter-only recursion and trusted
integer-formatting allowances. Native report generation and
replay are read-only. Do not regenerate frozen report bytes
through another serializer.

See FORMULATION.md and notes/ for the scientific ownership,
regulator and finite-order limitations. Original P8 remains
open; no user intervention is required for the current
research.
