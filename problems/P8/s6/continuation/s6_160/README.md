# S6.160: full two-loop derivative-action physical-source transport

Parent: S6.159, complete named fixed-order Phi pole, b2, vacuum and
stationary-H reference assembly. All prior scientific files and
reports remain unchanged.

This checkpoint uses the full generated action and the matched
physical source J F(Psi). It transfers the canonical parent
correlators through two loops without assuming a bounded
ordinary-Psi second MS residue. See FORMULATION.md and all six
notes for the precise prescription and limits.

The native report pins 18 own source, test and proof files.
Run from the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_160
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_two_loop_physical_source_map.verify

Ordinary, CLI, direct science and native generation use unmodified
SymPy. The full-regression-only exact arithmetic adapter is not
used to generate the native report. Written analytic arguments
and replay are not formalization or independent peer review.

Original P8 is OPEN. Ordinary-Psi composite mixing, physical
truncation, V contours/cuts, G and B are not closed here.
