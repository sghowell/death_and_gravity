# Complete two-loop fermion-sector ownership, not a total error bound

P8 S6.137 derives the mixed two-loop functional after integrating
the fermions and removing reducible source terms by the Legendre
transform. Its noncommuting background expansion gives six scalar
and three gauge primitive rows through four external Phi fields.

The report distinguishes one bounded paired quadratic row,
one partially bounded quartic row, seven unevaluated primitive
rows, and the separate counterterm/finite-conversion tasks.
It does not mark original P8 complete.

Read FORMULATION.md and all six notes. Read-only native replay:

    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_two_loop_fermion_ledger.verify

The immediate ancestor is frozen S6.136. Source hashes and
every report field are compared by read-only replay. Native
and ordinary scientific computation uses unmodified SymPy.
