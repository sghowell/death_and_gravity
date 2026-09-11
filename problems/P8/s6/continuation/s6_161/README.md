# S6.161: complete through-two-loop low-energy elastic cut

The unchanged parent S6.160 transports the named fixed-order
canonical references with the full derivative action and
matched physical source. This checkpoint bounds the complete
one-loop scattering amplitude on the physical [4,6] window
and hence its second elastic absorptive coefficient.

Read FORMULATION.md and all six notes. The native report pins
18 own source/proof/test files. From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_161
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_two_loop_elastic_cut.verify

Native generation, direct science, ordinary and CLI replay
use unmodified SymPy. The separately audited exact-arithmetic
adapter remains restricted to full regression. Written
analysis and independent tests are not formalization or
independent peer review.

The formal cut-subtracted margin is not strict physical
positivity without the omitted-order and dispersion controls.
Original P8 is OPEN.
