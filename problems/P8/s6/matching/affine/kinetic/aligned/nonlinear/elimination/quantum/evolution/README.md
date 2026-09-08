# S6.54: differentiated vector integral and clock source

See [FORMULATION.md](FORMULATION.md), the
[uniform derivative proof](notes/derivative.md) and
[clock-source and scale proof](notes/clock-source.md).

The unchanged Gaussian preparation and fixed matching prescription
give bounded C1 energy/pressure and a homogeneous clock source.
This is not all-order state admissibility or a corrected bounce.

The read-only entry point is `p8_vector_evolution.verify --check`
with all P8 source roots on PYTHONPATH. It checks all local source
hashes and fully rebuilds S6.53 and its frozen ancestry. Ordinary
tests do not use the broad exact GCD adapter.
