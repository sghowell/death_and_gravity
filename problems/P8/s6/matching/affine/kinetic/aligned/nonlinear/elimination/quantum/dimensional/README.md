# S6.52: dimensional local vector matching controls

See [FORMULATION.md](FORMULATION.md) and the
[dimensional, counterterm and Green-pole proof](notes/proof.md).

This checkpoint distinguishes the actual clock-sensitive lapse poles,
the matched finite ordinary-Proca diagnostic, and the still-unmatched
finite clock-mass curvature terms. It does not close original P8.

The read-only entry point is `p8_vector_dimensional.verify --check`
with the P8 source roots on PYTHONPATH. It checks every local source
hash and rebuilds S6.51 and its frozen ancestry. Ordinary tests do not
use the broad exact GCD adapter.
