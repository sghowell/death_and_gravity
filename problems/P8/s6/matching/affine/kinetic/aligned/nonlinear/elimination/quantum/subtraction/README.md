# S6.51: an explicit finite adiabatic vector subtraction

See [FORMULATION.md](FORMULATION.md) and the
[exact subtraction and integrated-error proof](notes/proof.md).

The result bounds the stated subtracted Gaussian energy and pressure
on the original rolling interval. It does not identify a full finite
covariant matching prescription, all-order Hadamard state or V/G/B
completion. Original P8 remains open.

The read-only entry point is `p8_vector_subtraction.verify --check`
with the P8 source roots on PYTHONPATH. Replay checks all local source
hashes and rebuilds the frozen ancestry. Ordinary tests do not use
the broad exact GCD adapter.
