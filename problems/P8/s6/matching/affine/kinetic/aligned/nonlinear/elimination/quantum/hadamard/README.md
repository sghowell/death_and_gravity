# S6.55: controlled Hadamard vector preparation

See [FORMULATION.md](FORMULATION.md), the
[all-order construction](notes/construction.md),
[constrained Proca comparison](notes/hadamard.md) and
[finite observable transfer](notes/transfer.md).

This is a new explicit state with the unchanged action and local
counterterms. It does not relabel the frozen fourth-order state
or establish a full quantum bounce or V/G/B.

The read-only entry point is `p8_vector_hadamard.verify --check`
with all P8 source roots on PYTHONPATH. It checks all local source
hashes and fully rebuilds S6.54 and its frozen ancestry. Ordinary
tests do not use the broad exact GCD adapter.
