# Frozen bounce-frequency diagnostic

See [FORMULATION.md](FORMULATION.md) for the definition and
strict scope. The frozen symbol retains the tree operator,
all local Gaussian terms, exact massive bubble and fixed
tadpole box. Its inverse has a positive-real pole in the
specified band; this is not an actual curved-bounce verdict.

The read-only report entry point is
`python -m p8_bounce_frequency.verify --check`, with all
frozen P8 source roots on the import path. The tests provide
that path without installation.

No earlier checkpoint, state, finite prescription or physical
action is overwritten.
