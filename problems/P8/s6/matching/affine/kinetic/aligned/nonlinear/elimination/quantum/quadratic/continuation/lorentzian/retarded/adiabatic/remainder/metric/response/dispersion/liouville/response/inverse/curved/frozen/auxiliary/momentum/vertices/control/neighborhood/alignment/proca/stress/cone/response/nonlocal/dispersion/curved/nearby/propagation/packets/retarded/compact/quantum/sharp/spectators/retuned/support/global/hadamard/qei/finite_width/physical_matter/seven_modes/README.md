# Full seven-mode physical quadratic reference state

See [FORMULATION.md](FORMULATION.md) and the proofs in
[notes/model.md](notes/model.md), [notes/covariance.md](notes/covariance.md),
[notes/preparation.md](notes/preparation.md),
[notes/infrared.md](notes/infrared.md),
[notes/hadamard.md](notes/hadamard.md),
[notes/canonical.md](notes/canonical.md), and [notes/scope.md](notes/scope.md).

The read-only entry point is
`python -m p8_proca_seven_modes.verify --check`, with the repository's
P8 source roots on the module path. The certificate rebuilds S6.102
and all its ancestors, pins every own source/proof/test byte, and
rejects changes to every report field. Scientific SymPy is unchanged.

The deliverable is the actual seven-mode positive global free state,
not a completed stress renormalization, interacting theory or UV verdict.
