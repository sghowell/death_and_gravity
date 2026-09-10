# P8 S6.143: finite outer MS conversion

This checkpoint supplies S6.138's missing
finite outer interaction-reference conversion.
The exact dimensional spectral integral,
Laurent finite part and nonzero mass-ratio
correction are all retained.

See [formulation](FORMULATION.md),
[dimensional derivation](notes/dimensional.md),
[finite part](notes/finite_part.md),
[remainder](notes/remainder.md),
[parent conversion](notes/conversion.md)
and [scope](notes/scope.md).

The certificate is
certificates/polynomial-vacuum-fermion-outer-ms.json.
Use the supported ordinary, CLI and full
repository replay commands. Only full
regression uses the audited exact GCD
adapter; native, direct science, ordinary
and CLI use unmodified SymPy.

Other matching work and the complete two-loop
error remain open. Original P8 is not closed.
