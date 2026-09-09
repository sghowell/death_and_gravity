# Actual nearby quadratic scalar CCR

This checkpoint quantizes the specified reduced quadratic scalar
system and proves its state-independent compact-probe response.
It does not quantize the interacting parent.

Read FORMULATION.md and the proofs in notes/canonical.md,
notes/growth.md, notes/state.md, notes/response.md and notes/scope.md.
The canonical density variables keep the volume and action
normalizations explicit. Polynomial Fourier estimates make the
positive Gaussian covariance a state on actual evolved linear
fields, including spacetime smearings.

The read-only verifier rebuilds the pinned quantitative compact
classical-response parent. Scientific checks use native exact
arithmetic; ordinary and complete regression runs are separate.
The report is generated only after source freezing and then
replayed independently. The CLI supports --check and never writes.

The compact commutator is nonzero outside the physical matter cone
of this exact quadratic theory. That is not a cutoff-valid UV
exclusion: omitted operators, interacting corrections, temporal
frequencies and the parent's own high-frequency response are not
bounded by this checkpoint.
