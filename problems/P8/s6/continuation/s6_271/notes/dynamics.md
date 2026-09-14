# Entire unitary, changed-readout and two-cutoff comparisons

The complete real Weyl symbol g_ext,c is the actual time-dependent scalar
g0 plus a smooth compact phase symbol. Its Weyl operator A_W,c is bounded
and self-adjoint by the explicit kernel estimate. Its non-scalar kernel is
Schwartz. The inherited smooth real-time dependence and uniform phase
derivative estimates imply norm continuity on |u| <= T = 10^-2000.

The bounded time-dependent Dyson construction therefore gives an exact
unitary U_W,c with U_W,c(0) = I. This is not a finite occupation-space
truncation or finite Taylor series. The scalar g0 contributes its full
phase. The smoothing-kernel seminorm argument used for the parent also
preserves Schwartz space: each seminorm of the smoothing term is bounded
by a constant times a Hilbert norm, and the scalar part acts exactly.
The common Schwartz strong equation follows with the inherited time
regularity. No assertion of an unregularized singular-Hamiltonian domain
is made.

## Unitary operator-norm ordering comparison

Let U_C,c be the parent's exact calibrated-coherent evolution, with
||A_C,c|| < B = 10^1010. Unitarity on both sides of the full Duhamel
identity gives, for either sign of u,

    ||U_W,c(u)-U_C,c(u)||
      <= integral_0^|u| ||A_W,c(s)-A_C,c(s)|| ds
      <= |u| epsilon_H < 10^-944.

Here epsilon_H is the actual rational ordering error from
notes/comparison.md, bounded strictly below 10^1056. There is no unnecessary
exponential amplification factor: both propagators in Duhamel are unitary.
This is an OPERATOR norm statement, not only a Gaussian-seed estimate.

Also ||A_W,c|| <= B + epsilon_H, and

    ||U_W,c(u)-I|| <= |u|(B+epsilon_H) < 10^-943.

Neither estimate drops either Hamiltonian's scalar phase. The original
reference picture is recovered by U_full,W,c = U_ref U_W,c, with the
same U_ref as before. Multiplication by this shared unitary preserves the
comparison norms.

## Change both state and readout

For the same normalized original seed psi0, let psi_W,c = U_W,c psi0 and
psi_C,c = U_C,c psi0. The difference of volume expectations satisfies

    |<psi_W,c,F_W,c psi_W,c> - <psi_C,c,F_C,c psi_C,c>|
      <= ||F_W,c-F_C,c|| + 2||F_C,c|| ||psi_W,c-psi_C,c||
      <= epsilon_F + 4 |u| epsilon_H < 2 * 10^-200,

because ||F_C,c|| < 2. Thus both the readout operator and the evolved state
are changed; a state-only comparison would be insufficient.

The physical readout is U_ref F_(W or C),c U_ref^*, paired with the physical
state U_ref psi_(W or C),c. Conjugating both gives exactly the same
expectation and the same bounds. Fixed-window coherent covariance is not
silently used instead of this explicit conjugation.

## Compare the two defined Weyl cutoffs

The parent already bounds the calibrated-coherent c = 1 versus c = 2 state
difference by less than 10^-1970 and its volume mean difference by less
than 10^-1230. Applying the triangle inequality through those comparisons,

    ||psi_W,1-psi_W,2|| <= 2 |u| epsilon_H + old_state_error < 10^-943,

    |mean_W,1-mean_W,2|
      <= 2(epsilon_F+4|u|epsilon_H) + old_mean_error < 10^-199.

The exact rational sums are checked. The old, much stronger coherent-only
numbers are NOT relabeled as Weyl bounds. No arbitrary other cutoff is
included in this two-cutoff statement.

## Leakage and original constraints

Let E_out be the original positive coherent-POVM effect for labels outside
the core. It satisfies 0 <= E_out <= I, but is not a phase-space support
projection. The original seed probability is below the parent's exact gamma
tail bound. Therefore

    <psi_W,c,E_out psi_W,c>
      <= original_tail + 2||psi_W,c-psi0|| < 10^-942.

This is small but not zero. It does not define a new conditioned seed.

The actual finite translation action is symplectic and preserves the
original covariance, so in the whitened coordinates it is also orthogonal.
It preserves both radial cutoffs and the entire reconstructed symbols.
Weyl metaplectic covariance then makes both Weyl operators commute with
the full unitary translation group. The original common zero-charge
sector is reducing for the full evolution. No generated spatial means or
residual charges are deleted, and a mean-zero descriptor subspace is not
mistaken for a closed Lie algebra.

Independent noncommuting finite-dimensional fixtures check the unitary
Duhamel inequality, readout change and physical-picture conjugation.
They test the comparison mechanism, not the actual full P8 quantum motion.
