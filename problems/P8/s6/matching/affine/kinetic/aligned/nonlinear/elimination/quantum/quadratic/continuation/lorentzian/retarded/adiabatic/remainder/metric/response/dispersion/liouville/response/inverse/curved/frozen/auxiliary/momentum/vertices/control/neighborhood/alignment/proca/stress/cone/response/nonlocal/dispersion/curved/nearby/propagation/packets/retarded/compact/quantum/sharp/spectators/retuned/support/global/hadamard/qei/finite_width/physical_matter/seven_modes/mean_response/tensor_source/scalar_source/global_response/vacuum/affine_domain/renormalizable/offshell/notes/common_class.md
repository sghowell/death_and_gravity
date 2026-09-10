# One nonempty function class controls both tree errors

Work in the canonical flat coordinates and take real Schwartz
Psi whose Fourier transform is supported in the Euclidean
four-momentum unit ball. Require

    integral |Psi_tilde(p)| d^4p/(2pi)^4 <=1,
    U=||Psi||_L2(spacetime).

The class is nonempty: scale any nonzero real-even smooth
compactly supported Fourier bump inside that ball.
The Fourier L1 condition bounds every derivative of Psi
in sup norm by one; Parseval bounds each derivative in L2
by U. This is a chosen estimate in fixed coordinates, not
a new Lorentz-invariance assumption about the theory.

R is cubic in derivatives of Psi. Derivatives do not enlarge
Fourier support and products add their supports, so

    support radius Psi <=1,
    support radius Phi=Psi+R(Psi) <=3,
    support radius J=Phi^2 <=6.

On the last support, |p^2_Minkowski|<=|p|_Euclidean^2<=36.
Thus K=Box+2 has multiplier |2-p^2_Minkowski|<=38<D.
The exact heavy stationary solution exists on this window
and the heavy inverse is a bounded real multiplier.

A cubic monomial in R has L2 norm at most U: put one
factor in L2 and the other two in L-infinity. Consequently

    ||R||_infinity<=C_R, ||R||_2<=C_R U,
    ||Phi||_infinity<=1+C_R, ||Phi||_2<=(1+C_R)U,
    ||Phi^2||_2^2<=(1+C_R)^4 U^2.

Every monomial in the local substitution remainder has
at least six factors, all field derivatives of order at
most ten. Put two factors in L2 and the remaining factors
in L-infinity. Its spacetime absolute integral is at most
U^2 times its coefficient. The pointwise coefficient
majorant E_field(1) therefore also gives an integrated
bound E_field(1) U^2 on this particular class.

The boundary currents are Schwartz and have integrable
derivatives. Their integrals vanish. Finally,

    |S_polynomial[Phi,H_*]-S_target[Psi]|
      <= { E_field(1)
           +G^2*38^4/[8D^5(1-38/D)]*(1+C_R)^4 } U^2
      <10^-800 U^2   for U>0.

No pointwise bound was integrated against an infinite
volume factor. No Fourier support was silently left
unchanged after the field map. The two estimates have
been placed in a common, explicitly nonempty domain.

This comparison is an actual restricted classical
stationary-action result. It is neither a loop momentum
cutoff nor a construction for the entire rolling bounce.
