# A genuine individual-coefficient obstruction to naive Gaussian substitution

The complete source inequalities on0<N<=1 imply
(2h)^(-1/4)N^(-3/2)<=D=M/N<=2^(1/4)N^(-3/2)
and
(2h)^(-3/4)N^(-1/2)<=Q=N/U<=2^(3/4)N^(-1/2).
These are the original lapse weights whose full parent bindings
are checked in notes/source, not new simplified action coefficients.

For a centered Gaussian n with variance sigma^2>0, the density of
N=1+n on0<N<1 is bounded below by
c=exp[-1/(2sigma^2)]/(sqrt(2pi)sigma)>0.
The actual original variance has the evaluated positive floor proved
in notes/covariance. It is not a formal zero-variance limiting state.

The exact cutoff comparison integrals are
integral_epsilon^1 N^(-3/2)dN=2(epsilon^(-1/2)-1)
and integral_epsilon^1 dN/N=-log(epsilon).
Both diverge as epsilon decreases to zero. Multiplication by the
very small but strictly positive c changes neither divergence.
Therefore the positive-lapse Gaussian expectation of D is infinite,
and the second positive-lapse Gaussian moment of Q is infinite.

The first moment of Q is finite. Near zero it is bounded by a
constant times N^(-1/2), whose integral is finite. For N>=1,
R<=2 gives Q<=2N; the remaining Gaussian tail is integrable.
Finite first moment is not finite second moment.

These statements use positive Borel spectral functions of the
original single LINEAR Gaussian lapse. Setting their observables
zero for N<=0 leaves the divergence as N approaches zero from above.
It neither projects nor conditions the state.

This is a precise obstruction to evaluating EVERY exact nonlinear
coefficient by unrestricted Gaussian substitution. It is NOT proof
of divergence of the complete constrained Hamiltonian/action, an
inconsistent interacting state, a candidate no-go or original-P8
closure. The full Hamiltonian and all its other channels remain.
Nonlinear auxiliary reconstruction, full-channel cancellations,
a physical measure and the EFT validity domain have not been supplied
by this substitution and are not assumed absent.

The original local regular branch does not authorize extrapolation
of a Gaussian lapse to all full-action histories. Global positivity
of R alone does not keep every auxiliary or Legendre pivot regular.
These domain and quantum-construction questions remain open.

Independent finite-localizer diagnostics integrate the full coefficient
with successively smaller positive-lapse cutoffs. A reciprocal-square
variable exposes the kinetic square-root growth, and a logarithmic
variable exposes the canonical-matter second-moment growth. The same
full coefficient has a convergent first-moment integral. These
moderate Gaussian fixtures are diagnostics; the proof never evaluates
the actual tiny density prefactor as floating-point zero.
