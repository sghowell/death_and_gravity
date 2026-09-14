# Strict localization effect and state-dependent comparison

Consider 0<=chi<=1, smooth, compactly supported and equal1 on
a nonempty core in finite-dimensional phase space. Since the
coherent vectors have unit norm, Q_V(chi) is positive trace class
with trace (2pi)^(-d) integral chi. For a Gaussian window the
coherent transform is, after its invertible squeezing and a
nonzero Gaussian factor, a Bargmann entire function. A nonzero
L2 vector therefore cannot have coherent transform zero on
a whole open phase set: analytic continuation would make the
entire function, and then its vector by the isometry, zero.

Apply this fact to an open set on which chi>0 and to an open
set outside its compact support. For every nonzero psi,

    <psi,Q_V(chi)psi> > 0,
    <psi,(I-Q_V(chi))psi> > 0.

The nonzero positive compact operator has its norm as an
eigenvalue. No eigenvalue equals1 by the second inequality;
thus 0<||Q_V(chi)||<1. Its eigenvalues tend to0, so
||I-Q_V(chi)||=1. The effect is not a projection. There is
no nonzero exact chart-supported state in this coherent-effect
sense. This does not introduce a forbidden classical joint
probability law for noncommuting q,p.

Whiten with a symplectic S satisfying V=SS^T/2, z=S w.
One choice is the positive symplectic square root of2V.
The translation group acts orthogonally and preserves V,
so its action commutes with that square root and preserves
||w||. The SAME Gaussian Husimi density becomes the standard
normal density in2d real dimensions, not the Wigner density
with covariance I/2.

For the sharp ball ||w||<=R set t=R^2/2. In standard coherent
coordinates alpha=(w_q+i w_p)/sqrt(2), rotational integration
makes Q(1_ball) diagonal in total occupation m. Integrating
the coherent monomial of total degree2m first over the unit
sphere and then radially gives

    lambda_m = gamma(d+m,t)/Gamma(d+m)
             = 1-exp(-t) sum_{k=0}^{d+m-1} t^k/k!.

The multiplicity is binomial(d+m-1,m). The exact difference
lambda_m-lambda_(m+1)=exp(-t)t^(d+m)/(d+m)! is positive, so
lambda_0 is the top eigenvalue. For finite R>0 every eigenvalue
lies strictly between0 and1. The vacuum itself has eigenvalue
lambda_0; its leakage is

    tail_d(R)=exp(-t) sum_{k=0}^{d-1} t^k/k!.

The normalized phase-ball trace is t^d/d!. The sharp ball is
a diagnostic effect, not a replacement for the smooth cutoff
needed in the calibrated Hamiltonian.

If chi is1 inside R1 and0 outside R2, operator positivity gives
lambda_0(R1)<=<psi0,Q(chi)psi0><=lambda_0(R2). For a radial
chi the vacuum is itself an eigenvector, by the same angular
calculation. No large-probability statement follows without
the actual WHITENED radii and full dimension.

Two useful dimension-dependent bounds follow directly. Dropping
exp(-s)<=1 in the lower incomplete gamma integral gives
1-tail_d(R)<=t^d/d!. For t>d, exponential Markov applied to
a unit-rate gamma variable gives
tail_d(R)<=exp[-t+d+d log(t/d)] by optimizing the positive
Laplace parameter. For fixed radius and growing d the lower
bound already warns against an automatic continuum limit.

There is nevertheless a useful operator-on-state comparison.
If a bounded symbol delta_a is zero on the whitened core radiusR
and |delta_a|<=M, the coherent compression inequality gives

    ||Q_V(delta_a)psi0||^2
      <= <psi0,Q_V(|delta_a|^2)psi0>
      <= M^2 tail_d(R).

This uses the exact full seed Husimi density, not a fictitious
sharp state projection. For calibrated cutoffs delta_a includes
all gradient and Hessian contacts. It vanishes in a common core
only because each cutoff is identically1 there.

For later times an additional estimate is needed. Let bounded
generators A(u),B(u) act on the SAME Hilbert space, with
delta(u)=A(u)-B(u), ||delta(u)||<=M, ||B(u)||<=B0, and
||delta(u)psi0||<=M sqrt(tail) uniformly. Exact Duhamel gives

    (U_A(t)-U_B(t))psi0
      = -i integral_0^t U_A(t,s) delta(s) U_B(s,0)psi0 ds.

Split U_B(s)psi0 into psi0 plus its difference, whose norm is
at most B0|s|. Unitarity then yields

    ||(U_A(t)-U_B(t))psi0||
      <= M sqrt(tail)|t| + M B0 |t|^2/2.

The corresponding statement for negative time uses the
appropriate suprema on the reversed interval. All quantities
refer to generators in units hbar=1; restoring hbar divides
Hamiltonian bounds accordingly. This is a comparison between
two defined bounded extensions, not a comparison with the
undefined original singular Hamiltonian. Initial Gaussian
leakage alone is not later-time support or regulator removal.
