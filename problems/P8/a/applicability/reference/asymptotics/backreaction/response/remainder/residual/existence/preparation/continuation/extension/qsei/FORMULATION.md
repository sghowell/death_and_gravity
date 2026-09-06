# P8-A.15 — fresh QSEI on the extended actual SEE geometry

This child pins the actual same-source A.14 solution and the general
positive-type/C1 scattering lemmas of A.12. It does not transfer A.12's
short-interval numerical bound or positive reference credit to a new
domain. All new constants and their hypotheses are recalibrated here.

## Frozen physical inputs

The metric and reference state are the actual smooth A.14 solution of
the full density and trace semiclassical Einstein equations (SEE), with
the same original radiation normalization, massless minimally coupled
real scalar, transported original vacuum and finite renormalization
prescription. Use the FK +--- curvature/sign convention, physical
epsilon=1, delta=10^-14, lambda=2sqrt(2) A eta_star^2 and gamma=0.
The physical scale is T0=A eta_star^2 with
`kappa hbar=2880 pi^2 delta T0^2`.

Dimensionless conformal time is x; proper time is tau=T0*s and ds=a dx.
Write h=a'/a and u=-a''/a. The A.14 future interval has length
L=10^-6, but the **unchanged** preparation window is L0=10^-10.
The external source is zero from x=x0+L0/2 onward. The new theorem's
sampling domain is the open interval

    I=(x0+L0/2,x0+L),

not (x0+L/2,x0+L), nor the whole prescribed old target. Every sampler
has compact support in I; real H2_0 samplers on compact subintervals
are allowed. The reference is homogeneous, but target states need
not be homogeneous, quasifree, zero mean, or solutions of the SEE.

## Certified analytic statement

For every Hadamard target state of this field on the fixed actual
smooth spacetime, and every such proper-time sampler f,

    integral f^2 E_target d_tau
       >= -2 hbar/(16 pi^2) integral |f_proper_second|^2 d_tau.

E is the effective energy density (EED) in the named prescription.
This is an absolute, not just reference-subtracted, QSEI. Its new
reference estimate is signed:

    E_reference >= -hbar/(40 pi^2 T0^4).

No positive reference EED is asserted on the longer domain. Weighted
metric gains from the actual A.14 fixed-point distance prove this
lower bound. Proper-time Poincare bounds absorb its possible negative
part into an additional `(2/5)L^4` in the second-derivative coefficient.
The freshly computed difference coefficient plus that penalty is
strictly below 2. No bound on higher potential derivatives is assumed.

## Actual focusing boundary

For the comoving normals to constant cosmic-time surfaces, in either
orientation, on segments contained in I, every admissible index trial
obeys `J[g]T0>=3/tau_hat-tau_hat/8>3/4>=-K T0`. The normal volume ratio
is at least 8/27. Thus the sufficient comoving index trigger cannot be
met on this certified interval, even using its actual curvature.
The ratio of the chosen sufficient QSEI coefficient to the square of
the available duration is at least 2/5, improved from the old short-slab
ratio but not a substitute for the actual index test.

This is not an optimal-coefficient theorem, an exclusion for other
hypersurfaces or unknown continuations, a realistic massive/interacting
field result, or a proof of cosmological incompleteness. Original P8(a)
and P8 remain open. The mathematical solution interval does not itself
assess stress fluctuations or high-derivative preparation EFT validity.
