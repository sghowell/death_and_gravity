# Exact Gaussian masters and dimensionless mass scaling

Let D=d=4-2epsilon and retain the same
exp(gamma_E epsilon) mu^(2epsilon) convention per loop.
After the numerator reduction, all terms are massive/massive/
massless vacuum integrals with integer powers alpha,beta,gamma.
For positive Schwinger exponents in a common convergent domain,

 I_D(alpha,beta,gamma)
 =exp(2gamma_E epsilon) mu^(4epsilon)/Q^2
  *(m^2)^(d-alpha-beta-gamma)
  *Gamma(d/2-gamma) Gamma(alpha+beta+gamma-d)
  *Gamma(alpha+gamma-d/2) Gamma(beta+gamma-d/2)
  /[Gamma(alpha) Gamma(beta) Gamma(d/2)
    Gamma(alpha+beta+2gamma-d)].

The two-loop Gaussian determinant is tu+v(t+u).
Integrating v first gives B(gamma,d/2-gamma).
Then t+u=T and t/T=z give the remaining Gamma and Beta factors.
This derives the normalization, including Gamma(d/2) in the denominator.
Continuation of this explicitly meromorphic expression defines the
same dimensional-regularization masters outside that convergence domain.

Use x=A-m^2, y=B-m^2 and
k dot l=(A+B-C-2m^2)/2, where C=(k-l)^2.
The trace tensors reduce to a finite list of these masters.
Nonpositive alpha or beta give scaleless polynomial massless integrals
after the remaining loop momentum is shifted, and vanish in this
regulator. For nonpositive gamma the formula is continued to the
factorized tadpole moments; it is not used as a convergent Schwinger
integral with a negative Schwinger exponent.

At m=mu=1 each master is Gamma(epsilon)^2 times an exact rational
function. Every shift is implemented with Gamma(z+n)/Gamma(z) and
finite products of exact SymPy integers/rationals. No floating
rational reconstruction is used. The independently checked case

 I_D(1,1,1)=Gamma(epsilon)^2/
            [(1-epsilon)(2epsilon-1)]

has its common loop/mass factors stripped. The four-dimensional
fractional-power sunset used for the finite mass difference is
separately convergent and does not require this continuation.

Restoring m, the raw physical p^2 coefficient has scale
(mu/m)^(4epsilon). The one-loop counterterm insertions have scale
(mu/m)^(2epsilon). Both become one at mu=m, but they must not be
identified before their fixed-mu mass derivatives are taken.
