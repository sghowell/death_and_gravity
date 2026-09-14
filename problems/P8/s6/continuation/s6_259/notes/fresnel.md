# Entire finite Fresnel integral and configuration density

At a finite cell let lambda>0 be its full action scale, including
kappa and any physical cell-volume factor. Let A be the complete
real symmetric 56-dimensional complement Hessian. Define the
oscillatory integral by

    integral d^56 z exp(i lambda z^T A z/2 - epsilon z^T z/2),

with epsilon>0, and take epsilon down to zero after diagonalizing
the finite real quadratic form. Each real eigenvalue contributes
the branch continued from its convergent Gaussian. Thus

    I(A,lambda)=(2 pi/lambda)^28 |det A|^-1/2
                 exp[i pi(n_positive-n_negative)/4]
              =-(2 pi/lambda)^28 / sqrt(det A).

The exact (26,30) inertia fixes the phase. Since (-i)^56=1 and det A>0,
the single endpoint determinant det(-i lambda A) is positive.
Its principal square root would give the wrong sign. Independent
convergent numerical Gaussians check the continued full 56-dimensional
limit; the phase is not inferred from numerical sampling.

The positive configuration density and separate inverse phase are

    rho_abs=(lambda/(2 pi))^28 sqrt(det A),
    phase_inverse=exp[-i pi(26-30)/4]=-1.

Their product with the entire integral is exactly one. The positive
constraint density in phase space is not itself a complex phase.
For the opposite CTP branch use the conjugate continuation; the
original retained boundary state and final trace are unchanged.

This gives a normalized conditional insertion into any separately
specified retained finite functional. It does not choose that
functional, an interacting state, or a regulator for its dynamics.
All cancellations here occur at the same finite regulator, before
any continuum determinant is defined.

## Flat measure is a different finite prescription

A flat Lebesgue configuration integral without rho_abs retains both
lambda^-28 and det A^-1/2. At fixed lambda, the normalized determinant
relative to R=1 is

    det A(R)/det A(1)
      = R^(65/2) (4sqrt(R)+5)^3 (2R-1)^3 (4-R^(3/2))/2187.

Its first and second logarithmic R derivatives at R=1 are 116/3
and -1228/27. The corresponding flat Gaussian log weight has minus
one half of those derivatives. This coefficient dependence does
not vanish simply because the connection variables were auxiliary
in the classical equations.

More generally,

    d log I=-(56/2)d log lambda-(1/2)Tr(A^-1 dA),

with the full mixed second variation retaining both d1d2 A and
-A^-1 d1 A A^-1 d2 A, as well as d1d2 log lambda.
If lambda contains a varying physical volume, its derivatives cannot
be dropped when evaluating lapse or metric variations. The fixed-lambda
coefficient derivatives above are not such a physical stress calculation.

For any invertible real change z=B zprime at the fixed history,
Aprime=B^T A B, bprime=B^T b, and dz=|det B| dzprime.
Hence sqrt(|det Aprime|)=|det B|sqrt(|det A|) and
bprime^T Aprime^-1 bprime=b^T A^-1 b. Sylvester inertia is unchanged.
The normalized insertion therefore transforms consistently even for
a background-dependent basis; no flat coordinate-density assumption
is hidden by switching to an orthonormal frame.

Do not set a formal delta(0) to zero to infer a physical continuum
stress, and do not cancel separately regularized functional determinants
without a prescription that establishes the required identities.
