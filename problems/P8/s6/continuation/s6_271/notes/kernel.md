# Direct Gaussian frame and explicit Weyl operator bound

Use the normalized unit-CCR coherent functions in d configuration variables

    phi_(q,p)(x) = pi^(-d/4) exp(-|x-q|^2/2 + i p.(x-q/2))

and the defining Weyl kernel

    [OpW(a)psi](x) =
      (2pi)^(-d) integral exp(i(x-y).xi) a((x+y)/2,xi) psi(y) dy dxi.

The coherent resolution of the identity has measure
dmu(z) = (2pi)^(-d) dz. Its analysis map W is an isometry. This infinite
phase-space measure is not itself a probability measure.

## Exact matrix-element normalization

For one pair, change variables to X = ((x+y)/2,xi) and eta = x-y.
The exponent of the two coherent factors times the Weyl phase is quadratic
in eta with coefficient -1/4. Its linear and constant coefficients B,C obey

    C + B^2 =
      -(X_q-(q+q')/2)^2 -(X_p-(p+p')/2)^2 + i Phi,

    Phi = (q-q')X_p - (p-p')X_q + (p q' - q p')/2.

The eta Gaussian integration supplies 2 sqrt(pi). Together with the
original (2pi)^(-1) Weyl and pi^(-1/2) seed factors, its prefactor is 1/pi.
Taking products proves the full identity

    <phi_z, OpW(a) phi_z'> =
      pi^(-d) integral a(X) exp(-|X-(z+z')/2|^2) exp(i Phi) dX.

For a = 1 it gives the coherent overlap with amplitude
exp(-|z-z'|^2/4) and the retained symplectic phase; at z = z' it gives 1.
The code checks the exact complete exponent, center, phase and prefactor.
Independent Gaussian quadratures compare the defining Weyl kernel and this
formula, including off-diagonal complex phases, not just the diagonal.

## Product integration by parts

The phase gradient in X is an orthogonal symplectic rotation of z-z'.
Apply the full product over the 2d phase coordinates of (1-partial_j^2)
to a(X) exp(-|X-m|^2), with m = (z+z')/2, and integrate by parts.
The oscillatory exponential acquires product_j(1+Delta_j^2), where
Delta is that orthogonal rotation.

For the normalized one-dimensional Gaussian pi^(-1/2) exp(-x^2), the
L1 norms of derivatives of orders 0, 1 and 2 are bounded by 1, 2 and 4.
Indeed E|x| <= (E x^2)^(1/2), and
|4x^2-2| <= 4x^2+2, with E x^2 = 1/2.
The complete product rule in each coordinate therefore has the weights

    derivative order of a:   0    1    2
    upper-bound weight:     5    4    1.

The first weight includes both the identity term and the second derivative
of the Gaussian; the middle includes both product-rule contacts.
Multiplying these bounds includes EVERY mixed multiindex alpha_j in {0,1,2}.
Thus the coherent kernel satisfies

    |K(z,z')| <=
      sum_alpha c_alpha ||partial^alpha a||_infinity
      * product_j (1+Delta_j^2)^(-1).

No supremum-of-symbol assertion is substituted for this kernel bound.

## Both Schur integrals and extension

At fixed z, translation and the orthogonal symplectic rotation preserve
Lebesgue measure. Since integral_R (1+t^2)^(-1) dt = pi,

    integral product_j(1+Delta_j^2)^(-1) dmu(z') = (pi/2)^d.

The other Schur integral has exactly the same bound. Applying the elementary
Schur test on the coherent L2 space and the isometry W gives

    ||OpW(a)|| <= (pi/2)^d sum_alpha c_alpha ||partial^alpha a||_infinity.

The coefficient sum is 10^(2d). Its total-order coefficients are exactly
those of (5+4t+t^2)^(2d); at d = 48 the code retains all 193 coefficients.
Since pi < 4 and 2^48 < 10^15,

    (pi/2)^48 * 10^96 < 10^112.

This proves the explicit fixed-dimension estimate used here. The largest
total derivative order needed is 192, and all multiindices are retained.

First apply this proof to smooth compact symbols. It also applies directly
to Schwartz symbols; integrations by parts have no boundary contributions.
For a scalar plus a Schwartz symbol, treat the scalar identity exactly.
On Schwartz vectors the frame reconstruction agrees with the defining Weyl
operator, so the estimate provides its unique bounded extension.
Real symbols give bounded self-adjoint operators. Smooth compact and
Schwartz phase symbols have Schwartz Weyl kernels. This is an explicit
proof for the needed class, not an unspecified-dimensional
Calderon-Vaillancourt citation.

## Primary convention audit and independent checks

The author-hosted Didier Robert text, *Propagation of Coherent States in
Quantum Mechanics and Applications*, section 1, pp. 6-9, was consulted for
the defining Weyl kernel and coherent-frame approach:

https://www.math.sciences.univ-nantes.fr/~robert/proc_cimpa.pdf

The inspected PDF's displayed Eq. (29) prints a factor 2^(2d), inconsistent
with its normalized diagonal-state test when combined with Eq. (28).
The later displayed frame inversion also omits the normalization factor.
Those printed prefactors are NOT assumed here. The preceding calculation
derives both the cross-Gaussian prefactor and frame measure independently.
The numerical tests explicitly reject an extra dimension-dependent factor.

No later theorem from that text is required for the exact constant used in
this continuation. Numerical integration fixtures, including the full
real-axis Schur integral, corroborate the written proof but do not turn it
into a formally verified theorem or a full P8 evolution simulation.
