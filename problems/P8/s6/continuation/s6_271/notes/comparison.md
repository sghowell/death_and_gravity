# Heat remainder to a genuine operator norm and positive Weyl volume

Work in the actual unit-CCR symplectic whitening of the original pure seed.
The normalized same-seed coherent map satisfies Q(a) = OpW(exp(D)a),
D = Delta_w/4. This identity was established with its Gaussian normalization
in the parent chain. It is not a claim that Q preserves its window under
arbitrary symplectic changes.

For either complete extended symbol f_ext = f0 + chi_c(f-f0), compare

    F_C or A_C = Q[(1-D)f_ext],
    F_W or A_W = OpW(f_ext).

The scalar f0 is present in both operators and in each actual propagator.

## Exact heat identity

For the commuting heat generator D,

    d/dt [exp(tD)(1-tD)f_ext] = -t exp(tD)D^2 f_ext.

Integrating from zero to one gives

    r_f = exp(D)(1-D)f_ext - f_ext
        = -integral_0^1 t exp(tD)D^2 f_ext dt.

Here f_ext is scalar plus smooth compact; its derivatives and the heat
convolution are in the class needed for these identities. The scalar drops
out of D^2. The resulting remainder is Schwartz, so the operator theorem
in notes/kernel.md applies directly.

Every partial derivative commutes with the Euclidean heat semigroup, whose
convolution kernel is positive and has total mass one. It therefore
contracts the supremum norm. Expanding both complete Laplacians in D^2,
including all ordered coordinate pairs, gives for total order n <= 192

    ||partial^alpha r_f||_infinity
      <= (m^2/32) J_(n+4)(A),   m = 96.

The factor is (1/16) from D^2 times (1/2) from the t integral, with m^2
ordered fourth derivatives. This includes every cross derivative.
By notes/derivatives.md, J_(n+4) <= J_4 through order 196.

## Actual operator estimate

Applying the complete product-Schur estimate, not a symbol-sup shortcut,

    ||OpW(r_f)|| <= K (m^2/32) J_4(A),   K = 10^112.

The exact same-seed identity identifies this as ||F_C-F_W|| or ||A_C-A_W||.
The code evaluates the rational quantities using

    J_4(A) = 2A * 2056^4 * (4!)^3 / R^4,   R = 10^20.

For the full physical volume A = E = 10^-255, the result is less than
epsilon_F = 10^-200. For the complete interaction A = H = 10^1000,
it is less than epsilon_H = 10^1056. The two exact errors have identical
relative ratios to E and H, respectively, and that identity is checked.

The high power in epsilon_H is not hidden: the evaluated time interval
T = 10^-2000 is what makes its entire unitary effect small. These are
fixed-family bounds, not uniform estimates in a limiting regulator family.

## Concrete positive volume, not a positive Weyl map

The parent calibrated-coherent volume satisfies

    ||F_C-I|| < 2E.

Since both compared operators are bounded and self-adjoint,

    ||F_W-I|| <= 2E + ||F_W-F_C|| < 10^-199,
    F_W >= [1-2E-||F_W-F_C||] I > I/2.

This is an operator inequality on the full Hilbert space and, by symmetry,
on the original zero-charge reducing sector. The source is the full
normalized physical volume, not a quadratic surrogate.

It would be false to replace this argument with "positive symbol implies
positive Weyl operator." Independent one-pair tests use the positive
Gaussian symbol exp[-a(q^2+p^2)], a > 1. Its first excited Weyl eigenvalue
is (1-a)/(1+a)^2 < 0, obtained from the radial Wigner integral.
The actual near-identity volume positivity proved here is compatible with
that counterexample.

Additional Gaussian-radius fixtures compare exact Weyl and calibrated
coherent spectra and check a separately evaluated one-pair Schur error
bound. They diagnose the comparison mechanism; they neither use a finite
occupation truncation to define the theorem nor simulate the full P8 state.
