# Constant-field determinant, all flavors and finite domain

Use the Dirac determinant in a gauge-invariant dimensional
regulator. At this order the gauge field vanishes. Each
color of each Dirac flavor contributes

    -m_j^4[log(m_j^2/nu^2)-3/2]/Q.

The active masses are mF+y Phi and mF-y Phi. The other
twelve masses are mF, independent of Phi. Therefore their
contribution to the vacuum constant must not be discarded.
At nu=mF the total constant is 63mF^4/Q; the active part
alone would incorrectly give 9mF^4/Q.

Before subtraction the entire dimensional pole/reference
is

    Ibar {Nc[(mF+y Phi)^4+(mF-y Phi)^4]+12Nc mF^4}/Q.

Its constant, second and fourth derivatives at zero are
42mF^4 Ibar/Q, 12NYmF^2 Ibar/Q and 24NY^2 Ibar/Q.
The second agrees with the independent two-point trace.
No four-dimensional norm bound is applied to this
unsubtracted divergent determinant.

The finite derivatives at zero are f(0)=4NYmF^2/Q and
v4=-64NY^2/Q. Odd derivatives cancel by active flavor
exchange, which is the candidate's exact Phi parity.

Put r=y Phi/mF. The dimensionless active coefficients
at degrees zero, two and four are 9,12,-16. For even
k>=6 they are

    c_k=96Nc/[k(k-1)(k-2)(k-3)(k-4)] > 0.

To derive this for every k, differentiate
z^4(log(z^2/nu^2)-3/2) five times: the result is 48/z.
Every further derivative is
48(-1)^(k-5)(k-5)!/z^(k-4); the two flavor signs add
at even k, and the determinant sign reverses the negative
derivative. Divide by k! and multiply by 2Nc.
The independently differentiated degrees through fourteen
are regression checks, not the proof of the infinite tail.

Both logarithms have convergent real series on |r|<1.
There c_k decreases for even k>=6, with c_6=2/5.
The remainder is nonnegative and bounded by
(2/5)r^6/(1-r^2), and by 8r^6/15 on |r|<=1/2.
This bound covers every higher field degree of this
one-loop determinant. It does not cover higher loop orders
or any point where a fermion mass crosses zero.
