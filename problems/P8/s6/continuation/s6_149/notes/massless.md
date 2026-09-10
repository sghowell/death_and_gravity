# Dimensional tensors, proper fermion insertions and fixed-scale checks

Let e=epsilon, d=4-2e, Q=16pi^2 and A=exp(gamma_E e)Gamma(1+e).
The vacuum theta graph contains two massive fermion lines and one
boson line, with weight one half. At m=1 its scalar trace is
tr(S_k S_l)/2; the gauge trace is
-tr(gamma_mu S_k gamma_mu S_l)/2. The latter minus sign includes
the two gauge vertex phases. Dimension-symbolic Gaussian reduction
gives, in Gamma(e)^2 units,

    V_s=4/[(1-e)(2e-1)]+1/(1-e)^2,
    V_g=-4/[(1-e)(2e-1)]+2/(1-e).

Scaleless masters vanish only in the prescribed dimensional regulator.
The two proper fermion self-energy counterterms combine into the
variation of the one-loop fermion determinant. Its vacuum expression is

    V1,F=2N exp(gamma_E e)mu^(2e)Gamma(e-2)m^(4-2e)/Q.

A kinetic rescaling changes the effective mass by (eta-z)m.
The constant Jacobian is a scaleless trace. Hence its variation is
(eta-z)m partial_m V1 at fixed mu, with eta,z exactly from S6.139.
In scalar NYm^4/Q^2 and gauge NaC_Fm^4/Q^2 units, the normalized
proper counterterms are c_s=6/(1-e), c_g=-12/(1-e).
The complete massless-exchange forests at mu=m are

    [A^2 V_s+A c_s]/e^2=3/e^2-4/e-19+O(e),
    [A^2 V_g+A c_g]/e^2=-6/e^2+2/e+18+O(e).

The pi-squared terms cancel only after the full proper forest
and its finite epsilon-times-pole products are combined.
MS subtracts the two poles, not the finite constants.

For general L=log(m^2/mu^2), the finite polynomials are

    scalar: -3L^2+14L-19,
    gauge: 6L^2-16L+18.

Their fixed-mu second mass derivatives at mu=m give the previously
computed quadratic mass coefficients -56 and +40. Independently,
the full regulated derivative multiplies the raw vacuum by
(4-4e)(3-4e), but the proper-counterterm part by (4-2e)(3-2e).
Both routes agree with S6.147. Differentiating only the final
mu=m constants would give the wrong result.

The scalar massless chord is only an auxiliary reference.
The gauge chord is actually massless. The scalar vacuum has N=6
active color states; the gauge vacuum has N=42.
