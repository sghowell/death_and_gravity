# Subtracted dispersion with an unfixed local sector

Write rho_i for the positive densities derived from the full cut.
The exact identity

    1/(s-z) = 1/s+z/s^2+z^2/s^3+z^3/(s^3(s-z))

isolates the three-subtracted dispersive part D_i. It does not define
the physical finite coefficients of the subtracted terms. Ward
identities, tensor contacts and the fixed covariant prescription constrain
their admissible local combination; they are not independently adjustable
counterterms chosen in this continuation.

Since rho_i=O(s^2), the D integrand is O(s^-2) at infinity. On each compact
subset of C minus [4m^2,infinity), |s-z| has a positive lower bound near
threshold and is uniformly comparable to s at infinity. The integrand
and its z derivatives are locally dominated. Integrating, or applying
Morera with dominated Fubini, proves analyticity. Threshold beta behaves
as sqrt(s-4m^2), which is integrable. With only two subtractions the
large-s absolute integrand is instead asymptotic to a nonzero multiple
of 1/s for z!=0, and the corresponding absolute ultraviolet integral
does not converge.

For |z|<=2m^2 and s>=4m^2, |s-z|>=s/2. Thus

    |D_i|<=2|z|^3 integral rho_i(s)/s^4 ds.

Subtracting z^3 I_i exactly gives
z^4 integral rho_i(s)/(s^4(s-z)) ds and hence the stated next-order
bound with J_i. These bounds are finite and apply to complex z in
the entire closed low-energy disk, not only its real diameter.

Set v=sqrt(1-4m^2/s). Then
s=4m^2/(1-v^2), ds=8m^2 v/(1-v^2)^2 dv.
Each I or J becomes a polynomial integral on [0,1], yielding

    I2=3/(3584 pi^2 m^2), I0=3/(8960 pi^2 m^2),
    J2=31/(322560 pi^2 m^4), J0=1/(32256 pi^2 m^4).

All are positive; I2/I0=5/2. Independent high-precision integrals
reconstruct these constants directly from the spectral measure and test
both D inequalities for positive, negative and complex z.

The exact moment-order and spin types are validated before memoization.
Otherwise Python cache equality can permit a floating-point request to
reuse an integer result without validation. A prewarmed-cache negative
test retains this distinction.
