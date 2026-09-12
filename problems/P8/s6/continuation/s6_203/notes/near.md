# Complementary near-momentum integral

NEAR consists of |k|>=m and |P|>delta nu(k)/2. Hence nu<2|P|/delta and |k|<L=2Amax|P|/delta. If L<=m the region is empty. The actual upper radius Amax sqrt((2|P|/delta)^2-m^2), when defined, is smaller thanL.

The all-real row2E nu and the near volume yield

integral_NEAR ||E5(P,k)|| <=(16/3) E Amax^3 delta^-4 |P|^4/pi^2.

For each fixed |k|>=m the S6.202 Cauchy estimate bounds its qth homogeneous Taylor polynomial by E delta^-q |P|^q nu^(1-q), q0,...,4. This estimate on the polynomial coefficients remains valid when the evaluation point is outside the small analytic ball. It does not assert analytic continuation of the full endpoint to arbitrarily large complex P.

On m<=k<=L, use nu<=2k for q0, the exact nu^0 for q1, and nu>=k/Amax for q2,...,4. The radial upper bounds, including measure1/(2pi^2), are respectively

L^4/(4pi^2),
L^3/(6pi^2),
Amax L^2/(4pi^2),
Amax^2 L/(2pi^2),
Amax^3 log(L/m)/(2pi^2).

All five orders are retained, especially the logarithmic fourth order. SubstitutingL and adding the raw row gives a coefficient of |P|^4 plus the log term. Bound log(L/m)<=L/m onL>m and pi^2>9. Thus the complete near remainder is at mostC4|P|^4+C5|P|^5, with

C4=E delta^-4 [4Amax^4+(26/3)Amax^3]/9,
C5=E Amax^4/(9delta^5m).

Numerically C4<1.265e52, C5<1.325e54, and C4+C5<2e54. These decimal displays round up; the executable constants are exact rationals.

The original both-created-mode indicator is retained throughout and only decreases the absolute majorant. No moving boundary or Fourier test is Taylor-differentiated. Independent hyperbolic-radial quadrature checks every Taylor order and the complete raw-plus-polynomial near integral at low through very large external momenta.
