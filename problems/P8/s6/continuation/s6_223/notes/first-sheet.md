# First-sheet zero and exact rational enclosure

Write c=1-y^2 and W=W2. On0<=y<=1,

30-20y^2+3y^4=13+(1-y^2)(17-3y^2)>=13.

Thus W is positive except at the measure-zero endpointy0. On every compact subset of the slit plane, the radial integrand has an integrable uniform majorant, so A2=1/30+integral pW/(4m^2+cp) is analytic. For p=a+ib,

Im[p/(4m^2+cp)]=4m^2 b/[(4m^2+ca)^2+c^2b^2].

Its integrated sign is strictly the sign of b. Hence there is no nonreal zero. For real p>-4m^2 the derivative is the integral of4m^2W/(4m^2+cp)^2 and is strictly positive. Direct integration, with dominated convergence at the threshold, gives

A2(0)=1/30,
A2(-4m^2)=-172/225,
m^2 A2'(0)=3/56.

Therefore exactly one real zero p=-r m^2 occurs with0<r<4. It is simple; no additional real zero can occur on(-4m^2,infinity). The open-cut banks have a nonzero imaginary part as proved separately, and the threshold is not a zero. This accounts for the entire finite first sheet. The reciprocal tends to zero at infinity.

For an exact enclosure set x=r/4 and

M_j=integral_0^1 y^2(30-20y^2+3y^4)(1-y^2)^j dy>0.

Since M_j<=M_0, with N=16,

A2(-r m^2)=1/30-r/120 sum_(j=0)^N M_j x^j - E,
0<E<=r M_0 x^(N+1)/[120(1-x)].

All moments and endpoints are exact rationals. At rlo=2899/5000 the finite expression minus its tail bound is positive. At rhi=5799/10000 the finite expression is negative. Consequently rlo<r<rhi. The certificate retains the exact interval endpoints and all17 moments; no sampled plot or numerical root is used for this conclusion.

Similarly,

m^2 A2'(-r m^2)=1/120 sum_(j=0)^N (j+1)M_j x^j + E1,
0<E1<=M_0 x^(N+1)[N+2-(N+1)x]/[120(1-x)^2].

This derivative increases with r. The lower series bound at rlo and the upper bound at rhi therefore enclose its value at the root. Taking reciprocals proves16<R/m^2<17 for R=1/A2'(-r m^2). R is the positive pole weight of1/A2. The residue of1/F2=-1/A2 is-R.

An independent test expands (1-y^2)^j binomially, integrates monomials as exact rationals and uses a shorter N=12 enclosure. It proves the same coarse root and weight bounds without production moment integration. A high-precision root is used only to test the dispersion formula numerically.
