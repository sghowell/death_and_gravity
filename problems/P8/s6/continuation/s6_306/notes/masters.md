# Exact whole mixed-box collapse and finite calibration

For A_b(v)=1-bv(1-v), the frozen S294 box remainder is

K(a,b)=integral_0^1 dv integral_(u,h>=0,u+h<=1) du dh
 u(u+2h)/[u^2 A_b+h(n-a(1-u-h))]^2.

First take b<4 and a<n so all denominators are positive. Put

y=h/u^2, w=u+h=u+yu^2.

For fixed y>=0 the map u to w is monotone, with derivative1+2yu;
the original triangle becomes0<=w<=1,0<=y<infinity. The first Jacobian
is dh=u^2dy. The denominator is u^2[A_b+y(n-a+aw)]. The full original
numerator and both Jacobians leave exactly

K=integral dv dw dy/[A_b+y(n-a+aw)]^2
 =integral dv/A_b * integral_0^1 dw/(n-a+aw)
 =J0(b)[ln n-ln(n-a)]/a.

The a0 limit is J0(b)/n. Positivity justifies the changes of variables
and integrations on the initial domain. Analytic continuation of the
entire identity gives the normal-sheet physical boundary. It does not
license taking absolute values of an undeformed singular real integral.
Thus every one of the six ordered boxes in S294 has

W(a,b)=T(b)-aK(a,b)=T(b)+J0(b)ln(1-a/n).

## Independent full finite coefficient

Write c_nu=EulerGamma-ln(4pi nu^2). Expanding the entire S294
D=4+2e box gives

D0=J0/(2(n-a)e)+[J1/2+c_nu J0/2+W]/(n-a)+O(e).

For the primary Box16 formula put
x=(sqrt(1-4/(b+i0))-1)/(sqrt(1-4/(b+i0))+1),
h=2/(sqrt(n)+sqrt(n-4)), q=n-a.
Then the residue identity is J0=-2x ln x/(1-x^2). After translating
the dimensional convention and subtracting the displayed J1 and c_nu
terms, its finite remainder is

W=-x/(1-x^2){2ln x[ln(q/sqrt(n))+ln(1-x^2)]
 +2ln^2 h+Li2(x^2)+Li2(1-xh^2)+2Li2(1-x)+Li2(1-x/h^2)-pi^2/6}
 -J1/2.

The analytic-continuation correction in the multiargument dilogarithms
vanishes here because h is positive and its powers introduce no extra
phase. Physical x is on the upper bank of the negative axis, so the
Li2 arguments1-xh^(+/-2) are on their lower banks. Arbitrarily joining
principal logarithms would not be valid in a general mass assignment.

Ten finite comparisons at subthreshold and physical b, with two distinct
a values, independently agree with the projective answer. The full J1/2
is essential; dropping it is an explicit negative control. The normal
sheet also fixes the nonzero imaginary parts. Two direct integrations
of the original u,h integrand supply a different numerical check. Their
sharp h~u^2 A/n corner is partitioned explicitly. Numerical comparisons
supplement, and are not the proof of, the general analytic identity.
