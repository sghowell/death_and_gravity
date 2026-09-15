# Whole offshell mixed amplitude and literal graph reduction

Original V/G/B/P8 stays OPEN. These are the complete selected minimal graphs.
Use D=4+2e and the S283 raw i/(16pi^2) master convention. All squared
masses are positive: external mu, heavy n>4mu. Channel a=s,t,u.

## Entire off-shell pair numerator

For external squared momenta ai,aj, Lagrangian masses mi,mj,
z=qi.qj, di=(qi-k)^2-mi, dj=(qj+k)^2-mj, d0=k^2,
delta_i=ai-mi, delta_j=aj-mj, direct full stress contraction gives

N=4z^2-4mi*mj/(D-2)+2z(delta_i+delta_j)
 +(2z+delta_i+delta_j)d0-(2z+delta_j)di-(2z+delta_i)dj.

For the light-light pair ai=aj=mu,z=(a-2mu)/2:
N=V_D(a)+(a-2mu)(d0-di-dj).
For the H-light pair ai=a,mi=n,aj=mj=mu,z=-a/2:
N=n[a-4mu/(D-2)]-n*d0+a*dH+n*dPhi.
The full metric-cubic contact trace is
delta-4m/(D-2)-d0+di, hence its graph contributes
[4m/(D-2)-delta]B0(ai;0,m)+A0(m).
The S292 on-shell formulas follow at a=n; off-shell terms cannot be discarded.

Define:
L(a)=C0 with(0,mu,mu), external(mu,mu,a), the whole light cusp;
U(a)=C0 with(0,n,mu), external(a,mu,mu), the whole H-light cusp;
T(a)=C0 with(n,mu,mu), external(mu,mu,a), the active-light matter triangle;
Bll(a)=B0(a;mu,mu), Bh(a)=B0(a;0,n),
Blon=B0(mu;0,mu), Bmix=B0(mu;n,mu).
Then the complete proper amputated H Phi Phi correction is
delta_vertex/g=P(a)/(16pi^2 kappa), with

P(a)=-V_D(a)L(a)-2n[a-4mu/(D-2)]U(a)
 -(a-2mu)Bll(a)+2n Bmix
 -4mu(D-4)/(D-2)Blon
 +[-a-n+4n/(D-2)]Bh(a)+2A0(mu)+A0(n).

External four-Phi residues contribute -2Sigma'_Phi times the full
matter heavy tree. No external H sqrtZ is inserted on an internal line.
The entire pure-GR heavy selfenergy is
Sigma_H(a)=[(4na-4n^2/(D-2))Bh(a)-2nA0(n)]/(16pi^2 kappa).
Expansion of g^2/(n-a-Sigma_H) gives POSITIVE g^2 Sigma_H/(n-a)^2.

There is additionally the two metric-cubic contact bubble:
-g^2*2D/(D-2)*Bh(a)/(16pi^2 kappa), each crossed channel.
Literal phase: (ig)^2*i^2 * [i B0/(16pi^2)] = +i g^2 B0/(16pi^2);
eta.P.eta=-2D/(D-2). This bubble is required, not optional.

With S'_mu=16pi^2 kappa Sigma'_Phi and S_H=16pi^2 kappa Sigma_H,
the full nonendpoint/nonbox bracket per channel is
[2P(a)-2S'_mu]/(n-a)+S_H(a)/(n-a)^2-2D/(D-2)Bh(a).

## Mixed boxes and exact denominator cancellation

Route d1=k^2,d2=(p-k)^2-mu,d3=(k+q)^2-n,d4=(r+k)^2-mu,
p^2=r^2=(p+q)^2=(r-q)^2=mu, a=q^2, b=(p+r)^2.
The entire numerator is V_D(b)+(b-2mu)(d1-d2-d4), not a scalar-only guess.
Thus the literal box is
V_D(b)D0(a,b)+(b-2mu)[T(b)-2U(a)].
Twelve labelled graphs reduce to twice the six ordered(a,b), a!=b.
The explicit enumeration below proves the factor two; the full cut inventory supplies an independent normalization check.
The box common prefactor is +g^2/(16pi^2 kappa).

Sum of denominator-cancelled terms is exactly
4 sum_a[(a-2mu)T(a)+a U(a)].
Together with P this gives the U coefficient
-4[a^2-4nmu/(D-2)]/(n-a).

The full minimal g^2/kappa amplitude is the sum of these graph groups
PLUS the entire S290 matter-graviton endpoint (quartic C set0),
including both active-light and active-heavy triangles, H-metric
bubble, all metric contacts, covariant OS terms, and unfixed R H anchor.
Its constant R Phi^2 anchor has zero single-insertion b20 but is not
a higher-EFT matching prescription. Flat onepoints use the original
whole density condition. Quadratic metric tadpoles are scaleless.
No pure-H self-interaction or extra H^2 Phi^2 contact is introduced.

## Separated UV and IR checks

Separated UV coefficients (not combined raw on-shell Laurent coefficients):
P_UV=(2a-4mu-4n)/e,
S'_mu,UV=-4mu/e,
S_H,UV=-4n(a-n)/e.
These yield -4/e per channel before the contact-contact bubble; that bubble
has +4/e and cancels it. Cusp IR poles must never be counted as UV.
The only remaining selected-sector UV issue is the S290 covariant R H
metric mixing counterterm; its finite matching coordinate remains explicit.

For a!=n, U(a) has only one on-shell massive line and is IR integrable.
L(a) is the usual two-on-shell light cusp. The mixed box soft limit is
D0(a,b)=Gamma(1-e)(4pi nu^2)^(-e)*J_b(e)/[2e(n-a)]+finite.
Its six ordered pieces reproduce the cross-channel F_b part of S278;
the proper light cusp plus four external residues supply the remaining
2F_a-mu. Therefore the whole sector factors the full heavy tree,
not only the resonant term.

## Positive Feynman box domain below threshold

Let u=x2+x4,v=x2/u,h=x3,x1=1-u-h. Jacobian u,
0<u<1,0<h<1-u,0<v<1. The entire denominator is
Delta=u^2[mu-b v(1-v)]+h[n-a(1-u-h)].
At subthreshold a,b near0 or2mu and n>4mu, Delta is positive.
At a=0 the h integral is elementary and gives exactly
D0(0,b)=Gamma(1-e)(4pi nu^2)^(-e)/n *
[J_b(e)/(2e)-integral_0^1du dv u*(u^2 A_b+n(1-u))^(e-1)].
The general a/b subtraction and explicit derivative bound are proved in notes/forward.md.
Fixed-angle superficial degree0 is NOT a fixed-transfer Regge proof.

## All twelve labelled boxes

Fix external label0 and enumerate the three unoriented cycles of four
distinct labels. In each cycle choose one of four edges for H; the
opposite edge is h and the remaining two are Phi. This gives12 graphs,
with no interchange of differently colored H/h edges. The invariant b
is the pair at the h edge and a is the neighboring Phi-edge pair.
The two complementary assignments give the same ordered(a,b). Hence
each of(s,t),(s,u),(t,s),(t,u),(u,s),(u,t) occurs exactly twice.

Independent full symmetric-tensor contractions in D4,5,6 retain all
loop components and offshell mass defects. They agree with the general-D
Gram expression before any propagator is canceled. The two cubic-contact
trace is eta.P.eta=-2D/(D-2). It is this full trace, not a chosen helicity,
that fixes the additional bubble and its UV cancellation.

## Selected graph UV and matching boundary

The flat nonendpoint g^2/kappa sum is UV finite only after combining the
whole offshell terms and contact-contact bubble. This does not renormalize
away a cusp IR pole. The remaining selected UV structure in the entire
matter endpoint is the covariant R H mixing from S6.290. Its complete
D-dimensional tensor is subtracted before its finite limit. Its finite
coefficient is independent matching data, not fixed by the flat UV cancellation.

All displayed raw integrals carry the common(4pi nu^2)^(-epsilon) factor.
The T and U API kernels used by the finite formula are evaluated at
epsilon=0, where this factor is1. The full-D tests retain it explicitly.
