# Current-only absolute domain and uniform proper-face hierarchy

Fix the real directions, TT fields, E,u, Born A0 and transfers. In
S=S_c(sigma_ab), a null contribution is
 a*(n_a.eps_c.n_a)/(n_a.n_c),
and similarly for b. For unit complex spatial TT eps_c, the numerator
has modulus at most sin(theta)^2. Dividing by1-cos(theta) bounds
the angular coefficient by1+cos(theta)<=2. This uses transversality,
not an absolute bound on separate singular Einstein diagrams. After
this cancellation the null term is linear in its energy.

Only the massive recoil and this kinematic current are continued on
 |a'-a|,|b'-b|<=d, d=10^-6,
about every real positive a,b with a+b<=1/8. The current contains no
hard propagator and no pure-soft current denominator. In the exact
recoil h=E-(a+b)/2, v=a*n_a_space+b*n_b_space, and
Eprime^2=h^2-v.v/4. The real radicands satisfy Eprime^2>=45/32 and
rprime^2>=13/32. Directly,
 |dh|<=d, ||dv||<=2d,
 |d(Eprime^2)|<=4d+dW+2d^2<5d.
Nonzero-disc square-root continuation and factorization give
 |dEprime|<5d, |drprime|<10d,
 |Eprime|>1, |rprime|>1/2, |h+Eprime|>2.

The S313 recoil calculation applies with its cW replaced by d:
for t=rprime/Eprime, |dt|<20d; for f=t/[4(h+Eprime)], |f|<1/8
and |df|<3d. Thus the recoil coefficient beta=f*(u.v)-1/2 has
|dbeta|<d and |beta|<1. The continued outgoing momenta obey
 |dp0|<5d, ||dp_space||<13d.
The real spatial norm is<=sqrt(3), so the continued norm is<7/4.
On the convex future-massive spatial tube |p0|>7/8 and the energy
gradient is<2. The real Doppler gap exceeds1/4; its perturbation
is below65d<1/8. Hence |p.n_c|>1/8. Signed incoming momenta obey
the same magnitude bounds by p->-p.

Every massive term has numerator modulus<4 and denominator>1/8,
so is<32. The null terms together are<=2*(W+2d), whence |S|<129.
Half-radius Cauchy gives
 |S_a|,|S_b|<=2*129/d<10^9,
 |S_ab|<=4*129/d^2<10^15.
No comparable absolute-radius domain is asserted for G2 or G3.
In particular these estimates do not cross a hard-amplitude cut.

On its much smaller radius eW, e=10^-12, S313 proves
 |G2-G20|<B W, B=(10^40*n^2+10^60)/kappa<10^-350.
Because G20 is constant in a,b, Cauchy gives |G2_a|,|G2_b|<=B/e.
The proven mixed derivative is below C2/W, C2=10^-326.
S319 independently supplies the real full-amplitude bound
 |G2|<=D2=6*n^2*(2*10^16)^2/kappa<10^-371.

Apply the ordinary product rule to the complete c face S*G2/sqrt(kappa).
All four mixed terms are retained. It follows that its first
derivatives are bounded by
 Lderived=(129*B/e+10^9*D2)/sqrt(kappa)<L=10^-735,
and its mixed derivative is bounded by Q/(a+b), where
 Qderived=(129*C2+2*10^9*B/e+10^15*D2)/sqrt(kappa)
          <Q=10^-723.
In the last bound W<=1 absorbs the constant product-rule terms
into1/W. The arithmetic uses only original parameters.

The uniform first derivatives give Lipschitz extensions on axes,
and the compatible corner limit gives |g_i|<L*w_i. Integrate the
mixed derivative over trimmed positive rectangles and then take
the compatible limits. Since
 int_0^a int_0^b 1/(u+v) dv du=I(a,b),
one obtains |g_ij|<Q*I(w_i,w_j) for every proper pair. This uses
local integrability and the scaled amplitude extension, not a
derivative evaluated at an undefined bare zero-energy pole.
