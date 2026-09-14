# Explicit cutoffs, complete ordering contacts and positive calibration

## Exactly two derivative-controlled regulators

Let b(y)=exp(-1/y) for y>0 and0 otherwise. Put
theta(y)=b(1-y)/(b(y)+b(1-y)) on0<y<1, extended as1 to y<=0
and0 to y>=1. Its full smoothness at the endpoints follows from
the flatness of exp(-1/y); the denominator stays positive there.

For x=1/y the complete nth derivative of b is exp(-x) P_n(x), with
P_0=1 and P_(n+1)=x^2(P_n-P_n'). For x>=0 the full positive
exponential series gives x^k exp(-x)<=k!. Summing absolute polynomial
coefficients times k! gives the exact outward ceilings
B=[1,2,36,1584,129600] through order4, including the constant
extensions outside the transition.

At least one of y,1-y is>=1/2, so the denominator is>=exp(-2)>1/9.
The latter follows from e<3. Differentiate denominator*theta=b(1-y)
and solve for the highest derivative, retaining all lower products.
With C0=1 the complete recurrence is
C_n=9[B_n+sum_(k=1)^n binom(n,k)*2 B_k C_(n-k)].
It gives C=[1,54,4860,672624,125411328].

Only c=1 and c=2 are claimed. For
chi_c(w)=theta((||w||^2-R^2)/(c R^2)), R=10^20,
the support radii are sqrt(2)R and sqrt(3)R<2R. Both cutoffs equal1
through the core R. Translation invariance of the actual seed makes
the original translation representation orthogonal after whitening;
therefore these radial cutoffs preserve the full three translations.

The inner radial map has first derivative norm<=4/R and second
derivative norm<=2/R^2 on the support; every higher derivative is0.
The complete partition rule gives bounds for each mixed derivative
of order0..4, in units R^-n:
K=[1,216,77868,43164576,32234502096].
For example the fourth derivative contains the fourth transition
derivative, all six first-first-second terms and all three
second-second terms. They produce256 C4+192 C3+12 C2.
The code checks both the full one-coordinate chain and an independent
mixed-coordinate fixture. Keeping only the leading fourth derivative
would be incorrect.

## Full holomorphic amplitude, then the real product rule

From notes/domain.md the full amplitude g has complex bound
H=10^1000 and F-1 has bound E=10^-255. The actual centers g0,F0
are retained, so g-g0 and F-F0 have modulus<=2H and2E.
For real phase derivatives through4 the complete Leibniz/Cauchy bound
for chi times either centered amplitude is

J_n=2 A R^-n sum_(k=0)^n binom(n,k) K_k (n-k)!4^(n-k).

Here A=H or E. The exact sums are
L=[1,220,79628,44120112,32940423888].
No omitted cutoff derivative, field Hessian term, cross covariance
or implicit N_i z^i_AB term is permitted by this estimate.
It bounds the entire amplitude, not its low-order Taylor truncation.

## Same coherent map and first Weyl calibration

The S268 normalized map Q_V(a)=W* M_a W uses the SAME full pure seed
and the isometric coherent transform W. It is unital, positive,
contractive on bounded symbols, and obeys Q(a)*Q(a)<=Q(|a|^2).
The normalization is (2pi)^(-d) at unit CCR, d=48; no displaced
coherent vector is a newly selected initial vacuum.

For full covariance V0=S0 S0^T/2, D=(1/2)V0_AB partial_zA partial_zB
becomes Delta_w/4 exactly, retaining all physical cross covariance.
Its coefficients must not be replaced by a guessed diagonal bound.
Thus, with phase dimension96,
||D F_ext||<=96 J2(F)/4 and
||D^2 F_ext||<=96^2 J4(F)/16.
The same inequalities apply to g_ext, with the retained scalar center.

Set A_c(u)=Q_V[(1-D)g_ext]. Its actual norm bound is
3H+96 J2(g)/4<10^1010=:B.
The displayed factor3H is conservative: g0 and the centered amplitude
are each bounded separately, so no background phase is discarded.
The real symbols give bounded self-adjoint operators. Norm continuity
in real u follows from the full smooth source on a common compact
support and the evaluated uniform domain.

For F_ext=F0+chi(F-F0), the real convex combination gives
||F_ext-1||<=E. The complete derivative bound gives
||D F_ext||<10^-275. Hence the real calibrated symbol
(1-D)F_ext is>1-E-||D F_ext||>1/2.
Positive unital Q_V proves that the CALIBRATED volume operator is
strictly positive and has norm distance from identity<2*10^-255.
The uncalibrated positive Q_V(F_ext) is distinct; it too remains
positive and close to identity. Positivity is not inferred merely
from Weyl quantization of a positive classical function.

## Precise ordering boundary

For bounded smooth symbols with the proved fourth derivative bound,
the exact heat identity is

exp(D)(1-D)f-f=-integral_0^1 t exp(tD) D^2 f dt.

The heat semigroup is sup-norm contractive, so the SYMBOL remainder
is<=(1/2)||D^2 f||. For the complete cut-off physical volume this
is<10^-310. The full interaction remainder is also finite, but is
not claimed instantaneously small from the short evolution time.

Since Q_V(a)=Op_W(exp(D)a), the expression identifies a Weyl SYMBOL
difference. It does not automatically bound the norm of the Weyl
operator with that symbol. No original Weyl-operator matching or
physical counterterm theorem follows. The two declared calibrated
operators and their dynamics are well defined independently of that
unproved identification.

The evaluated derivatives are for the two explicit chi_c only.
Arbitrarily oscillatory smooth cutoffs sharing the same support can
have arbitrarily large derivatives; no such wider class is included.

