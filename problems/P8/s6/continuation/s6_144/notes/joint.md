# Joint nonzero-soft remainder: no light-line Cauchy extrapolation

Let S_k=m^2+k^2, S_l=m^2+l^2 and
R=min(sqrt(S_k),sqrt(S_l))/360>=2. Multiply only the external-soft
momenta in the fermion kernel by a complex parameter. Its shifts have
one-norm at most 18R, so the Dirac resolvent Neumann ratio is below 1/2.
Each propagator norm is at most 2/sqrt(S). For arc counts (a,b)
with a+b=4, the trace is bounded by 64/(S_k^(a/2)S_l^(b/2)).
The full tail of degree at least one at physical parameter one is
bounded by twice this norm divided by R. This tail is UV convergent;
the constant MS subtraction cancels before it is formed.

Keep the physical light lines and W outside that Cauchy operation.
Write x=q^2, S=1+x, q_axis^2<=x. For z=a+ib on |z-2|<=1 their
denominator product is

    P=(S-z/4)^2+z q_axis^2.

Put a=2+t, b^2<=1-t^2, -1<=t<=1. Its real part is at least
(S-(2+t)/4)^2-(1-t^2)/16. Subtracting (S-3/4)^2 gives
(1-t)(8S-6-2t)/16>=0. Finally
(S-3/4)^2-S^2/16=3(S-1)(5S-3)/16>=0.
For z=0 the bound is immediate. Hence |D_+ D_-|<=16/S^2.
This uses the product directly, not a false lower bound on each
complex-shifted denominator. For a heavy shift of norm<18 and
M>=10000, the relative change is at most
18/sqrt(M)+324/M<=531/2500<1/2, so |H|<=2/(M+x).
Thus |W|<=V=L+g/(M-3)+4g/M.

Use S^-2<=x^(-7/4): for x<=1 the right side is at least one;
for x>=1, x^(7/4)<=x^2<=(1+x)^2.
Also 1/min(sqrt(S_k),sqrt(S_l))<=1/sqrt(S_k)+1/sqrt(S_l).
Each resulting joint integral is a massive/massive/massless sunset
with gamma=7/4, alpha+beta=5/2 and alpha,beta among 1/2,1,3/2,2.

Schwinger parameters t,u,v give determinant tu+v(t+u).
Integrate v first using B(gamma,2-gamma); then use t+u=T and t/T=z.
The exact integral is

 (m^2)^(4-alpha-beta-gamma)/Q^2
 *Gamma(2-gamma) Gamma(alpha+beta+gamma-4)
 *Gamma(alpha+gamma-2) Gamma(beta+gamma-2)
 /[Gamma(alpha) Gamma(beta) Gamma(alpha+beta+2gamma-4)].

Every displayed Gamma argument is positive. This checks the chord
diagonal, both one-large-momentum regions and the joint UV region;
it is not a hard-region-only estimate. The mass power is m^(-1/2).
The four numerator arguments lie in [1/4,7/4]. Splitting Euler's
integral at one bounds each Gamma by 4+6=10. The three denominator
values exceed 1/2, by the elementary half-integer recurrences and
pi>3. Thus each constant is below 80000, and the two terms from
the minimum below 160000.

Six boxes, the trace-tail factor 128*360, the light factor 16, and
the three channels with combined assignment weight one give

    E_delta <=2123366400000 N Y^2 V/(Q^2 sqrt(m)).

The jointly integrable majorant justifies the series difference and
integrations. The actual external variables can be continued from
the common Euclidean germ on this disc; no unregulated divergent
complex contour is shifted. Sign-flip/permutation symmetry or analytic
branches with nonzero E,p give the invariant forward function.
A Cauchy coefficient bound on the closed unit disc now bounds b2.
No lower-degree cancellation is needed for this full remainder bound.
