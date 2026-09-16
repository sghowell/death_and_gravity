# Whole Gamma-factor bound and original uniform error

Write F(a)=exp(-EulerGamma*a)/Gamma(1+a). Its convergent Weierstrass
product gives
lnF(a)=sum_{j>=1}[ln(1+a/j)-a/j].
For everyy>=0,
0<=y-ln(1+y)<=y^2/2.
Both inequalities follow by differentiating functions vanishing at0:
the relevant derivatives are y/(1+y) and y^2/(1+y).
Since sum1/j^2=pi^2/6, for ALL a>=0,

-pi^2*a^2/12<=lnF(a)<=0,
exp(-pi^2*a^2/12)<=F(a)<=1,
0<=1-F(a)<=pi^2*a^2/12.

These are whole-function bounds, not a truncated small-a guess.
The first nontrivial total-energy correction is
F(a)=1-pi^2*a^2/12+O(a^3).
Its two-real interpretation is the excluded wedge
integral_0^1 dx/x integral_(1-x)^1 dy/y=pi^2/6,
multiplied by the two-boson factor a^2/2.
An INDIVIDUAL-energy cut without a total constraint instead gives
expDelta*x^a and would miss F(a).

The analytic detector factor obeys
|A_E|^2=x^a*|A_hard,nu|^2 in the defined leading-soft sector:
the real part of the unchanged virtual exponent fixes this power.
Thus I_lead(E)/|A_E|^2=expDelta*F(a). At first Newton order this is
1+Delta, exactly S296. The Gamma effect begins at second Newton order.

On the inherited compact domain mu1,hard COM scalar energy[5/4,2],
S296 gives0<=K0<=(s-4)^2/5<=144/5 and|Delta|<112/kappa.
Therefore0<=a<4/(5*kappa) usingpi>3.
At the original kappa=10^800,pi^2<12 yields
1-F<1/kappa^2 and exp|Delta|<3. Taylor's integral remainder for the
exponential then bounds

|expDelta*F-(1+Delta)|
<=exp|Delta|*(1-F)+Delta^2*exp|Delta|/2
<18819/kappa^2<20000/kappa^2=2*10^-1596.

Multiplication by0<x^a<=1 preserves the same ABSOLUTE normalized bound.
No bound on|lnx| is needed after the prescribed regulator removal.
This is a conversion-factor error, not a complete physical cross-section
error: the unknown hard factor and all non-leading radiation remain.
