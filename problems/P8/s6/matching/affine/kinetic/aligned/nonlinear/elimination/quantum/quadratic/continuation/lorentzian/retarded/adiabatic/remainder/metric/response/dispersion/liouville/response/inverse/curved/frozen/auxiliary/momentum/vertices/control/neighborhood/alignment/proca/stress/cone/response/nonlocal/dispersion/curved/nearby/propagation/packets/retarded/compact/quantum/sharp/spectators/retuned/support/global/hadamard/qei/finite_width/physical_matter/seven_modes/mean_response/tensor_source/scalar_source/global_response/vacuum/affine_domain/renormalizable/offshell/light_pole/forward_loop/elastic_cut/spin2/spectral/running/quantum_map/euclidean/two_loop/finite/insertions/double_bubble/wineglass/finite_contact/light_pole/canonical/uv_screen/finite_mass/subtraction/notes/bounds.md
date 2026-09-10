# A coefficient bound from a genuinely complex error domain

Let h=delta rho. The domain proof gives h holomorphic
on the radius-five disc about two with |h|<=M_rho.
On the radius-four disc, Cauchy's derivative estimate
with margin one gives |h'|<=M_rho.
The strict routing inequalities leave an open margin
around the closed domain; alternatively the boundary
estimate follows by increasing inner radii to the limit.

For |s-2|<=1 and x in [0,6], the segment between s and
x stays in the radius-four disc. Hence the removable
divided difference obeys

    |(h(x)-h(s))/(x-s)|<=M_rho.

Each integral has length six. The sum of the two
channel quotient integrals is bounded by 12M_rho.
This establishes holomorphic control, not merely
a pointwise principal-value estimate on the real axis.

On the unit disc, 1<=|s|<=3 and 3<=|6-s|<=5.
The real logarithmic parts have magnitude below two,
since exp(2)>5. Their arguments have magnitude below
one: the discs remain in the positive real half-plane.
Thus each of Log(s) and Log(6-s) has magnitude below
three. Together with pi<4, either physical local
logarithm difference has magnitude below ten.

Both logarithm products contribute at most 20M_rho.
Including the quotient terms and the 1/pi prefactor,

    |delta F_low^germ(s)|<=32M_rho/pi, |s-2|<=1.

Cauchy's unit-radius coefficient bound now gives

    |delta b2_low|<=32M_rho/pi=288 K eta.

It would be invalid to start only with a small real-axis
spectral error and differentiate its singular Cauchy
integral twice as an ordinary integral. The larger
complex domain and removable representation are the
additional ingredients that justify this estimate.

## Exact numerical enclosure, without floating decisions

At the unchanged m=10^200, epsilon=4*10^-194 and

    eta=10*10^-194+2*10^-388,
    288eta<10^-190.

For log(2), integrate the positive geometric expansion
of 2/(1-x^2) from zero to r=1/3:

    log(2)=2 sum_{j>=0} r^(2j+1)/(2j+1).

After N=8 terms the omitted positive tail is strictly
less than 2r^(2N+1)/[(2N+1)(1-r^2)].
The exact rational partial sum and this rational upper
bound enclose log(2); no decimal logarithm decides a sign.

Combining this enclosure with the finite-mass error
gives a strict normalized band

    -4 < b2_low/K < -3.

The unchanged positive K has the rational upper bound
K<10^-1620 inherited from S6.128. Consequently
|b2_low|<4*10^-1620. This is the coefficient of the
specified finite-cut subtraction alone, not a bound
for the full new-model amplitude or its regular part.
