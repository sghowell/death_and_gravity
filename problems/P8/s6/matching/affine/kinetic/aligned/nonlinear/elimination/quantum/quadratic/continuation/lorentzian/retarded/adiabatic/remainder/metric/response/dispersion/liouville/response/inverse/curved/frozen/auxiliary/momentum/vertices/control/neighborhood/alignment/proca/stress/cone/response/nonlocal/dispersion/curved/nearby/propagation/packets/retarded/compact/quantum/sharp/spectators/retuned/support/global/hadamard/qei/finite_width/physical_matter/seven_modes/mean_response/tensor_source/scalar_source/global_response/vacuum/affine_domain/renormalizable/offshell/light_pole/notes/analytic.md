# Actual analytic kernel and remainder

At s=1 let d1=x M+(1-x)^2, a=x(1-x), b=a/d1.
For 0<=x<=1, d1>0 and b>=0. The exact gap identity is

(1-x)/M - b = (1-x)^3/(M d1) >=0.

Consequently b is bounded by (1-x)/M. For complex |s-1|<=r<M,
Delta=d1[1-(s-1)b] stays on one analytic logarithm branch.
The subtracted logarithm is identified by its derivative and its value
at s=1; no invalid global combination of complex logarithms is used:

Pi_R(s)=g/(16pi^2) integral_0^1 [-log(1-z)-z] dx,
z=(s-1)b.

The convergent Taylor series gives
|-log(1-z)-z| <= |z|^2/[2(1-r/M)].
Using integral_0^1(1-x)^2 dx=1/3 yields

|Pi_R(s)| <= g |s-1|^2/[96pi^2 M^2(1-r/M)].

S6.110 independently proved pi>3, so the rational majorant is obtained by
replacing 96pi^2 with 864. At r=1 its coefficient is less than 10^-405.
The analytic factorization
D_OS(s)=(s-1)[1+Pi_R(s)/(s-1)]
extends across s=1 with second factor equal to one. The strict unit-disc
bound is below one, excluding any additional zero there. This is a proof
for the actual one-loop truncated inverse and includes nonreal s.

Also Pi'(1)=g/(16pi^2) integral b dx is strictly positive and smaller than
g/(32pi^2 M)<g/(288M)<3 times 10^-208. Thus the finite local kinetic
coefficient 1-Pi'(1) is positive. The rational majorant itself exceeds
10^-208; no sharper ceiling is asserted without a separate estimate.

The actual two-particle branch threshold is (mu+1)^2: at that value
Delta=[(mu+1)x-1]^2, with its zero inside the parameter interval.
The conservative radius M is a sufficient estimate, not an exact location
of the nearest singularity. Higher-loop multi-light cuts can occur closer
and are expressly outside this one-loop statement.

The exact-radius API accepts finite rationals only. Parameter analyticity
r<M does not alone imply zero exclusion; it reports the additional
factored-error test independently. In particular a radius extremely close
to M can fail this conservative zero-exclusion estimate.
