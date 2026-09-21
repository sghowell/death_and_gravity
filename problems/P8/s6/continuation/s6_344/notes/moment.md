# Exact finite heavy-mass moment and original sign

The primitive coefficient in light-mass units is

    c(n)=-4/15 integral_0^1 z^2*(1-z)^5/[1+(n-1)z]^4 dz.

It is strictly negative for every n>0. For n!=1 its exact value is

    -[4n^5-60n^4 log(n)+155n^4-240n^3 log(n)+80n^3
      -120n^2 log(n)-220n^2-20n+1] / [45(n-1)^8].

To verify the complete integral without numerical fitting, put
y=1+(n-1)z. The integrand becomes
(y-1)^2*(n-y)^5*y^-4/(n-1)^8. Expand only this finite polynomial,
integrate each integer power (using log(y) at power-1), differentiate
the full primitive, and evaluate both exact endpoints1,n.
The construction also works with the reversed orientation when0<n<1.

At n=1 the integral is-4/15*Beta(3,6)=-1/630; the rational-log expression
has precisely this removable limit. Also lim n^3*c(n)=-4/45.
Exact rational diagnostic masses1/2,2,3,10 check the endpoint formula.
For stable numerical evaluation near n=1 the positive integral or a
local expansion is preferable to subtracting large rational-log terms;
floating values are not used as proof.

For n>1, discard(1-z)^5, substitute y=(n-1)z and extend the positive
integral to infinity. Since integral_0^infinity y^2/(1+y)^4=1/3,

    0<-c(n)<4/[45(n-1)^3].

For the original n>10^6 this is<n^-3.
The original endpoint is
C+g^2/n=-2*g^2*(n^2-2n-2)/[n*(n-2)^2]<0.
Its magnitude is<3g^2/n because the difference between3 and its
dimensionless ratio is(n-4)^2/(n-2)^2>0.
Thus the dressed known coefficient is positive and

    |chi_triangle|<g^4/(48*n^4),
    |chi_triangle|/A0<g^2/(192*n)<10^-206,

using16*pi^2>144 and the already established original positive
Born A0>4g^2/n^3. These are coefficient bounds, not a bound on the
discarded physical Taylor remainder or on any unknown additional chi.
The mass dimension is consistent: c carries dimension-8 before
setting the light mass to1; g^2*c has dimension-6.

