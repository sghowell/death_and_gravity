# Finite box moment, exact limits and honest bounds

With t=z(1-z), define

    c_box(n)=-4/315 integral_0^1
       t*(18+18t+92t^2+363t^3)/[1+(n-1)z]^5 dz.

The integrand proves c_box(n)<0 for every n>0.
Put y=1+(n-1)z and expand the finite numerator polynomial.
The transformed primitive is the sum of its integer powers divided
by y^5, with log(y) only at power-1. Differentiate the whole primitive
and evaluate both exact endpoints1,n. This gives, for n!=1,

    -2/[315*n^3*(n-1)^9] *
    [3n^10-21n^9+91n^8+210n^7 log(n)-1859n^7
     +10920n^6 log(n)-21896n^6+28560n^5 log(n)
     +10920n^4 log(n)+21896n^4+210n^3 log(n)
     +1859n^3-91n^2+21n-3].

The exact equal-mass value is-58/945, both by the Beta moments of
t,t^2,t^3,t^4 and by the removable limit of this expression.
Also lim n^2*c_box(n)=-2/105.
Reflection z->1-z proves c_box(n)=n^-5*c_box(1/n), independently
checked in the closed form. Exact diagnostic masses1/2,2,3,10
check the complete endpoint formula without floating arithmetic.

For0<=t<=1/4,
18+18t+92t^2+363t^3<=2171/64.
Use t<=z, substitute y=(n-1)z and extend the positive integral
to infinity. Since integral y/(1+y)^5=1/12,

    0<-c_box(n)<2171/[60480*(n-1)^2], n>1.

For the original n>10^6 this is<1/(25n^2), so

    |chi_box|<g^4/(3600n^2),
    |chi_box|/A0<g^2*n/14400<10^187,

using16*pi^2>144 and the original A0>4g^2/n^3.
The second bound is LARGE. It is neither perturbative smallness
nor a physical scattering estimate. The local comparison can
redistribute sizable terms between the covariant flat representative
and its explicit curvature correction. An aggregate must use
a consistent convention and include possible cancellations.

The integral is the stable sign representation; the exact rational-log
formula contains cancellations near n=1. No numerical evaluation is
used as proof. Restoring the light-mass units, c_box has mass
dimension-10 and g^4*c_box dimension-6 as required.

