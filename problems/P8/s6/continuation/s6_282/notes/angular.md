# Complete angular integral

For every monomial, divide x^a/(h-x^2) and y^b/(h-y^2) independently.
The product decomposes exactly into a double resolvent, two single
resolvents and a polynomial. Summing the complete original N, rather than
fitting a finite set of angles, gives C00 J0+C11 J1+CL L0+CP:

    C00=2[h^4+4h^3+20h^2z^2+2h^2+8hz^2-4h+4z^4-4z^2+1],
    C11=-16z(h^2+2h+2z^2-1),
    CL=-3h^3z^2-h^3+15h^2z^2-19h^2-45hz^2+h-15z^2+3,
    CP=(45h^2z^2-210hz^2+100h+87z^2-14)/15.

Here J0=<1/[(h-x^2)(h-y^2)]>, J1=<xy/[(h-x^2)(h-y^2)]>,
L0=<1/(h-x^2)>. The normalized moments are <x^2>=1/3,
<xy>=z/3 and <x^2y^2>=(1+2z^2)/15. Conditional azimuthal averaging gives

    <y^2/(h-x^2)>=P2(z)(hL0-1)+(1-z^2)L0/2,
    <xy/(h-x^2)>=z(hL0-1).

The polynomial identity, every coefficient and these moments are checked
independently. No printed external moment formula is assumed.

Partial fractions and a Feynman parameter give
J0=[I(h,z)+I(h,-z)]/(2h), J1=[I(h,z)-I(h,-z)]/2, where

    I(h,z)=integral_0^1 da/[h-a^2-(1-a)^2-2za(1-a)].

The spherical squared-linear-denominator primitive establishes this
identity for physical h>1. Its analytic extension then follows on the
domains stated in sheets.md. For h>1,

    I=2atanh sqrt[(1-z)/(2h-1-z)]/sqrt[(1-z)(2h-1-z)].

For h=-tau^2<0, evaluate the negative real denominator directly:

    I=-2atan sqrt[(1-z)/(1+z+2tau^2)]/sqrt[(1-z)(1+z+2tau^2)].

The z=1 values are1/(h-1) and-1/(1+tau^2). The single resolvents are
atanh(1/sqrt(h))/sqrt(h) and-atan(1/tau)/tau, respectively.
Blind substitution into a principal-log expression can select the wrong
branch and is not the proof. At z=1, J0=1/[2h(h-1)]+L0/(2h),
J1=hJ0-L0. The full result exactly reduces to the accepted S280 F(beta).
