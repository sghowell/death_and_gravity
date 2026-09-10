# The complete renormalized MS proper kernel

The S6.139 local references are not substituted
for the momentum-dependent self-energy. Retain
the entire one-loop proper inverse correction,
including its mass and kinetic counterterms.

With x on the internal fermion denominator,
at mu=mF and Euclidean external momentum p,

    Delta_s=x mF^2+(1-x)+x(1-x)p^2,
    Delta_g=x mF^2+x(1-x)p^2.

Let Js=log(Delta_s/mF^2) and Jg similarly.
All displayed x integrals run from zero to one.
The scalar MS inverse correction is

    Sigma_s=(Y/Q)[mF integral Js
                   -i slash(p) integral (1-x)Js].

The Feynman-gauge correction is

    Sigma_g=(a Cf/Q){
        mF[-2-4 integral Jg]
        +i slash(p)[-1-2 integral (1-x)Jg]}.

Here Q=16 pi^2. The constants -2 and -1 arise
from the dimensional numerator multiplying
the pole. The code recomputes their finite
coefficients from d=4-2 epsilon and pairs
the actual S6.139 mass/kinetic counterterms.
At p=0 the result exactly reproduces the
frozen eta and z anchors.

This is the full finite proper MS kernel,
not a MOM Taylor approximation. Each proper
counterterm is included once in Sigma_MS,
and the marked outer line is -S Sigma_MS S.
The remaining overall four-Phi divergence
is a local quartic contact; the soft
projection below removes it before the
continuum outer integral is bounded.

The external scalar canonical conversion
does not multiply this new formal order-two
family again at the same order: changing its
leading field/coupling factors affects it at
order three. The distinct order-two field
conversion of the old order-one box remains
a separate matching-ledger contribution.
