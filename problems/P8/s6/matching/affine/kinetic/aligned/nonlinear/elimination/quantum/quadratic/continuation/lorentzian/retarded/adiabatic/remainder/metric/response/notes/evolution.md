# Exact prepared transport, physical readouts and the fixed scalar action

Let f=(2W8)^(-1/2)exp(-i Theta), Theta'=W8. Write the exact
selected mode and its derivative using coefficients A_mix,B_mix.
The unchanged transport matrix preserves
|A_mix|^2-|B_mix|^2=1. The source is zero on an initial
neighborhood, so delta A_mix(u0)=delta B_mix(u0)=0.

The unvaried initial B_mix is not zero. The frozen preparation
and evolution bounds give

    |B_mix(u,p)|<=C_mix/nu^6<1,
    C_mix=1813229+5347035781757616/m^4,
    nu^2=m^2+|p_com|^2/Amax^2, Amax=25/16.

In particular |A_mix|<2. Both source variations must retain
the resulting interference and its phase response.

## Variation of transport

On I, nu<=omega<=Amax nu. Let C be the old residual bound,
C_delta the new varied residual bound, and B_freq the bound
for |delta W|/omega per unit source norm. Then

    |delta Theta|<=B_freq Amax nu,
    ||K||_infinity<=2C/nu^9,
    ||delta K||_infinity
       <=(2C_delta+4C B_freq)/nu^9+2C B_freq Amax/nu^8.

The interval has length one and its transport exponential is
below two. Variation of constants gives, with retained slack,

    |delta A_mix|,|delta B_mix|<=D9/nu^9+D8/nu^8,
    D9=ceil[16(C_delta+2C B_freq)],
    D8=ceil[16C B_freq Amax].

Common constants are D9=10593275744081548 and
D8=432407255532885.

## Differentiated physical readout

The frozen reference bounds imply |f|^2<=1/omega and
|p_ref|<=2omega|f|. Since |A_weight|<=1 and |B_weight|<=2,
both the real reference readout J_ref and pair amplitude Z_ref
are bounded by 3omega.

For a unit joint source norm set

    E=B_freq+C_delta_c/m_min,
    T_f=B_freq(1+Amax nu),
    C0=2C_delta_A+2E+C_delta_B/2+2C_r.

Direct differentiation gives

    |delta f|<=T_f |f|,
    |delta p_ref|<=omega |f|(E+2T_f),
    |delta J_ref|,|delta Z_ref|<=omega(C0+6T_f).

Therefore a valid nu^2 readout constant is

    K_ref=ceil{Amax[(C0+6B_freq)/m_min+6B_freq Amax]}.

The actual sector ceilings are 49 and 50 for both outputs.
This derivation uses the newly computed contact bounds; it
does not assume the old mass-only bounds of two for them.

Differentiate the exact mixing identity, using its Wronskian:

    |delta(J_exact-J_ref)|
       <=[30 Amax nu(D9/nu^9+D8/nu^8)
          +6 C_mix K_ref/nu^4] ||(n,zeta)||_C10.

The last term is required even with zero independent initial
response. The same squeezed initial phase control is -15/4.

The comoving integrals of nu^-4,nu^-7,nu^-8 are
Amax^3/(8pi m), Amax^3/(15pi^2 m^4), Amax^3/(64pi m^5).
Use pi>3, sum all three polarizations, and use a^-3<=1.
The varied-reference tail contributes at most 3T/(54m^2).
After division by L^2, the joint nonlocal bound is below
10^-42 at L=10^24,m=1000. Add the S6.67 local component,
using C10 control of C4, to obtain the bound below 10^-38.

## The already-fixed tadpole action

On an open neighborhood of the clock, the S6.60 scalar cutoff
is identically one and the original added scalar density is

    P_delta=-p_sigma+(rho_sigma+p_sigma)(x+1)/2,
    x=-N^-2.

The c-number functions rho_sigma,p_sigma were computed once
from the selected original state and are held fixed here.
The relative coordinate density is N exp(3 zeta)P_delta.
Its gradient is (rho_sigma,-3p_sigma), and its Hessian is

    [ -(rho_sigma+p_sigma), 3rho_sigma ]
    [ 3rho_sigma,           -9p_sigma ].

Varying the physical output normalization gives

    delta rho_delta=(rho_sigma+p_sigma)n,
    delta p_delta=(rho_sigma+p_sigma)n.

Both the density Hessian and this normalized Jacobian are
checked independently from the existing S6.60 action.
The frozen S6.55 bounds on the complete selected profiles
control the additional C0-to-C0 term by
|rho_sigma|+|p_sigma|. Combining it with the vector response
keeps the full background-cancelled Gaussian correction below
10^-37 at the stated scales. No profile is reselected and no
new cancellation action is introduced.
