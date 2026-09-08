# Selected-state mixing and the continuous response bound

For each original physical sector, let
f=(2W)^(-1/2)exp(-i Theta), Theta'=W, and write the exact selected
mode as v=A_mix f+B_mix f*. Variation of constants imposes
v'=A_mix f'+B_mix f*'. The exact transport matrix is

    K=i rho/(2W) [ [-1,-exp(2i Theta)],
                    [exp(-2i Theta),1] ].

It preserves |A_mix|^2-|B_mix|^2=1. The selected initial B_mix is
not set to zero. The frozen S6.59 preparation estimates yield the
global bound

    |B_mix|<1,
    |B_mix(u,p)|<=C_mix/nu^6,
    C_mix=1813229+5347035781757616/m^4,
    nu^2=m^2+|p_com|^2/Amax^2, Amax=25/16.

The middle/high initial bands are bounded by the same sixth-power
envelope: B8/(4m)^2<=B6 and B10/(8m)^4<=B6. This follows from
their actual frozen constants and does not replace the selected
preparation. Since |B_mix|<1, |A_mix|<2.

The source and all its jets vanish on an initial neighborhood.
The exact state and eighth-order reference data therefore have
zero source variation at u0, including both mixing coefficients.
The nonzero unvaried initial mixing remains in the response.

## Varied transport

The reference frequency satisfies nu<=omega<=Amax nu. For a unit
C10 source, the bounds of the preceding note give

    |delta Theta|<=B_freq Amax nu,
    ||K||_infinity<=2C/nu^9,
    ||delta K||_infinity
      <=(2C_delta+4C B_freq)/nu^9
        +2C B_freq Amax/nu^8.

The interval length is one. The frozen residual bound makes the
transport exponential less than two. Variation of constants with
zero initial response and |(A_mix,B_mix)|_infinity<2 then bounds
both mixing variations by

    D9/nu^9+D8/nu^8,

where the implementation keeps an extra factor-of-two allowance:

    D9=ceil(16(C_delta+2C B_freq)),
    D8=ceil(16C B_freq Amax).

Common sector envelopes are

    D9=2503997646582128, D8=51021843772450.

## Readout, including the initial interference

Let J_ref be the real reference readout and Z_ref its quadratic
pair amplitude. The exact identity is

    J_exact=(|A_mix|^2+|B_mix|^2)J_ref
             +2 Re(A_mix conj(B_mix) Z_ref).

The reference bounds |d|<=3, |W'/W|<4 and m>=1000 imply
|f|^2<=1/omega, |p_f|<=2omega|f|, and
|J_ref|,|Z_ref|<=3omega. The differentiated normalization and
phase imply

    |delta J_ref|,|delta Z_ref|<=K_ref nu^2 ||n||_C10,
    K_ref<=10.

This follows directly from delta f/f=-delta W/(2W)-i delta Theta,
p_f=(-d-W'/(2W)-iW)f, the independently varied readout weights,
and the coefficient bounds. The code computes the corresponding
sector ceilings 9 and 10.

Differentiating the exact mixing identity, using its Wronskian
relation, gives per polarization

    |delta(J_exact-J_ref)|
       <=[30 Amax nu(D9/nu^9+D8/nu^8)
           +6 C_mix K_ref/nu^4] ||n||_C10.

The last term retains the initial interference and its varied
phase. It cannot be deleted when the independent initial response
vanishes. A squeezed-oscillator control with A_mix=5/4,B_mix=3/4
has phase-variation fixture -15/4 even with zero transport forcing.

## Momentum integration and normalization

Including the comoving measure, the integrals of nu^-4,nu^-7,nu^-8
are respectively Amax^3/(8pi m), Amax^3/(15pi^2 m^4),
Amax^3/(64pi m^5). The physical output factor a^-3<=1 only
reduces this bound. Using pi>3 gives the exact rational bounds
implemented in evolution.py. Sum all three polarizations.

The varied reference-minus-adiabatic contribution is bounded by
3T/(54m^2) with T=85166609. Divide both finite pieces by L^2.
At L=10^24,m=1000, their sum is strictly below 10^-43.
Adding the S6.65 finite local bound, using ||n||_C4<=||n||_C10,
keeps the complete homogeneous mass-response bound below 10^-39.

The output norm is C0 on I. The high source regularity is explicit;
these estimates are not a no-loss operator norm or a contraction
of the coupled metric equations.
