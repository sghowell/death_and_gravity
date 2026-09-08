# Physical output normalization and continuous local bound

The computed local Euler rows are normalized by the background
a^3, after their variation. Physical stress normalization must
itself be varied:

    delta rho_local=-K_N-3 zeta rho_local,
    delta p_local=K_Z/3-(n+3 zeta)p_local.

Here rho_local and p_local are exactly the existing S6.53 local
finite coefficients. These additional normalization terms cancel
the undifferentiated constant-spatial-rescaling response in both
physical outputs. They are not an optional change of subtraction.

Every coefficient of the two sources and their derivatives zero
through four is an exact rational function c(u), times its
appropriate m^(4-2j). Its denominator has the form
c0(1+u^2)^k with c0>0 and k a nonnegative integer. If its numerator
is sum_l a_l u^l, then throughout I=[-1/2,1/2],

    |c(u)| <= sum_l |a_l| 2^-l/c0.

The verifier checks this denominator form exactly. Summing these
nonnegative envelopes over both source rows bounds the operator
in the maximum C4 input norm. Divide by 64 pi^2 L^2 and use
pi^2>9. No time sampling enters this proof.

At L=M tau=10^24 and m=m0 tau=1000, the joint physical bound is

    30738297112762532011 /
    5804752896000000000000000000000000000000000000000000000000,

which is strictly below 10^-38. Both the separate component
bounds and their three adiabatic-order pieces are retained.

## The fourth-derivative control is not a spectrum

The finite local Euler-current fourth-derivative matrix is

    [ 1256/(6561 r^6), -76/(243 r^3) ]
    [ -76/(243 r^3),   -4            ], r=1+u^2.

Its determinant is -50992/(59049 r^6). The indefinite sign is
explicitly checked. Neither positivity nor extra ghost particles
are inferred by resumming this finite local derivative truncation.
The physical nonlocal Gaussian response and the interacting EFT
frequency range would have to be controlled before a spectrum or
quantum stability verdict.

This local estimate therefore does not complete the nonlocal
metric block, prove a no-loss coupled inverse, establish quantum
cones or cutoff, or satisfy the common-parent V/G/B obligations.
The original clock, state and fixed tadpole profiles are unchanged.
