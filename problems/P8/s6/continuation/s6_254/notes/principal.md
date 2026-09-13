# Complete fast symbol and uniform fixed-background remainder

In the central chart the full lapse square contains

    -2a Theta P^2 Qb/D -(L2/2+L0 a^2/(2P^2)) Pb
    -w.ps/Z -3sqrt(zeta)(r Theta/D+R_N dH)PL
    +3a^3 Theta c.sigma/D-a^3 d.sigma.

For L2 nonzero, the exact old-from-new canonical shear

    Pb_old=Pb_new-4a Theta P^2 Qb/[D(L2+L0 a^2/P^2)]

removes its Qb term. If old_phase=R new_phase, the new Hessian is
R^T H R+Omega R^-1 Rdot. The entire Rdot, including a_dot=a Hhat and
Theta_dot,D_dot,L0_dot,L2_dot, is included. Time derivatives of other
coefficients are not set to zero: they enter the smooth coefficient
matrix and later growth-chart derivative, but not this algebraic map.

Now use phase similarity x=M y with

    M=diag(1,1,1,P^1/2,P^3/2,P^1/2,P^1/2,P).

M depends on the fixed comoving momentum, not time. In order
(Qb,sigma1,sigma2,QL,Pb,ps1,ps2,PL), the full normalized generator tends to

    Qb'=alpha Pb,            sigma_i'=beta_i Pb,
    QL'=-g Qb,               Pb'=g PL,
    ps_i'=-charge_i Qb-d sigma_i,  PL'=-h QL,

where these primes denote the principal matrix action, not yet physical
time derivatives, and

    alpha=L2^2/(8J a^3), beta_i=L2 w_i/(4J Z a^3),
    g=r sqrt(zeta)/(D a^2), h=a^3 Yv/zeta,
    charge_i=a c_i/D, d=aY.

The complete8x8 matrix has characteristic lambda^4(lambda^4-Bfast),
Bfast=alpha g^2 h=r^2 L2^2 Yv/(8D^2 J a^4)>0. The four other physical
modes, two TT and two transverse Proca, have generator O(P) after their
ordinary momentum scaling and therefore zero P^3/2 symbol. Keeping them
all gives the complete16-phase characteristic lambda^12(lambda^4-Bfast).
At r=0 this fast characteristic is lambda^16, not a reference instability.
The slower reference cones remain the S251/S253 cones.

The stored exact remainder is obtained from the ENTIRE rational matrix,
not a discarded-term prescription. Put e=P^-1/2. Every entry of the
normalized generator minus the listed fast matrix is e times a rational
function of e with a polynomial numerator and denominator. At e=0 every
denominator factors into nonzero constants and powers of a,D,J,K,Z,L2.
All64 numerator/denominator degrees and denominator values are stored and
their domain factors checked. On a compact time interval with strict
margins, these functions and their first time derivatives are uniformly
bounded for e sufficiently small. This proves a C1 O(P^-1/2) normalized
remainder. Smooth full coefficient bindings, including the giant fixed
heavy mass and full source, change its finite constant, not this order.

For comparison only, a separate eight-phase congruence with weights
(P^2,P,P,1,1,1,1,P) yields the full raw determinant coefficient at P^10

    -r^2 L2^2 Y^2 Yv/(8 D^2 J Z^2 a^8).

An exact rank-one separation proves this whole determinant identity.
Its sign is NOT used as the time-dependent stability argument; the next
note proves growth for the actual varying full generator.
