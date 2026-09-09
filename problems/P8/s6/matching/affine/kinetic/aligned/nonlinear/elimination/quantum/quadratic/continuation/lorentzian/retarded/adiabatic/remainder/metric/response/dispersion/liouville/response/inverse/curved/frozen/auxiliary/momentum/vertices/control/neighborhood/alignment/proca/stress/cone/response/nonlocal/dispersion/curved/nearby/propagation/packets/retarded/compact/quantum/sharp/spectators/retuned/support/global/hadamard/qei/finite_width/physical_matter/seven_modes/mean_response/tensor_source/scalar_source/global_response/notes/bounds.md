# Two-sided global actual mean bounds

S6.104 proves J>=1/(5t), with a strictly positive coefficient
remainder, and a weighted mean-generator integral below 2
after removing only nonpositive momentum damping. Its
Gronwall factor is therefore below 9 on either half-line.

Let I_f be the sum of the first two compensated forcing
half entry integrals from notes/compensation.md. With
C_cov=(40*3^59)^2 and zero anchor Y=(xi,100t^3 d),

    ||Y||_infinity <= C_Y eta_S,
    C_Y=9 I_f C_cov <10^64.

Each lapse state kernel also has a uniform bound. Let
C_nstate be half the sum of those entry bounds, using
|u|/t^s<=1/2 for s>=1 and 1/t^s<=1 for s>=0.
The classical constraint contribution is bounded by
(3/10)C_Y |u|/t^3+(3/40)C_Y/t^11.
Consequently

    |n|<=C_N eta_S,
    C_N=(9/40)C_Y+C_nstate C_cov <10^64,
    |xi+n/(2h)|<C_B eta_S, C_B=C_Y+C_N/2<10^64.

The actual spatial-frame logarithmic scale, rather than
only its linearization, is bounded by
(C_Y+C_N)eta_S<2*10^64 eta_S on the chosen amplitude domain.

Using |beta|<=1/(10t^6), ell=1/(10t^6), and the fourth
actual forcing integral,

    |delta psi|<=C_psi eta_S,
    C_psi=(C_N/10+3C_Y/10)I_6+I_psi C_cov<10^64.

All rounded comparisons are evaluated as exact rational
inequalities. No floating exponential is used. The selected
global amplitude is 0<=eta_S<=10^-80, a subset of the
previous local family. Thus its lapse and scale responses
are uniformly tiny, but not zero.

Every nonzero lapse state entry has strictly decaying
effective power. Hence n tends to zero. The scale derivative
is integrable because the first row of the weighted mean
generator and both forcing kernels have finite absolute
integral. The matter-field derivative is also integrable.
Thus xi and delta psi have finite limits at both ends.
Reflection gives equal xi limits and opposite field limits.

Restore the ORIGINAL trace:

    |delta_p|/eta_S <= C_Y/t^3+36|u|Q^2/t, Q=3^59.

The second term uses V<=Q^2 eta_S, not an assumption that
V vanishes. This tends to zero but is only an O(1/|u|)
bound on the physical trace. The stronger t^-3 bound
belongs to d, not to delta_p. No stronger tail is claimed.
