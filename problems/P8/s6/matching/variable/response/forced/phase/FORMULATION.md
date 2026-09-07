# S6.27.PHASE — fixed-source physical response has no delta limit

Reviewed child of the frozen S6.26 actual-source theorem, with a
source-hashed read-only certificate. The original action, physical metric, source and four zero
initial data are unchanged. This is a failure of a specified fixed-source
limit, not an exclusion of every local EFT or of original P8.

Fix `0<r<=1/1000`, `0<=K=(tau*k_com)^2<=4` and a real source profile sigma
independent of delta. It is zero outside `[-r,-r/2]`, obeys `0<=sigma<=1`
almost everywhere, and has deficit

`integral_{-r}^{-r/2} (1-sigma(u)) du <=r/100`.

This includes the rectangular pulse and fixed smooth pulses with compact
support strictly inside that interval. All four physical Cauchy data
vanish at `u=-r`. For `0<delta<=min(10^-9,r^2/100)`, let gamma_g,delta be
the actual full coupled physical response to this conserved external TT
probe. At fixed r and K, the claim is

`limsup_{delta->0+} gamma_g,delta(r) - liminf_{delta->0+} gamma_g,delta(r)
 > 3*r^2/20`.

The response is bounded, but does not converge as the center algebraic
mass grows. Multiplying the profile by any fixed positive amplitude
multiplies this separation by that amplitude; this is a linear-response
statement, not a finite-amplitude backreaction claim.

The proof keeps the actual g observable and source weights. It derives an
oscillatory term with amplitude greater than `2*r^2/25`, and allows a
uniform `1200*r^4` error for the bounded pole remainder and full coupled
feedback. Two phase sequences give the stated separation. No numerical
Gamma evaluation, sampled delta oscillations or isolated impulse diagonal
is substituted for the argument.

For H1_0 members of this pulse class, there is no parametrically small
version of S6.26's declared RMS/stiffness-proxy hierarchy. That RMS scale
need not exist for a general admitted L-infinity pulse. The result does not exclude long or
approximately bandlimited preparation, selected correlated initial data,
delta-retuned sources, an explicitly delta-dependent effective response,
or a different matching parent. It is not a rolling-gap, cutoff, nonlinear,
quantum, C/D dictionary or UV verdict. Original P8 remains open.
