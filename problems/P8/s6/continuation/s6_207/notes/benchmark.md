# Nonzero exact massive flat benchmark

The benchmark is flat, with a = 1, k = r e3 and l = (-r+P)e3. It retains the full massive three-polarization Proca field and normalized tracefree tensor, vector and scalar spatial channels. It is not the curved current or a contact-matched result.

Let x = 1/r and define

Omega_k = sqrt(1+m^2x^2),
Omega_l = sqrt((1-Px)^2+m^2x^2),
S = Omega_k+Omega_l,
A = Omega_k Omega_l+m^2x^2,
B = 1-Px.

The complete normalized j = 0 densities are

tensor: (A+B)^2/(4 Omega_k Omega_l S),
vector: m^2x^2 S/(4 Omega_k Omega_l),
scalar: [(A-B)^2/3+2A^2/3]/(4 Omega_k Omega_l S).

The flat source time-jet coefficient at endpoint j multiplies the corresponding expression by cos(j*pi/2)/S^j. This static cancellation of odd endpoints is not transferred to the curved calculation.

The exact tensor j = 0 coefficients through degree four are

1/2,
-P/4,
m^2/4-P^2/8,
m^2 P/8-P^3/16,
-m^4/16+m^2 P^2/16-P^4/32.

The vector coefficients are

0, 0, m^2/2, m^2 P/4, -m^4/4+m^2 P^2/4.

The scalar coefficients are

1/12,
-P/24,
5m^2/24-P^2/48,
5m^2 P/48-P^3/96,
5m^4/32+3m^2 P^2/32-P^4/192.

All remaining retained endpoint coefficients are included in the report. There are 45 channel/endpoint/power entries. In particular, the fourth-time-jet tensor and scalar coefficients are 1/32 and 1/192. The mixed vector and longitudinal scalar contributions are not discarded.

A degree-four exact polynomial recurrence checks the square-root branches, multiplication and inverse phase. Independent full ten-field polarization sums, not these closed forms, supply the numerical Cauchy coefficients at two transfers. All 45 entries are checked at each transfer.
