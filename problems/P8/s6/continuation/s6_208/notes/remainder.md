# Complete far, near and low estimates

Let c_jr be the full S6.198 Leibniz/Cauchy rows, C = 1000 and R = 1/20000. Each normalized source-jet coefficient is bounded by

9 C^2 c_jr R^(-(j-r)).

A common complete-row bound is B_j = 9 C^2 sum_r c_jr R^-j. The five B_j are

9000000,
540000000000,
50400000000000000,
6552000000000000000000,
1101600000000000000000000000.

The factor nine retains every physical pair. Source derivatives are taken with the detector time fixed.

For FAR, r >= L(P) = 2U/epsilon, Cauchy gives endpoint remainder

2B_j epsilon^(j-5) U^(5-j) r^-4.

Its actual far radial integral is below B_j U^(4-j)/(18 epsilon^(4-j)).

For NEAR, m <= r < L(P), use the complete all-real row bound 4e27 nu < 8e27 r. Its integral is below 4*REAL*U^4/(9epsilon^4), REAL = 4e27. Retain every UV coefficient separately. For j+d < 4, its integral is bounded by

B_j 2^(4-j-d) U^(4-j) / [18(4-j-d)epsilon^(4-j)].

For j+d = 4, the logarithm log(L/m) is retained and bounded by L/m, giving

2B_j U^(5-j) / [18m epsilon^(5-j)].

All fifteen slots and thirty-five source time-jet entries remain. No odd endpoint or finite oversubtraction is discarded. The unexpanded low band remains below 2e38.

For every integer q <= 6,

U^q <= 8m^q(1+|P|^2)^3.

The exact weight identities and m >= 1 establish this uniformly. Plancherel and the complete source-jet norm give the stated X46 bound. The exact summed coefficient is approximately 1.422222222422222222560e48, strictly below 2e48. This includes the entire complementary region, not just the far Cauchy remainder.
