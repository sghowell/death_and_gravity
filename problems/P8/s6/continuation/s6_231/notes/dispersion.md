# Full-error optical implication and improved matching tradeoff

Assume the four explicit physical and error premises in FORMULATION.md. In particular the scalar pole mass is1 with positive canonical LSZ normalization, and the pole-subtracted forward amplitude is analytic at crossing point s=2 with the required bounded arc. Let v=s-2 at t0. Combining the two crossing-related cuts gives the twice-subtracted kernel

2v²/[pi (S-2)((S-2)²-v²)].

Its HALF second derivative at v0 is2/[pi(S-2)^3]. Thus physical b2=B''(0)/2 is bounded from below by2/pi times the positive absorptive integral over any elastic annulus above threshold. Properly subtracted isolated poles and other intermediate states are handled by the hypotheses; their unknown values are not silently discarded. No new cut through the crossing neighborhood is allowed.

For S>=16 and abs(x)<=1/2, the lambda part of the full tree is positive and

A_tree >= (9gamma/16)S(S-4)²-8gamma >= gamma S³/4.

The second inequality follows from expanding
(9/16)S(S-4)²-8-S³/4 in S-16; every coefficient is strictly positive. The central interval has width1, so

norm(A_tree)_L2² >= gamma² S^6/16.

This is a lower bound on the FULL angular norm. Its proof does not delete the potential or neglect mass. The complete amplitude is real at tree level; the exact amplitude may be complex.

On every S in[E²/2,E²], the stated FULL-amplitude error and the reverse triangle inequality give

norm(A_exact)_L2 >= (1-eta) norm(A_tree)_L2.

This step includes all corrections in one explicitly unproved norm bound. It is not justified merely by tree coefficients being small. Apply the exact positive optical inequality from notes/normalization.md and beta>3/4 for S>=16. Also(S-2)^-3>S^-3. Retaining only this annulus,

b2 >= (2/pi) integral Im A_exact(S,0)/(S-2)^3 dS
   >= [3(1-eta)² gamma²/(2048pi²)] integral S³ dS
   = 45(1-eta)² gamma² E^8/(131072pi²).

The integral from E²/2 to E² is15E^8/64. Every factorial and the distinction between B'' and b2 are retained. The first inequality may receive additional positive weight outside the annulus. The lower bound is deliberately conservative.

For the actual lambda and gamma, pi²<10 implies the STRICT bound

b2/(4lambda) > 9(1-eta)²(E/10^125)^8.

If abs(b2-4lambda)<=4lambda delta, the upper bound b2<=4lambda(1+delta) yields the necessary inequality1+delta>9(1-eta)²(E/10^125)^8.

At E10^125 with eta<=1/2, the right side is at least9/4. The exact lower bound before the pi estimate is90lambda/pi², hence b2>9lambda. This contradicts b2<=8lambda, or delta<=1. It follows that the stated S-matrix/dispersion premises, the angular error tolerance and the physical b2 tolerance cannot all hold. This result does not select the failing condition.

For comparison, on the full physical angular range and2<=sqrt(S)<=E, absolute values give

abs(A_tree) <= 4lambda E^4+(3gamma/4)E^6+8gamma.

At the named E this is below10^-46. The contradiction is therefore not just the large-amplitude real-part ceiling restated. Conversely this small tree size does not establish the full error premise: matching, higher operators, heavy channels, real counterterms and nonperturbative contributions remain open.
