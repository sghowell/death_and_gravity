# Integrated sunset Schwinger bounds

Write Qnorm=16 pi^2. The scalar two-loop Schwinger normalization is
Qnorm^-2 integral d^N alpha U^-2 exp(-F/U). The denominator
certificate bounds the exponential by exp(-sum alpha/4).
There is no large-mass expansion or physical momentum cutoff.

For each of the six individually UV-finite different-heavy sunsets,
N=5 and U has degree two. Set alpha=4 beta; the factor is four.
If J_box is the unit-cube integral of U(beta)^-2, homogeneity gives
integral over [0,t]^5 = t J_box. Layer-cake integration and
exp(-sum beta)<=exp(-max beta) bound the full integral by
J_box integral_0^infinity t exp(-t)dt=J_box.

All 120 edge orderings for each graph are enumerated. In the ordering
beta_0>=...>=beta_4, the chosen co-tree positions dominate every
other co-tree componentwise; this is checked for every ordering.
U is at least that monomial. The monomial bound has powers -2
at its two positions and zero elsewhere. Every tail exponent
a_j=sum_{i>=j}(power_i+1) is positive, with minimum one.

The unit ordered simplex integral is product_j 1/a_j.
An independent triangular substitution beta_j=t_0...t_j checks
its Jacobian and all monomial powers, and each endpoint integral
is evaluated independently. Summing gives J_box=14 for every
one of the six refinements. There are 720 checked edge orders.
Their combined Wick weight is one. Their unprojected absolute
bound is 56 g^2/Qnorm^2. On-shell Cauchy projection from radius
two to radius one gives 28 g^2/Qnorm^2 times |s-1|^2.

For the local sunset U=ab+ac+bc and P=abc. Its three proper
Phi4 contractions leave p-independent tadpoles. Two derivatives
with respect to s annihilate these local contributions and the
overall mass/kinetic counterterms. The remaining Schwinger
integrand is P^2 U^-4 exp(-F/U).

Crucially, retain P^2: replacing it by a coarse multiple of U^2
sum(alpha)^2 would lose the proper-boundary convergence.
Scaling alpha=4 beta again gives four. The homogeneous
unit-cube growth is again t, with the same layer-cake integral.
All six orderings have powers (-2,-2,2), positive tails (1,2,3),
and simplex integral 1/6; the cube upper bound is one.

Thus |Pi_local''| <= (L^2/6) 4/Qnorm^2. Integrating twice
along the segment from s=1 with the two fixed on-shell zeros
gives L^2/(3 Qnorm^2) times |s-1|^2. Uniform positive
tail exponents and the first-sheet gap justify differentiation
and regulator removal after this grouping.
