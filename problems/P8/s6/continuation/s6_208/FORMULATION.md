# Uniform actual endpoint remainder and exact repartition

Retain the same CD massive Proca sector, m = 1000, kappa = 10^800, prepared state, unit compact time slab, full unit-W8 comparison and fixed finite covariant prescription. Let r = |k|, x = 1/r, P be external spatial momentum, U = m+|P| and epsilon = 1/100.

The new analytic radius is rho(P) = epsilon/U. For real unit n and complex |x| <= rho, the scaled momenta n and -n+xP share the original outer/inner time discs 1/10000 and 1/20000. All four W8 frequency coefficients and both analytic Schwarz signs are retained.

On the inner disc, the normalized ten-field norm is below 20 and the complete sixteen-component pair norm is below 1000. Every normalized inverse summed phase is below one. These improved constants concern normalized analytic modes only; the existing all-real unscaled bound remains separate.

Let E5(P,k) denote the complete pre-current first-five-endpoint row and let U5(P,k) be its actual inverse-radius UV symbol from S6.207, retaining degrees d = 0,...,4-j for each endpoint j. Every source time jet through order j remains.

At the original mask chi_K = 1_(|k|<=K,|-k+P|<=K), define

Qnew,K = integral chi_K [
  1_(r<m) E5(P,k) + 1_(r>=m)(E5(P,k)-U5(P,k))
] d^3k/(2pi)^3,

with source/detector contractions and time/Fourier integrations understood. No test, Borel-state coefficient or sharp mask is differentiated.

Split the high band at L(P) = 2U/epsilon. Cauchy controls the far remainder. The complementary near region retains the full real row and every UV power/log term. The low band is unexpanded. Together they give

|Qnew|, |Qnew,K| < 2e48 ||D||L2 X46[Gamma],
|Qnew-Qnew,K| < 4e53 ||D||L2 X46[Gamma]/K, K >= m,

where X46^2 = sum_(r=0)^4 ||(1-Delta)^3 partial_t^r Gamma||L2^2. The largest near-tail spatial power is six, so this same norm controls both statements.

The exact actual current is repartitioned as

Jactual,K = (R_actual-unit,K + F_unit,K + Qnew,K)
           + (C_unit,K + U5,K).

The state correction and finite time remainder retain their S6.201 bound 2e48 and tail 2e52/K. The new known actual piece therefore has bound 5e48 M[D]Y[Gamma] and tail 5e53 M[D]Y[Gamma]/K, with the unchanged Y^2 = N61^2+X46^2.

This is not Qnew added to the older known endpoint piece: that would double-count finite terms. The former Q203+Pfin+PUV and new Qnew+U5 are both exactly E5 at the same finite mask. The fixed local target is not added as a matched quantum contribution.

Both canonical metric factors give endpoint bounds 8e-752 and 16e-747/K, and known-actual bounds 2e-751 and 2e-746/K. Complete UV/contact finite and divergent coefficients, subleading sharp-band artifacts and original fixed covariant matching remain open. No full response/inverse/background/stability/cutoff or original V/G/B closure follows.
