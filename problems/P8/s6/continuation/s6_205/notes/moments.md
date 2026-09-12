# All finite massive radial moments

For s=j+n-1 in4,...,7, extend the positive high-band majorant to all real k only for estimation. With nu=sqrt(m^2+k^2/Amax^2),

J_s=integral d^3k/(2pi)^3 nu^-s
   =Amax^3 m^(3-s) Gamma((s-3)/2)/(8pi^(3/2)Gamma(s/2)).

The exact values are

J4=Amax^3/(8pi m),
J5=Amax^3/(6pi^2 m^2),
J6=Amax^3/(32pi m^3),
J7=Amax^3/(15pi^2 m^4).

Using pi>3 yields rational upper denominators24m,54m^2,96m^3 and135m^4, each with numeratorAmax^3. The actual finite coefficients still integrate only over |k|>=m. No low-mode term, physical cutoff or original comparison partition is changed.

Summing D_jn times these rational bounds over every finite cell gives approximately3.143626e49, strictly below1e50. All source time jets are included in D_jn and every physical pair is included in the factor9.

The resulting pointwise Fourier bound is a sum of |P|^n, n<=4. Plancherel and Cauchy-Schwarz with the already stated X46 norm give1e50||D||L2 X46[Gamma]. No new regularity requirement or same-space inverse claim is introduced.

Independent trigonometric-substitution quadrature checks all four full moments. Exact-order guards run before memoization, including after valid entries have prewarmed the cache. Candidate UV orders are not passed to a finite-moment formula.
