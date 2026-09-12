# The far region and its original two-leg tail

Define FAR by |k|>=m and |P|<=delta nu(k)/2. Keep the original finite-cutoff indicator1_(|k|<=K,|-k+P|<=K). Taylor expansion is performed on the unregulated coefficient at fixed k, not on this indicator or its moving intersection boundary.

Only the fourth-order Taylor remainder is integrated here. The radial majorant is

J4=integral d^3k/(2pi)^3 nu^-4=Amax^3/(8pi m)<Amax^3/(24m).

Extending FAR to the entire radial integral is an upper bound, not a change to the regulated definition. Therefore the pointwise Fourier estimate is at most4e57 J4 |P|^5 times the detector norm and five source-jet norm. Cauchy-Schwarz/Plancherel and |P|^5<=(1+|P|^2)^3 give

|Rfar|<1e54 ||D||L2 X46[Gamma],
X46^2=sum_(j=0)^4||(1-Delta)^3 partial_t^j Gamma||L2^2.

Time is integrated on the same unit slab, with no global-time claim.

On FAR, nu<=2|k| yields |P|<=delta|k|. Thus |-k+P|<=(1+delta)|k|. If either original internal leg is outside cutoffK, then |k|>K/(1+delta)>K/2. For K>=m,

integral_(|k|>K/2) d^3k/(2pi)^3 nu^-4
<=Amax^4/(pi^2 K)<Amax^4/(9K).

The corresponding error is below1e58||D||L2 X46[Gamma]/K. Dominated convergence therefore removes this regulator for the FAR REMAINDER alone.

Independent angular-radial tests retain the actual removed union. For fixed k,P and k<=K, its angular fraction is(clamp(c,-1,1)+1)/2, c=(k^2+P^2-K^2)/(2kP); for k>K it is1. The far lower limit is max(m,Amax sqrt((2|P|/delta)^2-m^2)) when the square root applies. Quadrature splits at the moving band boundaries. Tests include zero transfer and both transfer regimes, in addition to the exact full J4 and half-band radial comparisons.

Two external metric factors give4/kappa times these estimates:4e-746 and4e-742/K. These small numbers belong to this particular finite weak piece. They do not establish smallness of a full coupled inverse or erase derivative loss.
