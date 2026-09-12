# Original two-leg regulator error

Keep the original high band |k| >= m and the exact mask

chi_K(P,k) = 1_(|k| <= K, |-k+P| <= K).

No one-leg replacement is used. For a finite cell let s = j+n-1, so s is in {4,5,6,7}.

When |P| <= K/2, removal by either created-mode leg implies |k| > K/2. Since nu >= |k|/Amax,

integral_(|k|>K/2) d^3k/(2pi)^3 nu^-s
 <= Amax^s 2^(s-3) K^(3-s) / (2pi^2(s-3))
 < Amax^s 2^(s-3) m^(4-s) / (18(s-3)K),

for K >= m. The last estimate retains all four radial orders.

When |P| > K/2, use the full finite moment and 1 <= 2|P|/K. This adds one external spatial power. Since n+1 <= 5, the same X46 norm controls this regime and the low-transfer regime.

Sum every D_jn times the low-transfer coefficient and twice the full-moment coefficient. The resulting rational constant is approximately 5.582731e54, strictly below 1e55. Thus

|P_fin - P_fin,K| < 1e55 ||D||L2 X46[Gamma] / K.

The complete angular removed-union fraction is independently integrated after the substitution |k| = Amax m tan(theta). Splitting at the geometric thresholds avoids integrating through unhandled changes in the intersection. Tests cover every finite cell, both internal cutoff scales, and external ratios 0, 1/4, 1/2, 1, 2 and 100. Separate radial tests check each order and the half-band estimate up to K = 10^12.

This is removal of a mathematical regulator, not a physical cutoff determination. The contact keeps its separate one-mode band and is not covered by this finite-polynomial tail argument.
