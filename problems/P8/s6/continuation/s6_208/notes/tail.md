# Original two-created-mode cutoff error

Keep chi_K = 1_(|k|<=K,|-k+P|<=K) in every endpoint term. The contact retains its original one-mode band. Neither a domain derivative nor a one-leg replacement is used.

On FAR, r >= 2U/epsilon implies |P| <= epsilon r/2. Thus |l| <= (1+epsilon/2)r, and removal by either leg forces r > K/2. Integrating the complete r^-4 Cauchy envelope gives, for each endpoint,

2B_j U^(5-j) / [9epsilon^(5-j)K].

On NEAR, r < 2U/epsilon and both created-mode magnitudes are below (2/epsilon+1)U = 201U. A removed pair therefore implies K < 201U. Multiplying the complete near absolute bound by 201U/K controls the entire removed intersection. The highest spatial power rises from five to six, still covered by X46.

On the unexpanded low band, removal implies K < m+|P| = U. The existing low bound therefore gains only U/K. No inverse-radius field estimate or singular momentum frame is applied at zero internal momentum.

The exact weighted far, near and low tail sum is approximately

2.858666666682666667349e53,

strictly below 4e53. This proves the stated original regulator error for K >= m.

Independent quadrature uses the exact angular removed fraction: one for r > K, and otherwise the clipped value [1+(r^2+P^2-K^2)/(2rP)]/2, with the P = 0 and r = 0 limits handled separately. Breakpoints at |K-P|, K and K+P are retained. All near powers, the logarithmic case, the far inverse power and the low massive weight are tested at separated scales and cutoff ratios.

This is a mathematical regulator-removal rate, not a physical cutoff or a bound on the still-unmatched UV/contact sector.
