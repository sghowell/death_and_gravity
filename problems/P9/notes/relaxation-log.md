# P9(a) relaxation log (2026-08-25)

All numbers for L = 1.5, Δ = 4, geometric grid (94 nodes), full Pantheon+ (1580 SNe) unless noted.
Ground truth anchors: class-minimum χ² = 1380.213 (T = 1384.213); flat-ΛCDM χ² = 1399.74;
a verified-feasible class member reaches H₀ = 69.53 (at r_lo), so the true maximum is ≥ 69.53.

| relaxation | relaxed χ² floor | H₀ bound | verdict |
|---|---|---|---|
| v1: per-SN secant/tangent sandwich, uniform grid N=50, BAO-only starting brackets | 335 | 99.5 | 1580 independent slacks; SN term vacuous until brackets ≲ ±5%; bound tightening stalls |
| v2: κ nodes (log running-average distance), per-SN interpolation with rigorous class-only error, 94 node links | 1293 | 77.6 | link slack (~0.13 mag/node at factor-3 brackets) lets κ overfit; stall |
| v2 + class-only κ-difference bounds | 1340 | 75.5 | low-z overfitting removed; mid-z links (D unconstrained by BAO, factor-2 brackets) still leak ~0.07 mag → hockey stick at z<0.05 |
| v2 + bracket-aware differences, 3 tangents, OBBT on κ | 1341 | 75.4 | same mechanism |
| κ-only (no ũ): exact SN, class-only first/second differences, BAO via 6 exponential sandwiches | 1344 | 76.1 (subset tightening) | class-only linear conditions lose the slope–history coupling; also: subset tightening is invalid as a *useful* outer set (subset χ² floor ≈ 250 vs T ≈ 1384 → BAO ±30σ) |
| λ–κ–ρ (LKR): exact class in λ = log₁₀ ũ, exact SN in κ, segment identities in scale-free exponential form | 1362.7 | 72.22 | converges in 3 full-sample passes; node-level slack 0; the whole floor gap is the per-SN interpolation slack (∝ h²), used coherently (first-order estimate 17.0; midpoint ×2 → 4.1, ×4 → 1.0) |
| LKR on the ×4 refined grid (373 nodes), certified chain | (running) | (running) | expected floor within ~1 of the class minimum |

Certified anchors (all rigorous): max_F H₀ ≥ 69.5313 (L=1.5, Δ=4), ≥ 68.3481 (Δ=1), ≥ 70.5647 (Δ=9);
≥ 69.7527 (L=1), 69.5835 (L=2), 69.7414 (L=3), 70.1621 (L=5), 71.5592 (L=10) at Δ=4.

## Lessons
1. Any relaxation with per-object slack that scales with an *absolute* bracket width cannot bootstrap
   from data-free brackets when the data noise per object exceeds the slack needed to bind.
2. Put each constraint in the coordinates where it is exact: the class in log ũ, the SN term in log D,
   BAO in D. The only necessary nonconvexity is then the map between log ũ and log D, and it can be
   written with *differences* (δ = log₁₀ D_{i+1}/D_i, ρ = log₁₀ ũ_i/D_i, Δλ) whose brackets come from the
   class alone and are narrow.
3. Subset-of-SNe outer sets are valid but useless for quantities that need the χ² threshold to bind:
   T is a full-sample number. Bracket tightening must be done on the full sample; keep the number of
   tightened quantities small.
4. Always compute a verified-feasible point with large H₀ alongside the bound: it measures the
   relaxation gap and prevents mistaking a loose bound for physics.
