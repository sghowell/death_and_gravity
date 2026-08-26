# Claims ledger

Rules: one row per claim; a claim's level changes only together with an
evidence path; REFUTED is terminal; dates are ISO. Levels: HEURISTIC,
CONJECTURED, VERIFIED_N, CERTIFIED, FORMALIZED, REFUTED (see README).

| id | problem | statement | level | evidence | updated |
|---|---|---|---|---|---|
| P4-1 | P4 | The Evans–Coleman CSS solution of Einstein–Euler (p = ρ/3) crosses its sonic point (a saddle) along the branch V₁ > 0 with V₀ = 0.112439401388092; V has exactly one zero; t = 0 limits a → 1.0653, m/r → 0.05945 | CONJECTURED | problems/P4/notes/numerics-report.md §B; tests/test_shoot.py | 2026-08-25 |
| P4-2 | P4 | Linearized CSS problem (KHA sonic-point gauge): exactly one non-gauge eigenvalue in [0,15]×[−14,14], κ₁ = 2.8105525488 (γ = 0.355801922); gauge eigenvalue 0.3556992037; next eigenvalue κ₂ = −1.5456 | CONJECTURED | numerics-report.md §C (argument principle); tests/test_perturb.py | 2026-08-25 |
| P4-3 | P4 | On the closed sound cone of the EC profile, V_z increases monotonically to 1/√3, ρ̂ decreases, w = 1/3 − V_z² decreases (GHJS hypotheses hold in these variables); V_R is not monotone | CONJECTURED | numerics-report.md §D; results/ec_monotonicity.json | 2026-08-25 |
| P9-1 | P9 | For the class C(G, L=1.5), M' ∈ [0,40], frozen inputs, and T = χ²(ũ*) + 4 with the recorded reference point: an explicit class member with χ² ≤ T has H₀ = 69.5313 at r_lo, i.e. max_F H₀ ≥ 69.5313 | CERTIFIED (exact class check + ball-arithmetic χ² enclosures; certificate = the point) | problems/P9/results/certificates/feasible_L1.5_D4.json; src/p9/certify_feasible.py | 2026-08-25 |
| P9-1b | P9 | Same statement for L ∈ {1, 2, 3, 5, 10} (Δ = 4) and Δ ∈ {1, 9} (L = 1.5): float-verified feasible points H₀ ≈ 69.5–69.8; rigorous certificates being generated | CONJECTURED | results/certify_feasible_all.log | 2026-08-25 |
| P9-2 | P9 | Certified upper bound H₀ ≤ 75.40 for C(G, 1.5), Δ = 4 via the v2 relaxation (dual certificates not yet stored) | HEURISTIC | results/run_L1.5_D4_v3.json | 2026-08-25 |
