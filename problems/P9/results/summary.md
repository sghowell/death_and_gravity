# P9(a) certified results

Class C(G, L) on the midpoint-refined grid (refine r), M' ∈ [0, 40], frozen DESI DR2 BAO + the SN sample of the row
(C̃ as recorded); T = upper(χ²(reference) + Δ). Upper bounds: chains of Arb-verified conic dual certificates;
lower bounds: certified feasible class members. Values in km s⁻¹ Mpc⁻¹. The certified bound is on ũ₀ = c/(r_d H₀);
r_d enters only through H₀ = c/(r_d ũ₀), so for every r_d in the row's box H₀ ≤ H₀max · r_lo/r_d exactly, and the
'at r_d = 147.09' column is H₀max · r_lo/147.09 (rescaling to any other r_d is the same one-line identity).
Baseline rows: Pantheon+, BGS D_V row dropped, r_d box = Planck 2018 ±2σ (FORMULATION §1); variants per §6.1.

| L | Δ | r | SN | D_V | r_d box [Mpc] | T | certified H₀ max (upper, at r_lo) | at r_d = 147.09 | certified feasible H₀ (lower) | passes | status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.5 | 4 | 2 | dessn5yr | no | [146.57, 147.61] | 1629.6473 | **≤ 71.2683** | ≤ 71.0164 | ≥ 71.0206 | 5 | done (lkr_L1.5_D4_r2_dessn5yr) |
| 1.5 | 4 | 2 | dessn5yr | yes | [146.57, 147.61] | 1629.6604 | **≤ 71.6855** | ≤ 71.4320 | - | 2 | running (lkr_L1.5_D4_r2_dessn5yr_dv) |
| 1 | 4 | 2 | pantheon | no | [146.57, 147.61] | 1600.6025 | **≤ 69.8430** | ≤ 69.5961 | ≥ 69.6497 | 6 | done (lkr_L1_D4_r2) |
| 1.5 | 1 | 2 | pantheon | no | [146.57, 147.61] | 1380.7252 | **≤ 68.9288** | ≤ 68.6851 | ≥ 68.3373 | 5 | done (lkr_L1.5_D1_r2) |
| 1.5 | 4 | 2 | pantheon | no | [146.57, 147.61] | 1383.7252 | **≤ 69.8418** | ≤ 69.5949 | ≥ 69.5227 | 5 | done (lkr_L1.5_D4_r2) |
| 1.5 | 9 | 2 | pantheon | no | [146.57, 147.61] | 1388.7252 | **≤ 70.7971** | ≤ 70.5468 | ≥ 70.5562 | 5 | done (lkr_L1.5_D9_r2) |
| 2 | 4 | 2 | pantheon | no | [146.57, 147.61] | 1381.1495 | **≤ 70.0717** | ≤ 69.8240 | ≥ 69.5876 | 5 | done (lkr_L2_D4_r2) |
| 3 | 4 | 2 | pantheon | no | [146.57, 147.61] | 1376.8755 | **≤ 70.6184** | ≤ 70.3688 | ≥ 69.7561 | 6 | done (lkr_L3_D4_r2) |
| 5 | 4 | 2 | pantheon | no | [146.57, 147.61] | 1370.2487 | **≤ 71.9731** | ≤ 71.7187 | ≥ 70.2001 | 6 | done (lkr_L5_D4_r2) |
| 10 | 4 | 2 | pantheon | no | [146.57, 147.61] | 1358.4565 | **≤ 76.2758** | ≤ 76.0061 | ≥ 71.6439 | 7 | done (lkr_L10_D4_r2) |
| 1.5 | 4 | 2 | pantheon | yes | [146.57, 147.61] | 1384.9350 | **≤ 70.0769** | ≤ 69.8291 | ≥ 69.7713 | 5 | done (lkr_L1.5_D4_r2_dv) |
| 1.5 | 4 | 2 | union3 | no | [146.57, 147.61] | 17.5083 | **≤ 72.5870** | ≤ 72.3304 | ≥ 72.5207 | 5 | done (lkr_L1.5_D4_r2_union3) |
| 1.5 | 4 | 2 | union3 | yes | [146.57, 147.61] | 18.0521 | **≤ 72.7052** | ≤ 72.4482 | ≥ 72.6384 | 5 | done (lkr_L1.5_D4_r2_union3_dv) |
