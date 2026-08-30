# P9(a) — Certified upper bound on H₀ from frozen late-time data

Status: formulation v1, 2026-08-25. This file fixes the statement that the
P9(a) entry point will establish. Changes after the first certificate is
produced go in the Revision log at the bottom, never by silent edit.

## 1. Statement (theorem template)

Fix the frozen inputs D of `data/MANIFEST.json`:

- BAO: the 12 DESI DR2 measurements of D_M/r_d and D_H/r_d at
  z ∈ {0.510, 0.706, 0.934, 1.321, 1.484, 2.33} with their 12×12 covariance
  C_BAO (the BGS D_V/r_d point at z = 0.295 is dropped; see §3.4).
- SN: the n_SN = 1580 Pantheon+ light curves with zHD > 0.01 and
  IS_CALIBRATOR = 0, with apparent magnitudes m_j = m_b_corr, redshifts
  (zHD_j, zHEL_j), and the corresponding 1580×1580 block C_SN of the
  STAT+SYS covariance.
- Sound horizon: r_d ∈ [r_lo, r_hi] = [146.57, 147.61] Mpc (Planck 2018
  r_drag = 147.09 ± 0.26 Mpc, ±2σ box).

Fix the class parameters: z_max = 2.5, X = ln(1+z_max), the grid G in
x = ln(1+z) defined by x_0 = 0, x_1 = ln(1.005), x_{k+1} = 1.1 x_k while the
spacing is below 0.02, then uniform spacing 0.02 up to X (94 nodes; v2 —
geometric at low z so that the interpolation bound of §3.2 is small there),
and a slope bound L > 0. The a priori box 20 ≤ H(z) ≤ 2000 km s⁻¹ Mpc⁻¹ is
part of the class.

**Class.** C(G, L) is the set of continuous functions ũ: [0, X] → (0, ∞)
that are linear on each segment [x_k, x_{k+1}] of G and satisfy, for every
segment,

    |ũ_{k+1} − ũ_k| ≤ L · h_k · min(ũ_k, ũ_{k+1}),   h_k = x_{k+1} − x_k,

together with the box. Here ũ(x) = c / (r_d H(z)) is the dimensionless
Hubble distance in units of r_d. Every member satisfies
|d ln H / d ln(1+z)| ≤ L on each segment, so C(G, L) ⊂ Lip_L(ln H). (Flat
ΛCDM has d ln H/d ln(1+z) = (3/2)Ω_m(z) ≤ 1.5, so it lies in the class for
L ≥ 1.5 approximately; for smaller L the class excludes it.)

**Observables** (flat FLRW):

    D̃_H(z) = ũ(x(z)),
    D̃_M(z) = ∫₀^{x(z)} ũ(x') e^{x'} dx'      (linear in the node values ũ_i),
    D̃_L(zHD, zHEL) = (1 + zHEL) · D̃_M(zHD).

**Fit statistic.** With M' a free nuisance parameter (absorbing the SN
absolute magnitude, the constant 25, and 5 log₁₀(r_d/Mpc)),

    χ²(ũ, M') = (d − P ũ)ᵀ C_BAO⁻¹ (d − P ũ)
              + rᵀ C_SN⁻¹ r,   r_j = m_j − 5 log₁₀ D̃_L(zHD_j, zHEL_j) − M'.

**Feasible set.** For Δ > 0 and an explicit reference point (ũ*, M'*) ∈
C(G, L) × [0, 40] recorded in the certificate,

    F(G, L, Δ) = { (ũ, M') : ũ ∈ C(G, L),  M' ∈ [0, 40],  χ²(ũ, M') ≤ T },  T := χ²(ũ*, M'*) + Δ.

Because T ≥ min χ² + Δ for any reference point, F contains the Δ-sublevel
set of the class minimum; a bound over F is therefore a bound over the
latter. Two reference points are reported: the (locally optimized, exactly
feasible) class minimizer, giving the profile-likelihood statement, and the
flat-ΛCDM best fit when it lies in the class, giving the statement "every
class member fitting at least as well as ΛCDM up to Δ has H₀ ≤ …".

**Claim to be certified.** For all (ũ, M') ∈ F(G, L, Δ) and all
r_d ∈ [r_lo, r_hi]:

    H₀ := c / (r_d · ũ(0)) ≤ H_max(G, L, Δ) := c / (r_lo · ũ₀^min),

where ũ₀^min is a certified lower bound on min{ ũ(0) : (ũ, M') ∈ F }.

The deliverable is the curve L ↦ H_max(G, L, Δ) for
L ∈ {0.5, 1, 1.5, 2, 3, 5, 10} and Δ ∈ {1, 4, 9}, with an exact certificate
for each (L, Δ = 4) row. (L = ∞ is not covered by the v2 relaxation, whose
interpolation bound needs finite L.)

## 2. What is and is not claimed

- Claimed: a statement about the sublevel set F of a fixed, explicitly
  defined statistic, over an explicit finite-dimensional class, for the
  frozen inputs. It is a theorem about those objects.
- Not claimed: frequentist coverage of Δ (no Wilks argument is invoked;
  Δ is part of the definition); anything about models that change r_d
  (early-time physics); anything about the SN absolute calibration (M' is
  free, and the SH0ES calibrators are excluded by construction); anything
  about the continuous class Lip_L(ln H) beyond C(N, L) ⊂ Lip_L (the
  extension to the continuous class is a separate discretization lemma,
  deferred to v2).
- Assumed: spatial flatness; the distance-duality relation
  D_L = (1+z) D_M with the Pantheon+ convention D_L = (1+zHEL) D_M(zHD);
  the released covariances taken as the definition of the likelihood.

## 3. Reduction to a certifiable convex problem

### 3.1 Linearity
All D̃_M(z), D̃_H(z) are linear functionals of the node vector
ũ = (ũ_0, …, ũ_N); the class constraints are linear inequalities; the
objective (minimize ũ_0) is linear. The BAO term is a convex quadratic in ũ.

### 3.2 The SN logarithm (only non-convexity) — v2 relaxation
r_j depends on ũ through log₁₀ of the linear functional D̃_j := D̃_M(zHD_j).

v1 (superseded) sandwiched each ℓ_j = log₁₀ D̃_j between a secant and a
tangent on a per-SN bracket. That relaxation is valid but useless in
practice: 1580 independent slacks let the relaxed χ² fall ~10³ below the
true value until brackets are already at the few-percent level, which
bound-tightening from BAO-only brackets never reaches (documented in
`results/`).

v2 moves the SN term to log-distance space, where it is exactly convex, and
confines the nonconvex link to the grid nodes:

- Node variables κ_i := log₁₀[ D̃_M(x_i) / (e^{x_i} − 1) ] for i ≥ 1 and
  κ_0 := log₁₀ ũ_0 (κ is the log of the running average of ũ; it is smooth
  and tends to log₁₀ ũ_0 as x → 0).
- Per-SN values: ℓ_j = κ(x_j) + log₁₀(e^{x_j} − 1), with κ(x_j) replaced by
  the linear interpolation of the two adjacent node values up to a
  class-only error bound e_k: on segment k, |κ − interp| ≤ (h_k²/8)·B_k/ln 10
  where B_k bounds |F''|, F = ln[D̃_M/(e^x − 1)], over the class. With
  q = e^x/(e^x−1), θ = ũ(x)(e^x−1)/D̃_M(x) ∈ [θ_lo, θ_hi] (from
  |d ln ũ/dx| ≤ L), and a = ũ'/ũ ∈ [−L, L]:
  F'' = q[θa − (θ−1)(q(1+θ) − 1)], bounded factor by factor on each segment;
  the first segment uses the exact form F = ln(ũ_0 + b φ(x)). This gives
  5e_k ≤ 0.002 mag everywhere for L = 1.5 (implemented in
  `geometry.kappa_interp_slack`; Monte Carlo test in `tests/test_kappa.py`).
  The interpolation-error bound holds for C^{1,1} functions, so kinks of ũ
  at nodes are covered.
- Node links: κ_i + c_i = log₁₀ D̃_i with D̃_i = (A_nodes ũ)_i linear in ũ,
  sandwiched between the secant and three tangents on a certified bracket
  [D_lo,i, D_hi,i] (94 sandwiches instead of 1580).
- Smoothness of κ, two families of linear constraints on κ_{i+1} − κ_i:
  (a) class-only: D̃_{i+1}/D̃_i = 1 + inc_i/D̃_i with inc_i ∈ ũ_i[J_−, J_+] and
  D̃_i ∈ ũ_i[I_−, I_+] from the Lipschitz bound (the ũ_i scale cancels), so
  κ_{i+1} − κ_i lies in an explicit interval independent of any bracket
  (a few mmag wide at z ≲ 0.05 for L = 1.5); (b) bracket-aware: the same
  ratio bounded with the current node bounds on ũ_i, ũ_{i+1}, D̃_i. Without
  (a)–(b) the 94 independent link slacks let κ overfit the Hubble diagram and
  bound tightening stalls; with them the relaxed χ² minimum approaches the
  class minimum as brackets shrink.
- BAO stays exact in ũ. Objective min ũ_0. The SN nuisance is boxed,
  M' ∈ [0, 40] (part of the definition of F).
- Bound tightening computes rigorous min/max of ũ_i, of D̃_i through ũ, and
  of κ_i directly (the last converted to D̃-brackets through 10^{κ_i + c_i});
  brackets are the intersection of all valid bounds obtained so far.

Brackets are certified by optimality-based bound tightening: rigorous
min/max of each ũ_i and D̃_i over the current relaxed set (which contains F),
iterated; the SN term no longer degrades with bracket width, so the
iteration converges. Bound-tightening passes are run on a redshift-stratified
subset of the SNe (its marginal χ² is ≤ the full χ², so the subset set
contains F); the reference point and the final bound use the full sample.

### 3.3 The certified problem
After relaxation the problem is

    minimize ũ_0  subject to  A y ≤ h  (class, brackets, sandwich),
                             q(y) := ‖w‖² ≤ T,  L_C w = B y − m,

with y = (ũ, ℓ, M') and C = L_C L_Cᵀ the Cholesky factorization of the
block-diagonal covariance diag(C_BAO, C_SN) (the certificate is stated for
C̃ := L_C L_Cᵀ, with ‖C − C̃‖_max recorded). This is a second-order cone
program.

### 3.3a What "the covariance" means in the certified statement
Every certified statement uses C̃ := (WᵀW)⁻¹ with W the float64 inverse
Cholesky factor of the selected covariance block (SN: 1580×1580; BAO:
12×12); W is the recorded object, treated as exact. Recorded deviations
(`data/MANIFEST.json`, "whitening"): max |W C Wᵀ − I| = 6.2×10⁻¹⁵ (SN),
2.2×10⁻¹⁶ (BAO), i.e. C̃ agrees with the released covariance to double
precision. All file values are parsed as correctly rounded doubles and
treated as exact rationals thereafter; the derived covariance cache is
sha256-pinned in the manifest and checked at load time.

### 3.4 Dropped BGS point
D_V = (z D_M² D_H)^{1/3} is concave in (D_M, D_H), so an upper constraint
on it is not convex. Removing the BGS row and column from the BAO Gaussian
is a relaxation (the marginal χ² is ≤ the full χ²), hence valid. A sandwich
treatment can restore it in v2.

### 3.5 Sound horizon
r_d enters only through H₀ = c/(r_d ũ_0) because M' absorbs 5 log₁₀ r_d
in the SN term and the BAO data are already in units of r_d. Hence the
maximum over the r_d box is attained at r_lo.

## 4. Certificate format and verification

A certificate for a row (N, L, Δ) is a JSON/plain-text file under
`certificates/` containing, as exact rationals:

1. the reference point (ũ*, M'*) and T = χ²(ũ*, M'*) + Δ (χ² recomputed
   exactly by the verifier);
2. the brackets [D_lo,j, D_hi,j] and tangent points used in the sandwich,
   together with the certificates of the brackets (same format, recursively,
   or a single outer LP certificate);
3. multipliers λ ≥ 0 (quadratic constraint), μ ≥ 0 (linear inequalities),
   ν (equalities L_C w = B y − m);
4. the claimed bound ũ₀^min.

The verifier checks, in exact rational arithmetic with no floating point on
the verification path, the polynomial identity

    ũ_0 − ũ₀^min = λ (T − ‖w‖²) + μᵀ (h − A y) + νᵀ (L_C w − B y + m) + s(y, w)

with s a nonnegative quadratic (its ũ-, ℓ-, M'-linear coefficients must
vanish exactly and its constant term must dominate ‖L_Cᵀ ν‖²/(4λ)). A PASS
establishes ũ_0 ≥ ũ₀^min on F_rel ⊇ F by weak duality, hence the H₀ bound.

Ledger levels for this problem:
- CONJECTURED: numerical curve from the cvxpy pipeline (`src/`).
- VERIFIED_N: the same curve from an independent implementation
  (different solver and code path, e.g. SciPy or Julia) agreeing to 0.01
  km s⁻¹ Mpc⁻¹.
- CERTIFIED: exact verifier PASS for the (L, Δ=4) rows; certificate files
  committed.

## 5. Known-answer tests (must pass before any ledger entry above HEURISTIC)

1. Flat ΛCDM restricted fit (Ω_m, h r_d) to the frozen BAO alone reproduces
   the DESI DR2 published values (h r_d = 101.54 ± 0.73 Mpc, Ω_m = 0.2975 ±
   0.0086) within the quoted uncertainties.
2. With Pantheon+ added, Ω_m ≈ 0.31 and β = c/(r_d H₀) ≈ 29.5–29.8
   (Knox–Millea folklore value 29.54 ± 0.41 from BOSS+Pantheon).
3. The ΛCDM best fit lies in C(50, 1.5) and χ²_min over the class is ≤ its
   χ².
4. For smooth classes (L ≤ 2) the bound lands in the range reported by
   model-independent reconstructions (H₀ ≈ 68–69.5, e.g. Lemos et al. 2019,
   Zhou et al. 2025); for L → ∞ it rises, reproducing the low-z loophole
   (Efstathiou 2021).

## 6. Plan

- S0 (done 2026-08-25): data frozen, manifest with sha256 and upstream commits.
- S1: numerical pipeline; ΛCDM tests; curve H_max(L, Δ). → CONJECTURED.
- S2: independent re-implementation of the relaxed problem. → VERIFIED_N.
- S3: exact certificate generation and verifier. → CERTIFIED.
- v2 items: BGS D_V sandwich; DES-SN5YR and Union3 variants; BBN-based r_d
  variant; discretization lemma to Lip_L(ln H); monotone (NEC) subclass
  ũ_{i} ≤ ũ_{i−1}.

### 6.1 v2 options implemented (2026-08-29): D_V row, SN samples, r_d box
All options are off by default: the baseline statement (§1) and its
certificate directories (`results/certificates/lkr_L*_D*_r*`, empty variant
suffix) are unchanged. The drivers `run_lkr_certified`, `certify_feasible`
and `replay` take `--sn {pantheon,dessn5yr,union3}`, `--dv` and
`--rd_box LO HI | planck | bbn`; the certificate tag gets the suffix
`_<sn>`, `_dv`, `_rd<lo>-<hi>` for each non-default part (`src/p9/variants.py`),
and state.json records the variant so that `replay` rebuilds the same inputs.

- **BGS D_V row** (`ClassSpec.use_dv`, `--dv`; BAO covariance = the full
  13×13, `data.load_desi(drop_dv=False)`). In log variables
  log10 D_V = (log10 z_b + 2 y_b + y_H)/3 with y_b = log10 D_M(z_b)
  interpolated from κ (± e_k, exactly as for the D_M rows) and
  y_H = log10 D_H(z_b). Since D_H(z_b) = (1−t) ũ_k + t ũ_{k+1} is a convex
  combination of the node values, concavity of log10 gives the class-only rows
  (1−t) λ_k + t λ_{k+1} ≤ y_H ≤ (1−t) λ_k + t λ_{k+1} + gap_k(t), with
  gap_k(t) = max_{|d| ≤ lh_k} [log10((1−t) + t·10^d) − t d] attained at
  d = ±lh_k (the bracket is convex in d; ≈ 7×10⁻⁵ for L = 1.5, h = 0.02),
  a subset of the coarser bracket [min(λ_k, λ_{k+1}), max(λ_k, λ_{k+1})] ⊂
  [λ_k − lh_k, λ_k + lh_k]. y_V is a variable with the equality
  3 y_V − 2 y_b − y_H = log10 z_b, bracketed by OBBT like y_b (base bracket
  [log10(ũ_lo z_b), log10(ũ_hi z_b)], since D_M ∈ ũ[·] z_b and D_H ∈ [ũ_lo, ũ_hi]),
  and the prediction P_V is sandwiched (four tangents below, chord above) on
  10^{y_V}. All rows are generated once in `lkr_rows.build` for both
  arithmetics (§3.7); the exact statistic, the class minimizer's gradient, the
  rigorous χ² enclosure and the ΛCDM fit use D_V = (z D_M² D_H)^{1/3} exactly.
  Tests (`tests/test_variants.py`): Verifier3 reproduces the solver duals with
  the row on; every row and every absorption bound holds at random class
  members; the 13-row BAO-only ΛCDM fit reproduces DESI DR2
  (Ω_m = 0.2975, h r_d = 101.54 Mpc, known-answer test §5.1).
- **SN samples** (`--sn`). The likelihood form is unchanged for every sample:
  r_j = m_j − 5 log10[(1+zHEL_j) D̃_M(zHD_j)] − M', M' ∈ [0, 40] absorbing the
  sample's magnitude/offset convention. DES-SN5YR: the Dovekie release
  (repository main, commit 442c248, 2025-11-14; 1820 SNe with
  zHD ∈ [0.025, 1.14]; m := MU, released for a fixed M_0; the STAT+SYS
  matrix is released as a precision matrix P). Union3: the 22 binned distance
  moduli of Rubin et al. 2023 with their inverse covariance (z_bin ∈ [0.05,
  2.26]; the μ's carry an arbitrary constant; zHEL := zHD := z_bin). For
  precision-released samples the recorded whitening is W = Lᵀ with
  P_sym = L Lᵀ, so WᵀW = P to rounding and C̃ = (WᵀW)⁻¹; no inverse is formed
  on the certified path (deviations recorded in `MANIFEST.json`, "whitening").
  Joint flat-ΛCDM fits to the frozen 12-row BAO + SN (h r_d in Mpc):
  Pantheon+ Ω_m = 0.304, h r_d = 101.07; DES-SN5YR 0.306, 100.91;
  Union3 0.303, 101.15 (β = c/(r_d H₀) = 29.66, 29.71, 29.64).
- **r_d box** (`--rd_box`). The certificate bounds ũ₀; r_d enters only through
  H₀ = c/(r_d ũ₀) and the (inactive) a priori H box, so for every r_d in the
  box H₀ ≤ H_max · r_lo/r_d exactly; `report_certified` adds the column at
  r_d = 147.09 Mpc and rounds upper bounds up, lower bounds down. The
  BBN-based box (`--rd_box bbn`, `MANIFEST.json` "r_drag_BBN"):
  r_d = 148.04 ± 1.24 Mpc from the DESI DR2 fitting formula (arXiv:2503.14738
  Eq. 2) with ω_b = 0.02218 ± 0.00055 (Schöneberg 2024, arXiv:2401.15054),
  N_eff = 3.044 and ω_bc = Ω_m h² from the DESI DR2 BAO(+BBN) ΛCDM fit;
  box = ±2σ = [145.56, 150.52] Mpc.

### 6.2 Variant results (2026-08-30; L = 1.5, Δ = 4, r = 2, r_d box Planck ±2σ)
Each row is a two-sided certified statement of the §1 form for its own inputs
and its own T = upper(χ²(reference) + 4) (`results/summary.md`; chains
`results/certificates/lkr_L1.5_D4_r2_<suffix>/state.json`, feasible points
`feasible_L1.5_D4_r2_<suffix>.json`). Values at r_lo = 146.57 Mpc; lower
bounds rounded down, upper bounds rounded up.

| SN sample | BGS D_V | max_F H₀ | class-min χ² (n_SN + n_BAO) |
|---|---|---|---|
| Pantheon+ (baseline) | no | [69.5227, 69.8418] | 1379.73 (1580 + 12) |
| Pantheon+ | yes | [69.7713, 70.0769] | 1380.94 (1580 + 13) |
| DES-SN5YR | no | [71.0206, 71.2683] | 1625.65 (1820 + 12) |
| DES-SN5YR | yes | [70.9831, 71.2130] | 1625.66 (1820 + 13) |
| Union3 (22 bins) | no | [72.5207, 72.5870] | 13.51 (22 + 12) |
| Union3 (22 bins) | yes | [72.6384, 72.7052] | 14.05 (22 + 13) |

Reading: the shape wall is sample-dependent at the 1–3 km s⁻¹ Mpc⁻¹ level
(the SN sample fixes the *shape* of D_M(z) that the class must reproduce;
the BAO rows then fix its scale in units of r_d), and adding the BGS D_V row
moves the wall by ≤ 0.3 (up for Pantheon+ and Union3, slightly down for
DES-SN5YR). With Union3 the wall reaches ≈ 72.6 — still below the
local-distance-ladder values quoted in the problem statement. The width of
each two-sided enclosure (0.07 Union3, 0.23–0.25 DES-SN5YR, 0.31–0.32
Pantheon+) is the relaxation gap at r = 2; the class maximum lies inside it.
All numbers are for the Planck ±2σ r_d box; the BBN box rescales the upper
bounds by 146.57/145.56 (§1: H₀ ≤ H₀max · r_lo/r_d).

### 3.6 v3: the λ–κ–ρ relaxation (the one used for certificates)
The v2 relaxation (§3.2) is valid but stalls at H₀ ≲ 75 for L = 1.5 because
the link between the SN side (κ) and the class side (ũ) is loose wherever
BAO does not pin D̃ (`notes/relaxation-log.md`). v3 removes the link
gaps by putting each constraint in the coordinates where it is exact:

- λ_i := log₁₀ ũ_i. The class is exactly |λ_{i+1} − λ_i| ≤ log₁₀(1 + L h_i)
  plus the box (this is the same class as §1: the min-form and the ratio-form
  of the slope bound are equivalent).
- κ_i as in v2; the SN term is exactly convex in κ (interpolation slack e_k).
- Segment identities D̃_{i+1} − D̃_i = a_i ũ_i + b_i ũ_{i+1} in scale-free
  exponential form: 10^{δ_i} − 1 = a_i 10^{ρ_i} + b_i 10^{ρ_i + Δλ_i} with
  δ_i = log₁₀(D̃_{i+1}/D̃_i), ρ_i = log₁₀(ũ_i/D̃_i), Δλ_i = λ_{i+1} − λ_i
  (i = 0: 10^{δ_0} = a_0 + b_0 10^{Δλ_0}). Each exponential is sandwiched
  (four tangents below, chord above) on a bracket: δ from the class-only ratio
  bounds, Δλ from the class, ρ from the class θ-range and then bound tightening.
- BAO: D̃_M(z_b) = 10^{y_b}, y_b interpolated from κ (± e), sandwiched on a
  bracket tightened by OBBT; D̃_H(z_b) = (1−t)10^{λ_k} + t 10^{λ_{k+1}} with
  the two node exponentials sandwiched on λ brackets.
- Objective min λ_0. Bound tightening: rigorous min/max of ρ_i (all nodes),
  of λ at the D_H nodes and node 0, and of the six y_b.

Diagnostics at L = 1.5, Δ = 4 (94-node grid): the relaxed χ² floor is
1362.7 vs class minimum 1380.2; the entire gap is the per-SN interpolation
slack e_k (used at 99.6% by every SN, coherently: first-order effect
2·5·Σ_j e_j |(C⁻¹r)_j| ≈ 17), node-level slack ≈ 0. Since e_k ∝ h_k², the
class grid is refined by midpoint insertion (refinement r: 2^r segments per
segment of G; r = 2 gives 373 nodes and a first-order slack ≈ 1.0). The
coarse class is a subset of the refined class, so certificates for the coarse
class (feasible points) remain valid for the refined one and bounds for the
refined class bound the coarse one.

### 3.7 Certificate chain (v3)
All rows of the relaxed conic program are generated by one code path
(`src/p9/lkr_rows.py`) parameterized by the arithmetic — floats for the
solver, Arb balls for the verifier (`verify3.py`) — so structure and order are
identical by construction. The chain: class-box brackets (definition) →
each tightening solve certified by weak duality with residual absorption →
outward-rounded brackets → … → final bound on λ_0 certified. T is the
rigorous upper enclosure of χ² at an exactly class-feasible reference point
plus Δ. Feasible points (lower bounds on max H₀) are certified by exact class
membership plus a rigorous χ² enclosure (`certify_feasible.py`).

## Revision log
- v1 (2026-08-25): initial statement (uniform grid N = 50; per-SN sandwich).
- v2 (2026-08-25): geometric-at-low-z grid G (94 nodes); κ (log-distance
  node) relaxation replacing the per-SN sandwich; reference point defined as
  the polished class minimizer (ΛCDM reference reported alongside); L = ∞
  dropped from the deliverable. Reason: the v1 relaxation does not converge
  (see §3.2).
- v3 (2026-08-25): λ–κ–ρ relaxation (§3.6) with midpoint-refined grid
  (r = 2, 373 nodes) as the certified relaxation; shared-row certificate
  chain (§3.7); M' box made explicit in §1. The class definition is unchanged
  up to the grid refinement (a superset).
