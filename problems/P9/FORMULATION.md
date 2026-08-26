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

Fix the class parameters: z_max = 2.5, X = ln(1+z_max), a uniform grid
x_i = iX/N (i = 0..N) in x = ln(1+z) with N = 50, and a slope bound L > 0.

**Class.** C(N, L) is the set of continuous functions ũ: [0, X] → (0, ∞)
that are linear on each segment [x_{i−1}, x_i] and satisfy, for every
segment,

    |ũ_i − ũ_{i−1}| ≤ L · h · min(ũ_{i−1}, ũ_i),   h = X/N.

Here ũ(x) = c / (r_d H(z)) is the dimensionless Hubble distance in units of
r_d. Every member satisfies |d ln H / d ln(1+z)| ≤ L on each segment, so
C(N, L) ⊂ Lip_L(ln H). (Flat ΛCDM has d ln H/d ln(1+z) = (3/2)Ω_m(z) ≤ 1.5.)

**Observables** (flat FLRW):

    D̃_H(z) = ũ(x(z)),
    D̃_M(z) = ∫₀^{x(z)} ũ(x') e^{x'} dx'      (linear in the node values ũ_i),
    D̃_L(zHD, zHEL) = (1 + zHEL) · D̃_M(zHD).

**Fit statistic.** With M' a free nuisance parameter (absorbing the SN
absolute magnitude, the constant 25, and 5 log₁₀(r_d/Mpc)),

    χ²(ũ, M') = (d − P ũ)ᵀ C_BAO⁻¹ (d − P ũ)
              + rᵀ C_SN⁻¹ r,   r_j = m_j − 5 log₁₀ D̃_L(zHD_j, zHEL_j) − M'.

**Feasible set.** For Δ > 0 and an explicit reference point (ũ*, M'*) ∈
C(N, L) × ℝ recorded in the certificate,

    F(N, L, Δ) = { (ũ, M') : ũ ∈ C(N, L),  χ²(ũ, M') ≤ T },  T := χ²(ũ*, M'*) + Δ.

Because T ≥ min χ² + Δ for any reference point, F contains the Δ-sublevel
set of the minimum; a bound over F is therefore a bound over the latter.

**Claim to be certified.** For all (ũ, M') ∈ F(N, L, Δ) and all
r_d ∈ [r_lo, r_hi]:

    H₀ := c / (r_d · ũ(0)) ≤ H_max(N, L, Δ) := c / (r_lo · ũ₀^min),

where ũ₀^min is a certified lower bound on min{ ũ(0) : (ũ, M') ∈ F }.

The deliverable is the curve L ↦ H_max(50, L, Δ) for
L ∈ {0.5, 1, 1.5, 2, 3, 5, 10, ∞} and Δ ∈ {1, 4, 9}, with an exact
certificate for each (L, Δ = 4) row, plus N ∈ {25, 100} as a sensitivity
check (not certified).

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

### 3.2 The SN logarithm (only non-convexity)
r_j depends on ũ through log₁₀ of the linear functional D̃_j := D̃_M(zHD_j).
Relaxation: introduce ℓ_j and replace ℓ_j = log₁₀ D̃_j by the sandwich

    secant_j(D̃_j) ≤ ℓ_j ≤ tangent_j(D̃_j)   on a bracket D̃_j ∈ [D_lo,j, D_hi,j],

where the bracket is itself certified: [D_lo,j, D_hi,j] are rigorous lower
and upper bounds of the linear functional D̃_j over a convex outer set
containing F (first the BAO-plus-class set, then, iteratively, the current
relaxed set). Since log₁₀ is concave, the sandwich contains the graph, so the
relaxed set F_rel ⊇ F and any bound over F_rel is valid for F. The gap is
O((ln(D_hi/D_lo))²) in log₁₀ D and shrinks under iteration.

### 3.3 The certified problem
After relaxation the problem is

    minimize ũ_0  subject to  A y ≤ h  (class, brackets, sandwich),
                             q(y) := ‖w‖² ≤ T,  L_C w = B y − m,

with y = (ũ, ℓ, M') and C = L_C L_Cᵀ the Cholesky factorization of the
block-diagonal covariance diag(C_BAO, C_SN) (the certificate is stated for
C̃ := L_C L_Cᵀ, with ‖C − C̃‖_max recorded). This is a second-order cone
program.

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

## Revision log
- v1 (2026-08-25): initial statement.
