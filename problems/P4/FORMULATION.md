# P4(b) — Evans–Coleman self-similar radiation-fluid collapse: existence and unstable-mode count

Status: formulation v1, 2026-08-25, written from `notes/literature-digest.md`
(equation references there). Items marked **[open]** must be settled by the
S1 numerics before the CAP statement is frozen; the theorem templates below
are the targets, not yet frozen statements.

## 1. Objects

**System.** Spherically symmetric Einstein–Euler with p = kρ, k = 1/3,
continuously self-similar (CSS) ansatz. Two equivalent reductions are used:

- (E) the Guo–Hadžić–Jang comoving/Schwarzschild self-similar 2×2 system in
  x = r̃(y) for (D, W) (digest §1.4, GHJ eqs. 2.35–2.38), with sonic locus
  B(x) = 0 — used for existence (autonomous, rational right-hand side,
  closest to the Buckmaster–Cao-Labora–Gómez-Serrano sonic-point template);
- (P) the Harada–Maeda comoving Misner–Sharp system (M, S, η, y) in
  z = r/(−t) with its algebraic constraint (digest §1.2, HM01 eqs. 30–35),
  whose linearization (HM01 eqs. 83–87) defines the eigenvalue problem; the
  Koike–Hara–Adachi polar-areal system (digest §1.1) is the cross-check for
  eigenvalues (their κ, with s = −ln(−t), is the same growth rate as HM01's λ).

**Definition (Evans–Coleman solution).** A CSS solution of (E)/(P) such that
(i) the centre is regular (W(0) = 1/3; y(0⁺) = 1/3; analytic in local
Cartesian coordinates), (ii) it crosses its unique sonic point analytically,
and (iii) V_R has exactly one zero (Hunter (a) type; digest §2). At k = 1/3
the sonic point is a saddle, so C¹ regularity there already implies
analyticity and the crossing is along one of exactly two analytic branches
(digest §1.5); the EC branch is the first-order root that is *not* the
Larson–Penston-type root (digest §5, GHJ Lemmas 4.7–4.10). **[open]**: fix
which root of GHJ's quadratic Q is the EC branch at ε = 1/3 by numerics.

**Domain.** The critical exponent depends only on the eigenvalue problem on
the closed backward sound cone [0, x_sp] (KHA integrate from the sonic point
to the centre; free boundary at the sonic point since no information enters
from outside; digest §7a). The existence theorem is therefore stated on
[0, x_sp] with analytic continuation a short distance past x_sp; extension
to t = 0 (ξ → ∞, digest §2) is a separate statement, deferred.

## 2. Theorem templates

**Theorem A (existence, CAP).** There exist x_sp in an explicit interval and a
real-analytic solution (D, W) of (E) on [0, x_sp + δ] with W(0) = 1/3, an
analytic crossing of the sonic point at x_sp along the EC branch, and exactly
one zero of V_R on (0, x_sp]; it is the unique such solution with x_sp in that
interval (local uniqueness from the interval-Newton/Krawczyk step). Its
profile lies within explicit interval bounds of the numerical profile.

**Theorem B (mode count in a rectangle, CAP).** Fix the gauge of digest §3
(HM01 comoving gauge, in which the pure-gauge eigenvalue is λ_g = (1−k)/(1+k)
= 1/2; or KHA's sonic-point gauge with λ_g ≈ 0.35699 — **[open]**: choose by
which linearization is better conditioned for validated integration). Let
E(λ) be the matching function whose zeros are the eigenvalues of the
linearized CSS problem with regular-centre and analytic-sonic-point boundary
conditions on [0, x_sp]. Then in the rectangle R = {0 ≤ Re λ ≤ R₁, |Im λ| ≤
R₂} with (R₁, R₂) = (15, 14) (the KHA search box) the zeros of E are exactly
λ_g and one simple real zero λ₁ ∈ [2.81055, 2.81056]; hence, conditional on
Theorem C, the PBH mass-scaling exponent is γ = 1/λ₁ ∈ [0.355801, 0.355802].

**Theorem C (exclusion outside the rectangle, analytic).** E has no zeros
with Re λ ≥ 0 outside R. **[open — research item]**: the GHJS argument for
Re λ > 1 and |Im λ| ≫ 1 uses monotonicity of the Larson–Penston profile
(ρ̂' < 0, ω̂' > 0, hence w'/w < 0 inside the cone). Whether the EC profile
at k = 1/3 has the analogous monotonicity on [0, x_sp] (|V_z| increasing from
0 to √k, density decreasing) is the first thing S1 must determine; if it
does, the GHJS energy method is the template; if not, a WKB/Frobenius-index
argument for |λ| → ∞ is needed.

The near-term deliverable is A + B (a computer-assisted existence theorem and
a rectangle-restricted mode count — publishable on its own, with C stated as
the remaining conjecture). The full universality statement requires C.

## 3. What is and is not claimed
- Claimed (A+B): existence and local uniqueness of an EC-type CSS solution at
  k = 1/3; the exact count of non-gauge eigenvalues in R; a rigorous
  enclosure of λ₁ and hence of γ.
- Not claimed: global uniqueness of the EC solution within the CSS family
  (only local uniqueness in x_sp); nonlinear codimension-one stability;
  non-spherical modes (Gundlach 2002: stable for 1/9 < k ≲ 0.49); the
  extension of the solution to t = 0 and beyond; kink modes (excluded by the
  analytic/Sobolev class; EC is kink-stable for k < 0.89, digest §2).

## 4. Reduction to certifiable pieces

A1. Sonic-point data: zeroth order from the algebraic conditions (digest
    §1.1 formulas at γ = 4/3, or GHJ Lemma 4.1 at ε = 1/3); first order from
    the quadratic Q after factoring the ghost root; validated Taylor
    coefficients (D_n, W_n) at x_sp with explicit growth bounds
    (GHJ Thm 4.18 / App. A style; Buckmaster et al. recursion with
    resonance check).
A2. Centre data: Frobenius expansion at x = 0 with W(0) = 1/3, one free
    parameter; validated coefficients and tail bound. **[open]**: analytic
    variable at the centre (x or x²; digest §1.2 notes non-integer powers
    in z that become integer at k = 1/3).
A3. Shooting: integrate from both ends to an interior matching point with a
    validated ODE solver (VNODE-LP as in GHJS, or Arb/CAPD Taylor models);
    the matching condition F(x_sp) = 0 is verified by interval Newton
    (existence + uniqueness in an interval).
A4. Zero count of V_R: sign certificate along the validated profile.
B1. Linearization in the chosen gauge; Frobenius exponents at the centre
    and at the sonic point as functions of λ (the second sonic-point
    exponent must be shown non-admissible, cf. GHJS Lemma 3.6).
B2. Validated Taylor data for the eigenfunction at both singular points for
    complex λ, with growth bounds; validated integration to a matching
    point; E(λ) = Wronskian, evaluated on interval boxes.
B3. Argument principle on ∂R (or box subdivision with a Wronskian-nonzero
    certificate, as GHJS's `intermediate_evalue_excluder`), plus a Krawczyk
    enclosure of λ₁ and identification of λ_g by matching its eigenfunction
    to the gauge generator (digest §3, HM01 eqs. 96–103).
C.  Analytic large-|λ| exclusion (research; see §2).

## 5. Ledger levels and known-answer tests
- HEURISTIC → CONJECTURED (S1): a high-precision (mpmath) reproduction of
  the EC profile at k = 1/3 in both (E) and (P) variables; the sonic point as
  a saddle; the two first-order branches; κ₁ = 2.81055255 to ≥ 7 digits
  (KHA99 Table 2); the gauge eigenvalue (1/2 in HM01 gauge; 0.35699 in KHA's
  sonic-point gauge); the next eigenvalues (KHA99: Re κ ≲ −1.4, then a
  complex pair); EC94's t = 0 numbers (a → 1.07, m/r → 0.0596) as a check on
  the continuation. Plus the monotonicity diagnostics for Theorem C.
- VERIFIED_N: the same numbers from the independent (P) vs (E)/(KHA) codes.
- CERTIFIED: Theorem A; Theorem B for the stated rectangle.
- FORMALIZED: not planned.

## 6. Plan
- S0 (done): literature digest; this formulation.
- S1: numerics (profile, spectrum, monotonicity diagnostics); freeze the
  branch choice, gauge, and rectangle. Output: `notes/numerics-report.md`,
  `src/p4/`.
- S2: Theorem A (CAP existence). Tooling decision: VNODE-LP (C++, GHJS code
  as scaffold) vs Arb/Julia Taylor models.
- S3: Theorem B (rectangle mode count).
- S4: Theorem C (analytic exclusion) — research track; may become its own
  problem entry.

## Revision log
- v1 (2026-08-25): initial, with [open] items to be closed by S1.
