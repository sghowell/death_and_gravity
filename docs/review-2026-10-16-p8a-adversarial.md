# Independent adversarial review: P8(a) rows A.8–A.23

Date: 2026-10-16. Reviewer: a separate agent with no part in authoring
P8(a), given the code and the ledger rows but not the authors'
conclusions. The review was read-only; the repository was unchanged
before and after.

Scope: CLAIMS rows P8-A.16 to P8-A.23 were reviewed in full (formulations,
proof notes, verify and audit code, independent re-derivations in scratch
scripts). P8-A.8 to P8-A.15 were triaged. A.1–A.7 were reviewed earlier in
the Empiricist trial and were treated as inputs here.

**Overall verdict:** the chain survives as correct but narrow conditional
FLRW theorems, built mostly on known QEIs. It does not resolve the open
null or non-minimal problems in general. "CERTIFIED" overstates what is
machine-checked.

## 1. Replays

All replays were run with `PYTHONDONTWRITEBYTECODE=1` and
`pytest -p no:cacheprovider`, so nothing was written.

- **A.1 verifier:** `PYTHONPATH=problems/P8/a/src .venv/bin/python -m p8a.verify --check`
  exits 0: "A.1 … replay passed; physical application OPEN".
- **A.1 tests:** `pytest problems/P8/a/tests -q`, 36 passed.
- **A.2–A.19:** all 18 documented `--check` commands pass, 1–8 s each.
- **A.20–A.23:** their READMEs give no complete command. Each was run with
  every `problems/P8/**/src` on `sys.path`, as their conftests do, and all
  4 pass.
- **Sub-package tests:** all 23 pytest suites pass, 1,541 tests in total.
  A.20/21/22/23 have 175/110/139/123 tests, matching the ledger.
- **Not run:** ruff and the full-P8 snapshot runner.

## 2. Mathematics: verdict per row

No mathematical error was found in A.16–A.23.

- **A.23 — PROVEN AS STATED, as a conditional theorem.**
  - Re-derived: the screen identity I+θ₀ = 2∫b²((h/b)′)², the operators P,
    U, V, the Cauchy–Schwarz majorants and every constant. τI ≤ 1.43051,
    the margin is 0.56949 > 1/2, and the exact rationals reproduce.
  - The theorem is equivalent to: no SEE solution meeting (P) and (Q)
    extends past λ = 2τ while obeying (F) and the future Wick cap. Those
    future caps are assumptions about the very region whose existence is
    in question.
  - "Incompleteness" is incompleteness of the FLRW chart only. The
    contracting flat-sliced de Sitter patch satisfies every hypothesis for
    any β_S when 2 ≤ |H₀|τ ≤ 3, yet extends to complete global de Sitter.
    There A is linear in x, A″ = 0, w = 0, and T_conf ∝ g is absorbed
    into Λ.
  - The only exhibited members of the hypothesis class (that patch and the
    β_S = 0 thermal branch) satisfy the NEC. No NEC-violating member is
    shown.
  - Presentation-level GAP: scope and novelty wording.
- **A.22 — PROVEN AS STATED.**
  - Verified: identity (5), the mass-shell polynomial (8) by direct
    integration, the Gram matrices (10)/(11),
    Γ(r) = 35r + 24/r + 21/r³ with Γ(4/3) = 14117/192, the transport of L,
    M, G, and the reference null contraction (4).
  - At ξ = 0 the coefficient ħ‖f_tt‖²/(12π²) coincides with the
    Fliss–Freivogel–Kontou DSNEC at boost η = 0 (conversion of their
    Eq. 49).
- **A.21 — PROVEN AS STATED.** Checked the mode normalization, the n = 1/3,
  m = −2/3 covariance, the stress and Wick-square expressions, and the
  segment bounds.
- **A.20 — PROVEN AS STATED.**
  - Re-derived the coefficient (3−4ξ)/48 → 7/144 (Fewster–Osterbrink 2008,
    Thm 4.3, massless) and the transport F = a^{-3/2}f.
  - Checked the anomaly reference: ρ = 3H⁴, trace 12H²(H²+Ḣ). The Weyl²
    variation vanishes on this background, and β_S covers the R²
    ambiguity. δ = ℓ_P²/(πτ²) is consistent.
  - The margin 0.180929 exceeds 9/50 by only 9.3×10⁻⁴.
- **A.19 — PROVEN AS STATED.**
  - Checked H² − bH⁴ = κQ/(3a⁴), a_end⁴ = 16λ(1−4λ), and
    EED ≥ (3H²/κ)(1−2bH²) > 0.
  - This is the classic anomaly-plus-radiation model (Fischetti–Hartle–Hu,
    PRD 20, 1757, 1979); the repo cites Fröb–Glavan–Meda 2025 for the
    exact equation. The endpoint sits at Planck curvature.
- **A.18 — PROVEN AS STATED.** The history tube gives H ≤ −1.9131/τ and
  d = (2.1, 9, 70, 800); the margin is 47079/350000.
- **A.17 — PROVEN AS STATED.**
  - Recomputed V = 16704, C_P + C_F ≈ 1.2566×10⁷, and margins 0.7994
    (σ = 0) and 0.4243 (σ ≤ 1).
  - Presentation-level GAP: the named calibration has no SEE witness.
    Classical radiation violates |Ḣ|τ² ≤ 4.
- **A.16 — PROVEN AS STATED.** Re-derived L_H and the 4601/1280 radiation
  integration by parts. The content is the known Maxwell QEI (Pfenning;
  Fewster–Pfenning, both cited in the repo) plus conformal transport.
- **A.8–A.15 — TRIAGED.** Replays pass. Formulations for A.11, A.12, A.14
  and A.15 were read.
  - Every "actual SEE" domain is shorter than a Planck time. With
    Q₂ = κħ/(8π²) = ℓ_P²/π, A.12's Q₂/τ² ≥ 1.6×10⁸ gives
    τ_avail ≲ 4.5×10⁻⁵ t_P, and A.15's Q₂/τ² ≥ 2/5 gives
    τ_avail ≲ 0.9 t_P.
  - A.13's positive pole is about one e-fold per Planck time: the known
    Horowitz runaway (PRD 21, 1445, 1980).
  - The approach cannot reach macroscopic times without order reduction
    (Parker–Simon, PRD 47, 1339, 1993).
  - The chain proves no incompleteness; A.12 and A.15 report the focusing
    test failing.
  - In what was read, nothing relies on results A.1–A.7 did not establish.
    A.17 onward deliberately bypass A.1's initial-Ricci premise.

## 3. Certification: claim-level GAP for every reviewed row

README.md defines CERTIFIED as a machine-checkable certificate. Here the
machine checks cover symbolic identities, rational constants and
polynomial bounds only. The analytic core is written prose: positive-type
and Hadamard restriction, H² density, Banach and Hadamard propagation,
and the global endpoint arguments.

- **Literal `True` proof gates:** A.23 has 16 of 28, including
  `initial_outgoing_expansion…`, `actual_SEE_reduces_plane_total_stress…`
  and `mixed_quantum_cross_term…`. A.22 has 16 of 18, including
  `Hadamard_pullback…`. A.21 has 7 of 17 and A.20 has 6 of 23.
- **"Rejected inputs"** are type and range checks: True, 0.1, "0", oo,
  nan, a Symbol, −1.
- **A.23's 123 tests** are 44 report-field mutation tests, 24 bad-input
  tests, 4 serialization tests, 35 re-runs of the identities, and 16
  others.
- **Identity counts include tautologies**, for example `a/a−1`, and
  `2Hd/a²−2a2/a` where Hd is defined as a·a2.
- **Snapshot counts** are repo-wide, not specific to the row.

## 4. Novelty of A.20–A.23

- **A.20:** the inequality is Fewster–Osterbrink's bound transported to
  FLRW. Fewster–Kontou 2019 (arXiv:1809.05047) already gives
  state-dependent QSEIs for ξ ∈ [0, 2ξ_c] in curved spacetime. The
  coherent-state obstruction is also known from both papers. What is new
  is a concrete FLRW theorem with explicit budgets; the step is
  incremental (compare Fewster–Kontou 2022's toy cosmology,
  arXiv:2108.12668).
- **A.21:** Fewster–Roman's mechanism (PRD 67, 044003) extended to
  ξ ∈ [0, 1/4] with w ≤ 0. It is consistent with Fliss–Freivogel–Kontou–
  Pardo Santos 2024 (arXiv:2309.10848), whose field-value class is
  two-sided and whose single-line bound also needs a UV cutoff. Minor.
- **A.22:** built on Fliss–Freivogel–Kontou 2023 (arXiv:2111.05772) and
  Fliss–Freivogel–Kontou–Pardo Santos 2024. The real improvement is small:
  the (Df)² weight makes a one-sided Wick cap rigorous, where FFKP use a
  sign-indefinite ∂²(f²) weight.
- **A.23:**
  - No prior paper with this exact statement was found; the search was
    not exhaustive.
  - Closest prior work: Freivogel–Kontou–Krommydas 2022 (arXiv:2012.11569),
    a Penrose-type theorem from SNEC; Moghtaderi–Hull–Quintin–Geshnizjani
    2025 (arXiv:2503.19955), SNEC limits on FLRW bounces; Wall 2013
    (arXiv:1010.5513), the GSL implies null incompleteness of spatially
    infinite cosmologies; Fewster–Kontou 2020 (arXiv:1907.13604), the
    index-form method.
  - It evades the null obstruction rather than solving it. In FLRW the SEE
    makes the total null stress T_KK a function of t alone, so no
    null-line average is ever needed. For ξ = 0, a worldline QEI for the
    null-contracted stress (Fewster–Roman) would already suffice. For
    ξ = 1/6, any spatial smearing that integrates out Wick-square
    gradients works.
  - The genuinely open null case, in generic spacetimes, is untouched.
    Nothing here contradicts Fewster–Roman, Fewster–Kontou,
    Graf–Kontou–Ohanyan–Schinnerl, or the SNEC/DSNEC papers.

## 5. Wording fixes recommended

1. Replace CERTIFIED with something like "WRITTEN PROOF; identities and
   constants CERTIFIED". Drop literal-`True` gates, mutation tests and
   type tests from evidence counts, and remove repo-wide snapshot counts
   from the rows.
2. A.17, A.18, A.20 and A.23: state each theorem as "no extension past τ
   (or 2τ) under the future caps", and add that the incompleteness is of
   the FLRW chart, which may be extendible (de Sitter example).
3. A.23: replace "no homogeneous quantum state assumed" with "the SEE
   forces a homogeneous total stress; only the state and its Wick square
   may be inhomogeneous".
4. A.22: replace "second-width shrinking limit diverges" with "the lower
   bound diverges as δ₋ → 0".
5. Credit prior art in the rows: Fewster–Osterbrink and Fewster–Kontou
   2019 (A.20); Fewster–Roman (A.21); Fliss–Freivogel–Kontou and FFKP
   (A.22); Freivogel–Kontou–Krommydas, Moghtaderi et al. and Wall (A.23);
   Pfenning and Fewster–Pfenning (A.16); Fröb–Glavan–Meda and
   Fischetti–Hartle–Hu (A.19); Horowitz (A.13).
6. A.11–A.15: state the proper durations in Planck units.
7. Say that existence of solutions in the hypothesis class is shown only
   for β = 0 (the A.19 and A.20 thermal branches). For A.23 the de Sitter
   patch covers every β_S.
8. Root README: "photon/global-flat-FLRW objective is complete" overreaches
   P8(a)'s "assuming only QEIs", because FLRW symmetry and the future caps
   are extra hypotheses.
9. `problems/P8/a/README.md` is stale because A.1's certificate hash-pins
   it. Unpin it, or document the newer rows elsewhere.
