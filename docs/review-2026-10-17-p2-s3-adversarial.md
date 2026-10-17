# Independent adversarial review: P2 stage S3 (part (iii), validated zero counts at P_B, l = 0–19)

Date: 2026-10-17. Reviewer: a separate agent with no part in authoring P2 S3.

Scope: branch `p2/s3-part-iii` at `db2a37eb`. The S3 commits are:
- `89a2336b` (12:29:48): the pre-registration, committed together with `zsol.py`, `part3.py`,
  `zeq_identity.py` and `test_part3.py`;
- `3e7fca05` (12:31:11): the C0 deviation;
- `c7d1bfbc` (12:33:04): the control gate and the per-piece |W| diagnostic;
- `44904fbb` … `a3b12803` (12:45–13:30): the certificates as they landed;
- `1abe9d71`, `09013c78`, `db2a37eb`: the report, FORMULATION v0.5 and the commit table.

The review read the author's files without changing any of them. The reviewer's scripts, outputs and
README are in `problems/P2/certificates/review_s3/`; the script names below refer to that directory.

**Overall verdict:** the count statement is proven as stated, and nothing is blocking.
- **The counts hold.** For each l = 0, …, 19 at P_B:
  - W_l is analytic on a neighbourhood of R′_l = [−a_l, a_l] × [−b, e];
  - W_l has no zero on ∂R′_l;
  - its winding number around ∂R′_l is 1 for l = 0 and 0 for l ≥ 1.

  With W₀(0) = 0, R′₀ contains exactly one zero, the simple zero at ω = 0.
- **How this was checked.**
  - Every analytic step checks, including the reviewer's own sympy derivation.
  - The code is rigorous on a line-by-line reading.
  - Every stored number was reproduced by independent means (§5), and no check disagrees anywhere.
- **Level.**
  - CERTIFIED fits the counts on R′_l.
  - The strip corollary ("S_δ contains no QNM of degree ≤ 19 other than 0") uses P2-0 for Im ω > 0,
    so it is at P2-0's level (WRITTEN_PROOF), exactly as in the P2-2 and P2-1w rows.
  - The proposed row files the corollary under CERTIFIED without that qualifier (S1).
- **The l = 0 enclosure is correct, and so is β ≤ 0.90938098.** The P2-K6 amendment is acceptable with
  small corrections (§7).
- **The pre-registration was followed.**
  - The C0 change is a legitimate, stricter correction of the control, and it is disclosed.
  - The same commit also replaced a test point. That change was unnecessary, weakens the test, and is
    not disclosed (S2).
- **Should-fix items:**
  - S1: the ledger wording;
  - S2: the undisclosed test change;
  - S3: two cited checks have no committed output;
  - S4 (optional): no test exercises the new box enclosure.

## 1. Verdicts

| item | verdict |
|---|---|
| A1. Pre-registration against execution | FOLLOWED. Region, method, constants, l-range, controls and kill criteria are as registered. The pipeline code was committed with the pre-registration. Afterwards only the C0 reference call, metadata and a 3-line diagnostic changed. The author's 2351 working records match the certificates one for one, with no re-runs or splits. C7 and the extra l = 0 check have no committed output (S3) |
| A2. The C0 deviation (`3e7fca05`) | LEGITIMATE and STRICTER; disclosed, with two exceptions. The test point (7, 40.0) became (7, 30.0) without disclosure (S2), and the commit message names the wrong failing point (W6). The original C0 would have caught a point-value error at 21 of the 24 points. At the 3 high-frequency points its overlap would have been vacuous |
| B3. Phase factor | NONE. W_zsol equals the frozen W_l exactly (factor ≡ 1). The zsol normalisation at each horizon is the frozen one |
| B4. Analyticity; the true zero counts; D4 | PROVEN. No pole within distance (κ₊ − b) of R′_l. Winding 1 together with W₀(0) = 0 makes ω = 0 the only zero in R′₀, and a simple one |
| B5. Conjugation symmetry | PROVEN, and checked explicitly on the left half for every l |
| B6. Rigour of the enclosures | SOUND. This covers the zsol series tails, the Lohner deviation recursion, the box bound M, the Cauchy aliasing and truncation bounds, the piece increments, the \|W\| lower bound and the integer rounding. zsol does not share arbw's large-\|ω\| failure mode |
| B7. The corollary on S_δ | CORRECT, at P2-0's level (S1). s3-report §6 does not overclaim relative to P2-3 |
| B8. The l = 0 Krawczyk enclosure | CORRECT. Re-run identical. A certified sign change across the stored interval with two codes. β ≤ 0.90938098 holds. The code is "independent of KAT-6" only in its Wronskian (W8) |
| C9. Recomputation | ALL REPRODUCE (§5): l = 0, 1, 2, 19 replayed in full; l = 0, 1, 2, 5 recounted rigorously by the reviewer's counting layer on the full boundary; all 2351 increments reproduced by an independent float code |
| C10. Controls | They exercise the whole pipeline for counting, discrimination and symmetry. No control exercises the box bound M (S4); the reviewer's `box_checks.py` does |
| C11. Tests | Fast: 76 passed, 13 skipped (18 s). P2_SLOW=1: 89 passed, 1 pre-existing warning (425 s) |

## 2. Findings by severity

### Blocking

None. Checked without finding a flaw:
- the identity W_zsol = W_l and the horizon normalisations;
- the pole clearance of R′_l and of every disc;
- the Frobenius and Taylor tail bounds;
- the Lohner recursion and the box bound;
- the Cauchy-disc error terms;
- the piece logic and the coverage of H;
- the symmetry argument;
- the integer rounding;
- the coverage of S_δ with P2-0 and P2-2;
- the Krawczyk step.

### Should-fix before the merge

- **S1. The ledger level of the corollary** (s3-report.md:148; FORMULATION v0.5, the "S3 result"
  bullet).
  - The proposed P2-3B row states "With P2-0 and P2-2: for every l ≤ 19, S_δ contains no QNM of
    degree l other than ω = 0", and files the whole row under "CERTIFIED".
  - That corollary needs P2-0 (WRITTEN_PROOF) for Im ω > 0. This includes the band 0 < Im ω ≤ e with
    |Re ω| ≥ a_l, which neither P2-2 nor R′_l covers.
  - The P2-2 row ("at P2-0's level") and the P2-1w row ("The Im ω > 0 part is at P2-0's level")
    already mark this, so P2-3B should too.
  - Two further changes to the row:
    - drop "not yet independently reviewed";
    - cite analyticity of W_l off E₊ ∪ E_c to FORMULATION §1.2 rather than to P2-D1′. P2-D1′ is the
      resonance equivalence, and it is not used here.
  - Use the text of §7.
  - Failure scenario: a reader of the ledger takes "no QNM of degree ≤ 19 in S_δ other than 0" as
    machine-certified. In fact its Im ω > 0 part rests on P2-0's written argument. A gap there would
    leave a zero with 0 < Im ω ≤ e and |Re ω| ≥ a_l excluded by no certificate.
- **S2. An undisclosed change in the deviation commit `3e7fca05`.**
  - The commit also changed `tests/test_part3.py::test_zsol_matches_arbw`, replacing its third point
    (l, Re ω/κ₋, Im ω/κ₋) = (7, 40.0, −0.55) by (7, 30.0, −0.55).
  - s3-report §3 says "The fix touched only the control".
  - The dropped point does not make the original C0 raise: unwrapped arbw at 256 bits gives a relative
    radius of 2.5e−21 there (`c0_detail.json` (a)). It also passes the new C0 at 256 bits, in 0.3 s
    (`c0_and_krawczyk.json` C0b). So the change was unnecessary, and it lowers the frequency the test
    reaches.
  - Restore 40.0, or record the change and its reason in s3-report §3.
  - Failure scenario: the fast suite is the routinely run check. A later change to zsol that breaks
    point values only above about 30κ₋ at moderate l would pass it, although the registered 40κ₋
    point catches it. More generally, a deviation commit that also edits a test without saying so
    cannot be audited from the report.
- **S3. Two cited checks have no committed output.**
  - s3-report §3 and the proposed row cite C7 ("132 segments recomputed; every increment overlaps")
    and "an extra, unregistered check at l = 0 [that] gives π(1.000 ± 2.9e−4) on the left half".
  - C7's "part3 replay: OK" lines exist only in gitignored logs in the author's worktree
    (`results/part3/verify_*.log`). The l = 0 left-half sum appears nowhere in the repository.
  - Commit both outputs, or remove "l = 0" from "the symmetry checked at l = 0 and 2" and cite this
    review's artifacts. These cover the same ground:
    - `indep_winding.json`: the left half traced for every l;
    - `recount.json`: l = 0, 1, 2, 5 certified on the full boundary;
    - `replay_author.json`: l = 0, 1, 2, 19 replayed.
  - Failure scenario: the ledger row cites two checks that nobody can inspect or replay from the
    repository. If either had in fact failed or been run differently, the merged record would not
    show it.
- **S4 (optional). Test coverage of the new rigorous component.**
  - zsol's box enclosure, the sup bound M, is the one genuinely new rigorous piece of S3.
  - No pre-registered control or test checks it against an independent value. C0 compares point
    values only.
  - Because M exceeds the sampled sup|W| by factors between 73 and 6e20 (`box_checks_A.json`), C1–C3
    would still pass with an M too small by orders of magnitude.
  - The replay tests recompute l ≤ 5 only; l ≥ 6 is only re-summed.
  - Suggested additions:
    - a tight-box containment test against S1's arbw at an l = 19 corner, like `box_checks.py B`
      (about 6 s);
    - one or two recomputed l = 19 segments.
  - Failure scenario: a later edit to `propagate` drops a term of the deviation recursion, so that M
    comes out 100 times too small. Some disc is then below its true sup (the tightest margin seen is
    73). Every control and test would still pass, because the actual interpolation errors are far
    smaller than the claimed bounds. But the claimed bounds, and so the increment balls, would no
    longer be proofs.

### Wording and nits

- **W1.** s3-report §2 says the Lohner bounds "give M/|W| ≤ 10¹⁸ for ρ ≤ κ₋ up to l = 19". The §1 table
  and the certificates show max M/|W_c| above 10¹⁸ for every l ≥ 3, up to 1.4e22 at l = 18. The figure
  is the pilot's, copied from the pre-registration. It is harmless, because K adapts (at most 96 of
  the allowed 160).
- **W2.** s3-report §3, C0 row: "(l = 0–19, both edges, up to the corners)". The points are at
  l ∈ {0, 1, 2, 5, 10, 19} only.
- **W3.** FORMULATION v0.5:
  - "∂R_l" should be "∂R′_l".
  - "dyadic b ≥ β" uses β for the strip depth (½+δ)κ₋. That clashes with D3's β = α/κ₋, the quantity
    S3 also bounds (β ≤ 0.90938098). The same clash is in s3-report §1, s3-preregistration §1 and
    `params.Point.beta`. Write "b ≥ (½+δ)κ₋".
  - "The pole checks are recorded per segment" is loose: the certificates store (c, ρ), and the check
    is enforced at compute time. `check_certs.py` re-verifies it.
- **W4.** The FORMULATION header still lists the ledger as "P2-par, P2-0, P2-2 and P2-K6". P2-D1′ and
  P2-1w are missing; the list is stale since v0.4.
- **W5.** The stored balls are Arb decimal strings whose radii are dominated by the rounding of the
  printed midpoint. These are valid enclosures, but:
  - **Krawczyk.** The "± 3.2e−31" is the 30-digit print; K2 places the zero about 4e−34 inside its
    upper end.
  - **Segment increments.** These are printed with 4–7 significant digits. The reviewer's float
    increments sit up to 0.89 radii from the printed midpoints, which is consistent with this.
  - **Format.** FORMULATION §4's "exact dyadic balls" was not followed. Only segment totals are
    stored, not the per-piece increments, so no certificate can be re-checked without recomputation.
    The same holds for KAT-6.
- **W6.** The commit message of `3e7fca05` says arbw "raised at the l = 19 corner". In the
  pre-registered order the first raise is at the l = 10 corner (56.57κ₋, −0.55κ₋); arbw raises at both
  high-frequency l = 19 points too (`c0_and_krawczyk.json` C0a).
- **W7.** The pilot is "about 0.5 core-h" in the pre-registration but "≈ 0.3 core-h" in s3-report §5.
- **W8.** "Independent of KAT-6" holds for the Wronskian code only. The construction, the parameters
  (r, R, h), the starting value and the `Point` module are all KAT-6's (`winding.krawczyk`). The
  report's own wording, "with a different Wronskian code", is the accurate one.
- **W9.** The existing P2-K6 row says "radius 6e−30" for the P_B zero. That is the pre-fix KAT-6
  value: 5.51e−30 at `3be6ad89`. The certificate regenerated at `9f1b8e6f` has 7.10e−30. Correct it
  together with the amendment (§7).

## 3. Pre-registration integrity (scope A)

### A1. The pre-registration against execution

| pre-registered (`89a2336b`) | executed | verdict |
|---|---|---|
| Target P2-3 at P_B, l = 0–19. P_A only as control C2(d) | as registered | ✓ |
| R′_l: a_l = `W_l_exact` of `lg_PB_l{l}.json`; b = (½+δ)κ₋ rounded up to 2⁻³⁰; e = κ₋/10 rounded down | `check_certs.py`, with independently computed parameters: a = `W_l_exact` for all 20 l; b/κ₋ = 0.5500000107; e/κ₋ = 0.0999999801; b < κ₊ | ✓ |
| Method and constants: zsol; Cauchy discs; ρ ≤ κ₋; planner ρ ≤ dist/(1.5√2); θ = 1/4 (box) and 3 (points); 128 bits; K rule capped at 160; x = 1/2; pieces ≤ κ₋/16; depth 14; halving down to κ₋/256 | the same constants in `part3.py:46–55`. Minimum pole clearance dist/(√2ρ) = 1.504 over all segments. No segment needed halving | ✓ |
| C0: overlap at 24 points | changed in `3e7fca05` (A2) | disclosed, stricter; S2 |
| Gate C0–C3 and C6 before production | `controls_PB.json` and `krawczyk_PB_l0.json` stamped `3e7fca05` (clean); 80 s wall; launch at `c7d1bfbc` 12:33:08, clean | ✓ |
| C4, C5, C7 after production | C4 committed (`symmetry_PB_l2.json`). C5 holds. C7 was run (author's gitignored logs) but not committed | S3 |
| Kill criteria K1–K3 | none fired; 4.7 core-h of 24 | ✓ |
| Execution: nice 10, ≤ 4 workers, JSONL records, launch record | nice 15, 4 workers; launch record committed. The 2351 JSONL lines (author's worktree, read only) equal the certificate pieces exactly: one ok record per segment, no failures, splits or duplicates | ✓ |
| Level CERTIFIED, with routine written steps | §7 | S1 |

**Git timeline.**
- `zsol.py` and `part3.py` first appear in `89a2336b`, together with the pre-registration.
- From then to `db2a37eb`, `zsol.py` is unchanged. In `part3.py` only three things change:
  - the C0 reference call (`_arbw_reference`, `zsol_vs_arbw`);
  - the `computed_at` metadata;
  - the 3-line |W| diagnostic of `c7d1bfbc`.
- `cauchy_segment`'s certification logic is otherwise identical.
- Every winding certificate records `computed_at` = `c7d1bfbc` with no `+dirty`.

**The pilot.** The pilot is disclosed (s3-preregistration §5), and it already ran the full windings at
l = 0, 1, 2. The expectations at those l were therefore known; only l = 3–19 were predictions. That is
acceptable.

### A2. The C0 deviation

**What the first run met.** `c0_and_krawczyk.py` C0a replays the original C0: S1's
`arbw.wronskian`, unwrapped, at 256 bits, at the 24 pre-registered points.
- arbw raises "Frobenius: N too small relative to |rho|" at three points: (10, 56.57, −0.55),
  (19, 71.20, 0.1) and (19, 101.71, −0.55).
- In the pre-registered order the first raise is the l = 10 corner (W6).
- At the other 21 points the 256-bit arbw balls are tight (relative radius ≤ 2.6e−27) and overlap
  zsol.

**"At 256 bits its ball is also too wide for an overlap to mean anything."** This is accurate.
- With the term-count wrapper but still at 256 bits, arbw's balls at the three points contain 0.
- Their radii are 4.9e4, 1.1e29 and 2.5e78 (`c0_detail.json` (b)).
- An overlap there is vacuous.

**Legitimate.**
- The change alters only how the reference is called: a minimum Frobenius term count of 2|ρ| + 8,
  applied through a wrapper, and a precision ladder from 256 to 4096 bits.
- It adds the pass condition "both relative radii < 1e−10".
- arbw itself is unchanged, and its own tail check still guards the result. zsol and the counting
  code are unchanged.

**Stricter.** Yes.
- The pass predicate is a strict strengthening ("overlap" and "both balls tight").
- Under the original predicate the three high-frequency points could only have passed vacuously.
- Now all 24 points are meaningful. The largest relative difference of midpoints is 4.0e−36, and
  arbw needed 512 or 1024 bits at those three points.

**Disclosure.** It is accurate apart from S2 (the test point) and W6 (the failing point named in the
commit message).

**Would the original have caught a real pipeline error?**
- A point-value error in zsol: at 21 of the 24 points, yes; at the three high-frequency points, no.
- An error in zsol's box enclosure: neither version tests it (S4).

## 4. Mathematics (scope B)

### B3. The phase factor

**The factor is identically 1.** Write Z₊ = e^{iω(r*−c₊)}R₊ and Z_c = e^{−iω(r*−c_c)}R_c. Then:
- Z₊Z_c = e^{iω(c_c−c₊)}R₊R_c;
- Z₊Z_c′ − Z₊′Z_c = e^{iω(c_c−c₊)}[R₊R_c′ − R₊′R_c − 2iω(r²/Δ)R₊R_c].

Hence

  Δ(Z₊Z_c′ − Z₊′Z_c) + 2iωr²Z₊Z_c = e^{−iωc₊+iωc_c} Δ(R₊R_c′ − R₊′R_c),

which is r₊r_c times the frozen W_l of FORMULATION §1.2, Jost factor included. The reviewer's sympy
script (`sympy_identities.py` S3) derives this independently of `zeq_identity.py` Z2. A negative
control with the wrong sign of 2iωr² returns False (`sympy_negative_controls.py`).

**The horizon normalisations agree.** zsol computes Z₊ as the unique solution that is analytic at r₊
with value 1.
- It is unique off E₊: the exponents of the s-equation at r_h are 0 and 2isωr_h²/Δ′(r_h), which is
  iω/κ₊ at r₊ and iω/κ_c at r_c (S5).
- This solution is e^{iω(r*−c₊)}R₊. Since A₊ = r₊²/Δ′(r₊) = 1/(2κ₊) (S4),

  e^{iω(r*−c₊)}(r−r₊)^{−iω/(2κ₊)} = exp(iω Σ_{i≠+} A_i log|(r−r_i)/(r₊−r_i)|),

  which is analytic at r₊ with value 1.
- r_c is handled the same way, with A_c = −1/(2κ_c).

**No zeros, poles or branch points are introduced.**
- zsol never uses c₊ or c_c numerically. S1's arbw multiplies by exp(−iωc₊ + iωc_c) computed from
  `Point.c_plus` and `Point.c_cos`.
- C0's agreement to 4e−36 confirms that the two normalisations coincide numerically.
- Had they not, the discrepancy would have been an entire, nonvanishing factor, which moves no zero.

### B4. Analyticity and the true counts

**Analyticity.** W_l is analytic on ℂ∖(E₊ ∪ E_c) (FORMULATION §1.2).
- In zsol's own terms, the exponent-0 recursion divides only by Δ′(r_h)m(m − iω/κ_h), which vanishes
  only on E_h.
- The majorant bounds are locally uniform in ω.
- **No branch points in ω.** The only ω-powers are (r − r₊)^{−iω/(2κ₊)} and (r_c − r)^{−iω/(2κ_c)},
  with fixed positive bases, and these are entire in ω. They cancel in W_zsol anyway (B3).
- **Singular points of the ODE.** These are r_n, r₋, r₊ and r_c only (Δ(0) = Q² ≠ 0). The Taylor paths
  run along the real segments [r₊ + t, r_P] and [r_P, r_c − t], which contain none of them. The
  Frobenius series are used only at |t| well inside their discs of convergence:
  - r₊ − r₋ = 0.097 at r₊;
  - r_c − r₊ = 4.77 at r_c.

  The verified majorant ratio q|t| < 1 enforces both.

**No poles near R′_l.**
- R′_l ⊂ {Im ω ≥ −b}, and `region()` certifies b < κ₊ (b/κ₋ = 0.55, κ₊/κ₋ = 0.8229).
- Every disc box is certified clear of E₊ ∪ E_c; the reviewer's recheck gives a minimum clearance ratio
  of 1.504.
- So the winding number equals the number of zeros with multiplicity.

**l = 0.**
- W₀(0) = 0 exactly, since R ≡ 1 solves the equation at ω = 0.
- With winding 1, ω = 0 is therefore the only zero in R′₀, and it is simple. This gives D4's simplicity
  independently of P2-0's written proof.
- The Z-form gives W₀′(0) = i(r₊² + r_c²)/(r₊r_c) in two lines (sympy S6):
  - Δζ₊′ = i(r² − r₊²) and Δζ_c′ = −i(r² − r_c²);
  - this matches D4 and C1.

### B5. Conjugation symmetry and the axis crossings

**The symmetry.**
- Under ω ↦ −ω̄ the coefficients of the s-equation are conjugated, for real r (Z4; S7).
- Uniqueness then gives Z_h(r; −ω̄) = conj Z_h(r; ω).
- The factored formula then gives W(−ω̄) = conj W(ω) (S7).

**The left half.**
- The reflection maps H (−ib → a−ib → a+ie → ie) onto the left half traversed backwards.
- So Δ_left arg W = Δ_H arg W, and the total change is 2Δ_H.

**The axis crossings.**
- On the imaginary axis W is real.
- At the crossing points −ib and ie it is nonzero, because the end pieces exclude 0.
- So Δ_H ∈ πℤ, and the winding number is Δ_H/π.

**Checked explicitly.**
- The left half was traced in floats for every l; |left − right|/π ≤ 1.7e−13 (`indep_winding.json`).
- It was certified on the full boundary at l = 0, 1, 2 and 5 (`recount.json`).
- The author checked it at l = 2 (C4).

### B6. Rigour of the enclosures

**zsol series** (`zsol.py:66–185`).
- `frob_z`'s recursion has three lags: Δ is quartic with Δ(r_h) = 0, B is cubic and C is linear.
- **Tail bound.** The bound |a_m| ≤ Σ_k g_k|a_{m−k}| uses three facts: n(n−1) ≤ m², n ≤ m, and
  |D_m| ≥ |A₁|m(m − |z|). The resulting g_k(m) = (|A_{k+1}|m² + |B_k|m + |C_{k−1}|)/(|A₁|m(m−|z|))
  decreases in m, so g_k(N) bounds every m ≥ N.
- **Majorant.** q is verified in Arb. K is taken over the last three coefficients, and the geometric
  tails are exact series sums.
- **`taylor_z`.** It uses four lags, with weights 1, 1/(N0+1) and 1/((N0+1)(N0+2)). The latter two are
  conservative for an index p ≥ N0 + 2.
- Both routines raise if q|t| ≥ 1.

**Large |ω|.** `frob_z` does not share arbw's failure mode.
- It sets min_terms ≥ 2|z| + 8 before summing and checks N > 2|z| + 2.
- `propagate` shrinks the Frobenius offset like 1/|z|.
- Any failure is a raise, never a wrong ball.

**The Lohner recursion** (`zsol.py:315–347`). With z′ = B_{k+1}⁻¹(y_{k+1}(ω) − yc_{k+1}),

  z′ = Bn⁻¹[T(ω)yc − T(ω_c)yc] + Bn⁻¹[T(ω_c)yc − yc_{k+1}] + Bn⁻¹T(ω_c)B_k z + Bn⁻¹(T(ω) − T(ω_c))B_k z.

These four terms are `dB`, `add`, `Cm·e` and `G·e`:
- `Bi` encloses the exact inverse of the point frame `Bn`;
- `Tb` ∋ T(ω) for every ω in the box;
- `vb` ∋ T(ω)yc;
- `yn` ∋ T(ω_c)yc.

Further:
- the initial deviation includes the whole Frobenius ball over the box;
- `dev = (e₀ + e₁)(1 + 1e−7)` bounds both components, because the frame columns have unit 2-norm up to
  rounding;
- `wform` is evaluated on the widened Z balls with ω ranging over the box.

So M ≥ sup|W| on the box. The choice of frame affects sharpness only.

**The box bound, tested** (`box_checks.py`).
- **A.** The reviewer's float |W| was sampled on a 9×9 grid and 32 circle points of every one of the
  2351 production boxes. It never exceeds the stored M. M/(max sampled |W|) ranges from 73 to 6e20.
- **B.** On 48 boxes it contains the float W at 16 points and S1's arbw balls at two corners (arbw at
  256 bits, rising to 512 or 1024 bits at the high-frequency corners):
  - half-widths κ₋/16, 1e−4κ₋ and 1e−8κ₋;
  - l = 1, 5, 10 and 19;
  - including the high-frequency corners.
- **Its width scales like ρ.** The relative half-width is 2.3e−7 to 6.1e−6 at ρ = 1e−8κ₋. So it is a
  genuine first-order enclosure, not merely a large ball.

**Cauchy discs** (`part3.py:224–248`).
- The DFT of samples at the exact nodes equals Σ_{m≥0} a_{k+mK}(ρ/2)^{k+mK}.
- W is analytic on the closed box, and M bounds it there, so |a_n| ≤ M/ρⁿ. The aliasing error is then
  ≤ Mx^K/((1 − x^K)ρ^k), with x = 1/2.
- The computed nodes are balls containing the exact nodes.
- The truncation term is M(s/ρ)^K/(1 − s/ρ), where s is the upper bound of |d| over the piece ball.

**Pieces** (`part3.py:256–276`).
- E encloses W on the whole piece. It is an axis-parallel rectangle excluding 0, so it lies in an
  open half-plane through 0.
- The continuous argument change along the piece therefore lies in (−π, π) and equals the principal
  Arg of W(z₁)/W(z₀).
- A ratio ball straddling the branch cut gives a ±π-wide arg, which fails the < π check: a failure,
  not a wrong value.
- Adjacent pieces share the same endpoint balls, so the increments telescope through fixed exact
  points.
- Every popped piece is either accepted or split into two pushed halves, so the accepted pieces
  partition the segment.

**The |W| lower bound (`c7d1bfbc`).**
- The bound max(|E|_lower, dist(Re E, 0), dist(Im E, 0)) holds at every point of E.
- It is a diagnostic only.
- Every stored per-l minimum lies below the reviewer's float minimum on H (§5 table). The C2(d) and C3
  controls show 0.0 because they ran before this change (`3e7fca05`, Arb's acb abs).

**The integer.**
- The true Σ_H/π is an integer (B5).
- An exact rational re-sum of the stored balls lies within 6.5e−4 of 1 (l = 0) or 0 (l ≥ 1)
  (`check_certs.py`).
- So the integer is determined.

### B7. The corollary on S_δ

Let ω ∈ S_δ be a zero of W_l with l ≤ 19.
- **Im ω ≥ 0.** By P2-0, (l, ω) = (0, 0).
- **−(½+δ)κ₋ ≤ Im ω < 0 and |Re ω| ≥ a_l.** P2-2 excludes ω. Its region is closed in both
  inequalities, so the line Im ω = −(½+δ)κ₋ and the vertical lines Re ω = ±a_l are covered.
- **|Re ω| < a_l and −(½+δ)κ₋ ≤ Im ω < 0.** Since −b < −(½+δ)κ₋ and e > 0, ω lies in the open
  rectangle. There the only zero is ω = 0 at l = 0, which contradicts Im ω < 0.

The corners and the real axis are covered by P2-0, P2-2 and the zero-free ∂R′_l. E₊ ∪ E_c does not
meet S_δ (FORMULATION §1.3), so D1 applies throughout.

The statement is correct, but it needs P2-0 for Im ω > 0, so its level is WRITTEN_PROOF (S1).

s3-report §6 does not overclaim relative to P2-3:
- it disclaims P2-3 itself (every l < l₀), P_A and P2-6;
- the row names R′_l, not D7's R_l.

### B8. The l = 0 enclosure

`krawczyk_z` is S1's `winding.krawczyk` with zsol's W:
- the square Z = z₀ ± r(1 + i), with r = 10⁻²⁰|z₀|;
- the Cauchy radius R = 10⁻⁴|z₀|, on a pole-free box;
- M = sup|W₀| on that box, from the zsol box;
- a central difference with h = 10⁻³⁰|z₀| and error Mh²/(R³(1 − (h/R)²));
- a spread of W′ over Z of at most 2Ms/(R − s)², with s = √2 r. This is the S1-review corner fix, and
  it is present.

**Why it is correct.**
- W(z) − W(z₀) = D(z)(z − z₀), where D(z) is the mean of W′ on a segment in Z. D(z) lies in the convex
  rectangle dW.
- So the Newton map sends Z into K(Z) ⊂ int Z, which gives existence.
- |1 − Y·dW| ≤ 5.3e−15 < 1 gives uniqueness in Z.

**What was checked.**
- **Re-run.** Identical to the stored certificate (K1).
- **Independent of the Krawczyk code (K2).** W₀ is real on iℝ. At the two ends of the stored interval,
  Re W₀(iyκ₋) is +2.99e−30 and −2.06e−33, each with a certified sign, by zsol and by S1's arbw. By the
  intermediate value theorem a zero lies inside.
- **The reviewer's float zero (K3)** agrees to 3.7e−15.
- **KAT-6 overlap.** The intervals overlap; the KAT-6 interval, of radius 7.1e−30, contains the S3
  interval.

**β ≤ 0.90938098.** This is correct.
- ω_* lies between −κ₊ = −0.8229κ₋ and −2κ₊, so it is not in E₊ ∪ E_c and is a QNM by D1.
- Hence α ≤ 0.9093809789948054555249337626…κ₋.

## 5. Recomputation and tests (scope C)

### C9. Independent recomputation

| check | independent of | result |
|---|---|---|
| `check_certs.py`: all 20 certificates. Parameters from the exact definitions (mpmath, 80 digits); geometry in exact rationals | all `p2` code | Every certificate tiles H exactly: from (0, −b), via the corners (a_l, −b) and (a_l, e), to (0, e). In each: a_l = the P2-2 threshold; b ≥ 0.55κ₋; 0 < e ≤ κ₋/10; every segment within ρ/2 of its centre; every box clear of E₊ ∪ E_c. Exact re-sums within 6.5e−4 of the integer. 2351 segments |
| `indep_winding.py`: float W from the unfactored system u = R, v = ΔR′ (Frobenius start, scipy DOP853) | arbw, zsol and part3 (different formulation and integrator) | Agrees with zsol at the 24 C0 points to ≤ 2.9e−10. All 2351 segment increments lie inside their stored balls. Full-boundary windings: 1 for l = 0, and below 8e−14 in modulus for l ≥ 1. The left half equals the right half for every l |
| `replay_author.py`: the author's `cauchy_segment`, re-run on every segment of l = 0, 1, 2, 19 (306 segments), with K re-derived | — (reproducibility) | Every increment overlaps. K, the piece counts and M reproduce exactly. Windings 1, 0, 0, 0. 2008 CPU-s |
| `recount.py`: the reviewer's own plan, DFT, error terms, pieces and Arg increments, on the FULL boundary (no symmetry); point values from S1's arbw; only zsol's box M reused | part3's counting layer, zsol's point values, and the symmetry | Winding balls: l = 0: 1.00000 ± 3.3e−7 (58 segments); l = 1: 0 ± 6.5e−7 (98); l = 2: 0 ± 5.7e−7 (140); l = 5: 0 ± 1.05e−6 (265, arbw at 256 bits). All certified; certified min \|W\| on the boundary 0.022, 0.44, 0.55 and 1.07 |
| `box_checks.py` A and B | — | §4 B6: no M violated; tight boxes contain the independent values |
| `c0_and_krawczyk.py`, `c0_detail.py` | — | §3 A2 and §4 B8 |
| `sympy_identities.py` | `zeq_identity.py` | All items true; negative controls false |

Per l (certificates, re-summed exactly; float values from `indep_winding.json`):

| l | a_l/κ₋ | segments | K max | max M/\|W_c\| | Σ_H/π (exact re-sum) | float full-boundary winding | float min \|W\| on H | certified lower bound |
|---|---|---|---|---|---|---|---|---|
| 0 | 5.40 | 21 | 64 | 9.0e11 | 1 ± 6.5e−5 | 1.0000000 | 0.02226 | 0.02184 |
| 1 | 10.31 | 31 | 80 | 2.4e16 | 0 ± 5.1e−5 | 1.4e−14 | 0.5273 | 0.5249 |
| 2 | 15.61 | 41 | 80 | 2.5e17 | 0 ± 1.1e−4 | 1.1e−16 | 0.7807 | 0.748 |
| 3 | 20.86 | 51 | 88 | 1.0e18 | 0 ± 1.0e−4 | 5.8e−16 | 1.019 | 0.9997 |
| 4 | 26.02 | 63 | 88 | 2.2e18 | 0 ± 1.7e−4 | −2.8e−14 | 1.246 | 1.186 |
| 5 | 31.19 | 73 | 88 | 6.4e18 | 0 ± 2.3e−4 | −1.6e−15 | 1.465 | 1.429 |
| 6 | 36.28 | 83 | 88 | 4.2e19 | 0 ± 2.3e−4 | 8.7e−16 | 1.679 | 1.599 |
| 7 | 41.37 | 93 | 88 | 1.9e19 | 0 ± 2.6e−4 | 2.5e−15 | 1.886 | 1.832 |
| 8 | 46.43 | 102 | 96 | 1.5e20 | 0 ± 3.9e−4 | −2.5e−14 | 2.090 | 1.998 |
| 9 | 51.52 | 113 | 88 | 7.2e19 | 0 ± 3.5e−4 | −3.3e−15 | 2.290 | 2.221 |
| 10 | 56.57 | 123 | 96 | 1.8e20 | 0 ± 3.2e−4 | 4.3e−16 | 2.486 | 2.367 |
| 11 | 61.60 | 133 | 96 | 1.7e20 | 0 ± 4.0e−4 | 4.3e−14 | 2.679 | 2.578 |
| 12 | 66.60 | 143 | 96 | 5.4e20 | 0 ± 4.1e−4 | 1.8e−14 | 2.870 | 2.719 |
| 13 | 71.66 | 153 | 96 | 8.9e20 | 0 ± 5.9e−4 | −9.8e−15 | 3.058 | 2.941 |
| 14 | 76.70 | 163 | 96 | 2.4e21 | 0 ± 4.9e−4 | −6.3e−14 | 3.244 | 2.997 |
| 15 | 81.62 | 173 | 96 | 7.5e21 | 0 ± 5.3e−4 | −2.9e−15 | 3.428 | 3.274 |
| 16 | 86.72 | 183 | 96 | 4.0e21 | 0 ± 5.8e−4 | −9.8e−15 | 3.610 | 3.329 |
| 17 | 91.70 | 193 | 96 | 1.1e22 | 0 ± 5.9e−4 | −7.6e−14 | 3.790 | 3.596 |
| 18 | 96.66 | 203 | 96 | 1.4e22 | 0 ± 6.5e−4 | 9.5e−15 | 3.968 | 3.700 |
| 19 | 101.71 | 213 | 96 | 1.1e22 | 0 ± 6.1e−4 | 7.4e−15 | 4.145 | 3.947 |

**Integrity of the stored run.** `compare_author_jsonl.py` read the author's gitignored JSONL working
files in place, read only.
- They contain exactly 2351 ok records, one per planned segment.
- No segment was split, re-run or duplicated.
- Each record equals the corresponding certificate piece.

### C10. Do the controls exercise the pipeline?

| control | what it exercises | assessment |
|---|---|---|
| C0 | zsol point values against an independent code | meaningful at all 24 points after `3e7fca05`; not the box bound |
| C1 (D4) | box M, point values, the DFT coefficient of order 1, against an exact value | good; but M is loose, so it would not detect a moderately small M |
| C2 (four known zeros, full rectangles, no symmetry) | the whole pipeline, as a positive control (a pipeline that always returned 0 would fail), including below −κ₊ and at P_A | good |
| C3 (1/0 pair, edge 0.0103κ₋ above a PS mode) | discrimination at a quarter of a typical disc radius from a zero | good |
| C4 (left half at l = 2) | the symmetry, by explicit computation | good; the reviewer extended it to every l (float) and to l = 0, 1, 2, 5 (rigorous) |
| C5 (S0 float counts) | — | weak, but harmless |
| C6 (Krawczyk against KAT-6) | zsol below −κ₊, between the poles | good; shares the Krawczyk code with KAT-6 (W8) |
| C7 (5% replay) | reproducibility | not committed (S3); the reviewer replayed four l in full |

The gap is the box bound M (S4). It is closed for this review by `box_checks.py`.

### C11. Tests

All runs at nice 10, in the review worktree at `db2a37eb`:

| command | result |
|---|---|
| `pytest problems/P2/tests -q` | 76 passed, 13 skipped, 18.3 s (`replay_pytest_fast.txt`) |
| `P2_SLOW=1 pytest problems/P2/tests -q` | 89 passed, 1 warning (the pre-existing m3 Riccati overflow of S2), 425 s (`replay_pytest_slow.txt`) |

**What `test_part3.py` covers.**
- Fast:
  - the sympy identities;
  - C0 at 3 points;
  - D4;
  - the Krawczyk enclosure against KAT-6;
  - the re-sum and coverage of all 20 certificates;
  - two recomputed segments of l = 0.
- Slow:
  - the C2(c) and C3 rectangles;
  - a 5% replay for l ≤ 5.

**What it does not cover.** No test checks the box bound against an independent value, and no test
recomputes anything at l ≥ 6 (S4).

`git status` showed no change to any tracked file after any run.

## 6. Compute

| item | CPU |
|---|---|
| pytest fast, twice, and slow | 0.13 core-h |
| `replay_author.py` (306 segments) | 0.56 core-h |
| `recount.py` (l = 0, 1, 2, 5) | 0.48 core-h |
| `indep_winding.py`, `box_checks.py` A and B, `check_certs.py`, `c0_and_krawczyk.py`, `c0_detail.py`, sympy, timing probes | 0.07 core-h |
| **total** | **1.24 core-h** (cap 3); at most 3 concurrent workers, nice 10 |

## 7. The ledger

**Can P2-3B be CERTIFIED?** Yes, for the counts.
- The certificate is under `certificates/part3/`.
- The decisive steps are machine-checked: the enclosures, the pole clearances, the coverage of H, the
  sums, and the identities (sympy).
- The written steps are routine:
  - the argument principle;
  - the conjugation symmetry (two lines plus Z4);
  - analyticity off E₊ ∪ E_c;
  - the identification of the normalised exponent-0 solution with the frozen Jost solution;
  - the validity of standard bounds (majorants, Cauchy estimates, the Lohner recursion), which this
    review checked line by line.

**The corollary is at P2-0's level.**
- The strip corollary is a consequence that needs P2-0 for Im ω > 0, and P2-0 is WRITTEN_PROOF.
- Precedent: the P2-2 row says "at P2-0's level", and P2-1w says "The Im ω > 0 part is at P2-0's
  level". The row must do the same.
- WRITTEN_PROOF for the whole row is not needed. Nothing non-routine is written in S3 itself.

**P2-3B (new row; replaces the text of s3-report §7):**

| P2-3B | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20) and for each l = 0, 1, …, 19, let R′_l = [−a_l, a_l] × [−b, e], with a_l the P2-2 threshold (`W_l_exact` of certificates/lg/lg_PB_l{l}.json; a_l/κ₋ from 5.40 at l = 0 to 101.71 at l = 19), b = 6280879/2²⁸ (b/κ₋ = 0.55000001 ≥ ½ + δ) and e = 4567911/2³⁰ (e/κ₋ = 0.09999998 ≤ 1/10). Then W_l is analytic on a neighbourhood of R′_l and has no zero on ∂R′_l, and its winding number around ∂R′_l is 1 for l = 0 and 0 for 1 ≤ l ≤ 19. So R′₀ contains exactly one zero of W₀, the simple zero at ω = 0 (D4), and R′_l contains none for 1 ≤ l ≤ 19. With P2-0 (Im ω > 0) and P2-2 (abs(Re ω) ≥ a_l), at P2-0's level: for every l ≤ 19, S_δ contains no QNM of degree l other than ω = 0 (l = 0). This is P2-3 for l ≤ 19 only, on R′_l rather than D7's R_l: P2-3 itself needs every l < l₀, and l₀ comes from S4 | CERTIFIED for the counts on R′_l (Arb, 128-bit). Cauchy-interpolation discs on the right half of ∂R′_l (conjugation symmetry), 2351 segments; every winding ball lies within 6.5e−4 of its integer. W comes from the phase-factored Jost solutions (zsol): an exact identity, checked in sympy, makes it equal to the frozen W_l. Cross-checked against the S1 code at 24 points (relative difference ≤ 4e−36, both balls below 1e−10). Controls: four known zeros counted (P_B l = 0, 1, 2; P_A l = 1); a 1/0 pair of rectangles sharing an edge 0.0103κ₋ above a PS mode; the left half computed explicitly at l = 2. Written steps, all routine: the argument principle, the conjugation symmetry, analyticity of W_l off E₊ ∪ E_c (FORMULATION §1.2), and the identification of the zsol solutions with the Jost solutions. The corollary on S_δ is at P2-0's level (WRITTEN_PROOF) for Im ω > 0. Independent review 2026-10-17: proven as stated. Every certificate was re-summed exactly, with its tiling, corners and pole clearance re-checked from independently computed parameters. An independent float code reproduces all 2351 increments and traces the left half for every l. l = 0, 1, 2, 19 were replayed in full; l = 0, 1, 2, 5 were recounted rigorously on the full boundary with the reviewer's counting layer and S1 point values; the box bound M was tested on all 2351 discs | problems/P2/src/p2/zsol.py, part3.py, zeq_identity.py; problems/P2/certificates/part3/; problems/P2/notes/s3-preregistration.md, s3-report.md; problems/P2/tests/test_part3.py; docs/review-2026-10-17-p2-s3-adversarial.md; problems/P2/certificates/review_s3/ | 2026-10-17 |

**P2-K6 amendment (accept, corrected).** Two changes to the existing row:
- **Statement:** "radius 6e−30" becomes "radius 7.1e−30". This is the value of the KAT-6 certificate
  since `9f1b8e6f` (W9).
- **Evidence:** append "Re-derived in S3 with a different Wronskian code (the phase-factored solutions
  of problems/P2/src/p2/zsol.py; the Krawczyk construction, its parameters and the Point module are
  KAT-6's): K(Z) ⊂ int Z for the square Z of half-width 10⁻²⁰·abs(ω_*), so W₀ has a unique zero in Z,
  enclosed in −0.909380978994805455524933762622i κ₋ ± 3.2e−31 (real part ± 6.6e−35), contraction
  5.3e−15 (problems/P2/certificates/part3/krawczyk_PB_l0.json). Independent review 2026-10-17: re-run
  identical; a certified sign change of W₀(iyκ₋) across that interval with both the S3 and the S1
  codes".

The level is unchanged: CERTIFIED for the zeros.

**Optional, for the owner.** The one-sided bound β ≤ 0.90938098 follows at P2-K6's level. It is the
upper half of the stretch goal P2-6. It can be added to the P2-K6 statement, or left in s3-report §4
until P2-6 is attempted.

## 8. Before the S3 branch is merged into main

1. **Ledger text (S1).** Use the §7 text for P2-3B and for the P2-K6 amendment. Either correct
   s3-report §7 and FORMULATION v0.5 ("at P2-0's level") to match, or mark them superseded by this
   review.
2. **The test change (S2).** Restore (7, 40.0, −0.55) in `test_zsol_matches_arbw`, or disclose the
   change in s3-report §3.
3. **Missing outputs (S3).** Commit the C7 replay output and the l = 0 left-half sum. Alternatively,
   cite `review_s3/`, and drop "l = 0" from the symmetry claim.
4. **Wording (W1–W9)**, as convenient. Three of them touch the ledger: W3 (the clash between β and the
   strip depth), W8 (what "independent" means) and W9 (the K6 radius).
5. **Optional (S4).** Add a tight-box containment test against arbw at an l = 19 corner, and a replay of
   one or two l = 19 segments.

Under the repository's rules the merge into main must be run by the user.

## 9. Follow-up check of the author's fixes (`e874b5f8`, `b19065b6`), 2026-10-17

**Setup.**
- `p2/s3-review` was fast-forwarded to `p2/s3-part-iii` at `b19065b6`: the author's two commits on top
  of this review's `123277fb`.
- Against `123277fb`, only four files changed: FORMULATION.md, s3-report.md, `part3.py` and
  `test_part3.py`.
- Two records were added: `certificates/part3/replay_PB.json` and `symmetry_PB_l0.json`.
- This review's files, CLAIMS.md and the README were not touched.
- New outputs are in `review_s3/`: `followup_checks.py` and its `.json`/stdout, `followup_verify_out.json`,
  and `followup_pytest_{fast,slow}.txt`.

**Verdict: PASS, merge-ready.** Every finding is resolved as intended, and nothing changes a certified
result. The ledger text of §7 stays correct; an optional update for the new evidence is below.

| item | result |
|---|---|
| S1 (the corollary's level) | **PASS.** The s3-report verdict, §1 and §7, and FORMULATION v0.5 now say the counts are CERTIFIED and the strip corollary is at P2-0's level. They name the band 0 < Im ω ≤ e with \|Re ω\| ≥ a_l, and cite analyticity under FORMULATION §1.2. s3-report §7 is the review's §7 verbatim: the P2-3B row and the P2-K6 block compared string for string. The v0.5 entry records that its wording was corrected after the review, so the edit is not silent |
| S2 (the test point) | **PASS.** `test_zsol_matches_arbw` tests both (7, 40.0, −0.55) and (7, 30.0, −0.55), and both pass. s3-report §3 discloses the change and withdraws "the fix touched only the control" |
| S3 (missing outputs) | **PASS.** `replay_PB.json` and `symmetry_PB_l0.json` are what they claim (F1–F3 below). The C4 and C7 rows of s3-report §3 cite them |
| S4 (box-enclosure tests) | **PASS.** Fast: the zsol box ball at the corner (a₁₉, −b), half-width 1e−4κ₋, must contain S1's arbw enclosures at two box corners. P2_SLOW: the same at κ₋/16 with four corners, and two stored l = 19 segments recomputed. All pass. The docstring's "about 60×" reproduces as 56×: the ball's relative half-width is 5.83e−3, and the corner values lie 1.039e−4 from the centre (F4) |
| W1–W8 | **PASS.** W1 (the 10¹⁸ figure) is corrected in §2. W2 (the C0 points) is stated in full. W3 is fixed in FORMULATION, plus a correct erratum in s3-report §1 for the pre-registration's two uses of β and its "S1, D1′" citation; I checked all four quoted places against the pre-registration. W4 (FORMULATION header) is fixed. W5 (format) is noted in §2 and §4. W6 (failing point) is in §3. W7 (pilot compute) is in §5. W8 (what "independent" means) is in §4 and the verdict |
| W9 (P2-K6 radius) | Left to the coordinator, as intended; it is part of §7's K6 amendment |
| New code paths | **PASS.** `cauchy_segment`, `plan`, `region`, `assemble`, `run`, `controls`, zsol, arbw, params, m3, lg_bounds and winding are unchanged; only `part3.py` changed in `src/`. The certificates read nothing new. The changes are: `verify` (records results, `--workers`, `--out` through `_replay_task`, which calls the unchanged `cauchy_segment` with the stored K at the same 128-bit precision); `symmetry_check` (`--lsym`, through the unchanged `left_half_sum`); argument parsing. None can alter a winding certificate. The new path was re-run with 2 workers for l = 0, 1, 5, 10, 19 and reproduces `replay_PB.json` string for string (F2) |
| Pre-registration | **PASS.** `notes/s3-preregistration.md` is byte-identical to `89a2336b`: same blob `4e46db04`, sha256 `053ceeaa…` |
| Tests | **PASS.** Fast: 77 passed, 15 skipped, 25 s. P2_SLOW=1: 92 passed, 1 pre-existing warning (the S2 m3 overflow), 452 s. The three new part3 tests ran |

**Spot-replays of the new records** (`followup_checks.py`):
- **F1, `replay_PB.json`.**
  - It holds 132 segments, exactly the deterministic sample of `part3.verify`.
  - Every stored ball and K equals its certificate piece.
  - The 20 segments at l = 0, 1, 2 and 19 equal this review's full replay string for string.
  - All 132 recomputed balls contain the reviewer's independent float increments, at most 0.87 radii
    from the printed midpoint.
  - It is stamped with the clean commit `e874b5f8`. That commit added the code, and the outputs were
    committed in `b19065b6`.
- **F2, the new verify path.** Re-run with 2 workers and `--out`, it produces identical recomputed
  balls and no problems on its 28 segments.
- **F3, `symmetry_PB_l0.json`.** The l = 0 left half, recomputed with the unchanged `left_half_sum`,
  gives exactly the committed "[3.142 +/- 5.00e-4]", and its quotient by π, [1.000 ± 2.89e−4],
  contains 1.
- **F4, the fast box test's numbers.** These are the numbers quoted in the S4 row above.

**Compute for the follow-up:** 0.19 core-h; at most 3 processes, nice 10. That is 197 CPU-s of checks
and 474 s of tests. The review total is now 1.43 core-h.

**The ledger text.** The §7 P2-3B row and the P2-K6 amendment remain correct as written.

Now that the C7 record and the l = 0 left half are committed, the controls sentence of P2-3B can name
them. This is optional. The paste-ready row with that change is below; it differs from §7 only in the
"Controls" sentence and in the review sentence:

| P2-3B | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20) and for each l = 0, 1, …, 19, let R′_l = [−a_l, a_l] × [−b, e], with a_l the P2-2 threshold (`W_l_exact` of certificates/lg/lg_PB_l{l}.json; a_l/κ₋ from 5.40 at l = 0 to 101.71 at l = 19), b = 6280879/2²⁸ (b/κ₋ = 0.55000001 ≥ ½ + δ) and e = 4567911/2³⁰ (e/κ₋ = 0.09999998 ≤ 1/10). Then W_l is analytic on a neighbourhood of R′_l and has no zero on ∂R′_l, and its winding number around ∂R′_l is 1 for l = 0 and 0 for 1 ≤ l ≤ 19. So R′₀ contains exactly one zero of W₀, the simple zero at ω = 0 (D4), and R′_l contains none for 1 ≤ l ≤ 19. With P2-0 (Im ω > 0) and P2-2 (abs(Re ω) ≥ a_l), at P2-0's level: for every l ≤ 19, S_δ contains no QNM of degree l other than ω = 0 (l = 0). This is P2-3 for l ≤ 19 only, on R′_l rather than D7's R_l: P2-3 itself needs every l < l₀, and l₀ comes from S4 | CERTIFIED for the counts on R′_l (Arb, 128-bit). Cauchy-interpolation discs on the right half of ∂R′_l (conjugation symmetry), 2351 segments; every winding ball lies within 6.5e−4 of its integer. W comes from the phase-factored Jost solutions (zsol): an exact identity, checked in sympy, makes it equal to the frozen W_l. Cross-checked against the S1 code at 24 points (relative difference ≤ 4e−36, both balls below 1e−10). Controls: four known zeros counted (P_B l = 0, 1, 2; P_A l = 1); a 1/0 pair of rectangles sharing an edge 0.0103κ₋ above a PS mode; the left half computed explicitly at l = 0 and 2; a deterministic 5% sample (132 segments) recomputed (certificates/part3/replay_PB.json); the box enclosure tested against the S1 code at the l = 19 corner (tests/test_part3.py). Written steps, all routine: the argument principle, the conjugation symmetry, analyticity of W_l off E₊ ∪ E_c (FORMULATION §1.2), and the identification of the zsol solutions with the Jost solutions. The corollary on S_δ is at P2-0's level (WRITTEN_PROOF) for Im ω > 0. Independent review 2026-10-17, with a follow-up check of the fixes: proven as stated. Every certificate was re-summed exactly, with its tiling, corners and pole clearance re-checked from independently computed parameters. An independent float code reproduces all 2351 increments and traces the left half for every l. l = 0, 1, 2, 19 were replayed in full; l = 0, 1, 2, 5 were recounted rigorously on the full boundary with the reviewer's counting layer and S1 point values; the box bound M was tested on all 2351 discs | problems/P2/src/p2/zsol.py, part3.py, zeq_identity.py; problems/P2/certificates/part3/; problems/P2/notes/s3-preregistration.md, s3-report.md; problems/P2/tests/test_part3.py; docs/review-2026-10-17-p2-s3-adversarial.md; problems/P2/certificates/review_s3/ | 2026-10-17 |

The P2-K6 amendment is unchanged from §7. The owner's open choice also stands: β ≤ 0.90938098 can go
into P2-K6 now, or wait for P2-6.

Under the repository's rules the merge into main must be run by the user.
