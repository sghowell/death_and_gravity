# P2 S2 report: part-(i) feasibility probe at P_B

Status: S2 exit, 2026-10-17. Revised the same day after the independent review
(`docs/review-2026-10-17-p2-s2-adversarial.md`, commit `03f949f5`; nothing
blocking). Branch `p2/s2-part-i`, from main at `9f1b8e6f`.

**Verdict: GO for S3** under the pre-registered rule (`s2-preregistration.md` §6).
The kill criterion is not met: l = 100 certifies, with margin 0.56 μ_HO. The
rule's projected l₀ = 2 is only the formal output of the pre-registered fit. That
fit is mis-specified (§3), and the verdict does not depend on it.

What the revision changed (the review's should-fix items):
- §6 carries the review's ledger text: P2-1w at WRITTEN_PROOF, and the P2-2
  amendment with its evidence path (S1, S8).
- §3 replaces the loss attribution and qualifies the projection (S2).
- §2 discloses how A0 and A3 were implemented differently from the
  pre-registration, and what the design-study figure measured (S3). It also
  corrects the A2 wording and adds the near-mode control A2′ (S4).
- `s2-m3.md` gains the Frobenius radii, the x-representation and the tail facts
  (S5), plus the W1 and W2 wording fixes.
- The float round trip in `tail_bound` is fixed (S6, `85fad9f1`). All m3
  certificates were regenerated at `85fad9f1` and re-verified (§2).
- New tests replay the l = 100 P2-2 certificate and re-run the near-mode control
  (S7).

## 0. Decision record (step 0)

- Problems document revision 1.2 is approved as proposed in `s0-gate.md` §7, with
  one change: the two-sided β enclosure is a stretch goal, not part of the entry
  point.
- The .tex carries the revision-log entry and the P2 Status, Entry point and
  bibliography changes. The PDF was rebuilt with tectonic. The one new TeX
  warning is a mild underfull line in a bibliography entry.
- FORMULATION v0.3 records the decision and s0-gate.md has a one-line decision
  record (commit `3aa4af27`).
- No CLAIMS.md row depends on the decision: the rows P2-par, P2-0, P2-D1′,
  P2-2 and P2-K6 are stated at both points.

## 1. Route and design

The route is **M3**: complex rotation x = x_P + e^{iπ/4}s of the tortoise
coordinate through the photon sphere. On the rotated line a QNM becomes a decaying
solution of −u″ + e^{iπ/2}(V_l − ω²)u = 0. Exclusion follows from

  Re⟨e^{−iφ}(−u″ + Qu), u⟩ ≥ ∫ c|u′|² + q|u|² ≥ μ‖u‖² > 0,

where:
- φ is a variable phase, c = cos φ − η|φ′| and q = Re(e^{−iφ}Q) − |φ′|/(4η);
- μ comes from a ground-state (Allegretto–Piepenbrink) pointwise inequality
  E = q + c′G′ − cG′² + cG″ ≥ μ, checked in Arb on s-cells.

Proof: `notes/s2-m3.md`. At the barrier top the rotation turns the inverted
oscillator into a harmonic one, so the trapping gap γ₀ is built in: the ideal
margin is 1 − β/γ₀ = 0.673 of μ_HO = 2(l+½)Ω_Pγ₀.

The brief named M2 as the leading route and M3 as the alternate. M3 was run
first, and the choice was recorded in the pre-registration before production
(§2 there). M2's O(h)-quality centre needs the complex-coalescence Weber map of
route M1, which is open analysis; M3 needs only elementary functions. M2 was not
implemented.

## 2. Certified windows (P_B; A0–A4 all pass)

The window is 𝒲_l = [ω_l − 8γ₀, ω_l + 6γ₀] × [−β, 0], with β = 0.55κ₋. In
scaled units it runs from 8 oscillator quanta below the barrier top to 6 above.

Certificates are in `certificates/m3/`, at code commit `85fad9f1`. They were first
produced at `c4159e7f`, then regenerated after the review's S6 fix, with the new
A2′ control.
- μ_cert, every verdict, every ground-state node and every tail bound are
  identical to the `c4159e7f` values.
- The only other changes are the commit stamp and the float design eigenvalue,
  which differs in its last digits (about 1e-14 relative; it sets only the
  bisection target).
- The float design step is not bit-reproducible: a re-run of the A2′ design gave
  ground-state nodes differing in their last bits. Every replay (A4, and the A2′
  test) therefore uses the stored nodes.

| l | Re ω/κ₋ | μ_cert | μ_HO | μ_cert/μ_HO | float λ at low corner/μ_HO | min c | a_l/κ₋ ≤ top |
|---|---|---|---|---|---|---|---|
| 20 | [87.370, 110.935] | 0.32455 | 0.61436 | **0.528** | 0.661 | 0.119 | 106.64 ≤ 110.94 |
| 50 | [234.935, 258.500] | 0.85991 | 1.51342 | **0.568** | 0.681 | 0.120 | 255.51 ≤ 258.50 |
| 100 | [480.877, 504.442] | 1.68539 | 3.01185 | **0.560** | 0.679 | 0.121 | 503.43 ≤ 504.44 |

What each acceptance criterion found:

- **A0 (geometry, independent of l).**
  - The validated contour satisfies Re r*(r(s)) = x_P + s/√2 within its balls.
  - The sectors between the real axis and the contour are singularity-free. The
    proof uses the Re r* identity and Σ A_i = 0: the lower-left sector is free
    because Re r*(0) = −2.29664 > x_P = −2.43797, and the two special vertical
    segments of the upper-right sector were validated.
  - Tail bounds: |r − r₊| ≤ 0.00956 for s ≤ −64 and |r − r_c| ≤ 0.351 for s ≥ 32.
  - **Deviation from the pre-registration** (`s2-preregistration.md` §4 item 1
    and §5 A0, lines 68–71 and 90–94).
    - A0 was registered as a fan of validated rays: ray spacing below half the
      certified analyticity radius, and certified horizon funnels.
    - The implementation uses a different argument: the Re r* identity with
      Σ A_i = 0, which confines any singularity to Re x = −2.29664 or Re x = 0,
      plus two validated vertical segments (`s2-m3.md` Lemma 1).
    - The argument is valid and simpler, and the review checked it. But it is
      not the registered check, and the first version of this report presented
      it as if it were.
    - The unused fan code (`Ray.analyticity`) was removed at `85fad9f1`.
- **A1.** E ≥ μ_cert > 0 on all cells and on both tails. c ≥ 0. The Jost decay
  condition holds.
- **A2 (negative control; weak).** Boxes of half-width γ₀/50 around each l's PS
  fundamental mode (imaginary parts −1.6836, −1.6833 and −1.6832 in κ₋ units) do
  not certify.
  - Branch (a) uses the window's ground state: μ = −0.068, −0.097 and −0.285,
    i.e. −0.110, −0.064 and −0.095 μ_HO.
  - Branch (b) designs a ground state at the box. The float margin there is
    ≤ 0, so the Riccati sweep overflows (`design_ground_state`, m3.py:212) and
    the checker never runs. It counts as "not certified" by convention; the
    record now says `checker_run: false`.
  - So A2 catches only gross errors. The window's ground state was designed
    1.1κ₋ above the mode, and a checker that overestimated the margin by less
    than about 0.06–0.11 μ_HO would still return "not certified" on these boxes.
    The first version's "the method therefore cannot certify a true QNM"
    overstated one negative test.
- **A2′ (near-mode control; added after the review, its §5).** The same checker,
  with a ground state designed at Re ω = Re ω_PS, Im ω = −0.88γ₀, where the float
  margin is still positive. μ/μ_HO; boxes of real half-width γ₀/50 about
  Re ω_PS:

  | l | float λ at the design point | IN (±γ₀/50 about the mode) | STRAD, Im ∈ [−1.05, −0.95]γ₀ | OUT, Im ∈ [−0.90, −0.86]γ₀ (mode-free) |
  |---|---|---|---|---|
  | 20 | +0.101 | −0.102, not certified | −0.131, not certified | **+0.017, certified** |
  | 50 | +0.116 | −0.113, not certified | −0.143, not certified | **+0.007, certified** |
  | 100 | +0.119 | −0.142, not certified | −0.172, not certified | −0.022, not certified |

  - Every box containing the mode fails, with a near-optimal design, by a clear
    margin. The numbers reproduce the review's run
    (`certificates/review_s2/sharp_negative_control.py`).
  - At l = 20 and 50 the mode-free box 0.1γ₀ above the mode certifies. So the
    checker's verdict switches within about 0.1γ₀ of a true QNM, and A2's
    failure is not an artefact of a poor design.
  - At l = 100 OUT does not certify: the spline loss of §3 (0.07 μ_HO) exceeds
    what is left of the float margin.
  - The float margin vanishes at the modes: −0.020, −0.004 and −0.001 μ_HO at
    the PS mode for l = 20, 50 and 100.
  - The mode location comes from the float S1 solver, so A2′ is a control, not a
    certificate. Its cost was 203, 129 and 102 s of CPU.
- **A3.** μ_cert ≤ the float margin, for every l.
  - **Deviation from the pre-registration** (§5 A3, lines 100–101). A3 was
    registered against "the float best margin λ_min(l)" of the window.
    - The code computes the float λ only at the design point, the low-left
      corner (a₁, −β) (`certify_l`; the comparison is in `main_run`).
    - The worst corner is the low-right corner (a₂, −β): float 0.648, 0.661 and
      0.667 μ_HO, against 0.661, 0.681 and 0.679 at the low-left corner.
    - A3 still passes: μ_cert lies below the float λ at all four corners and
      along the bottom edge (review §2, S3).
  - **The design-study figures** in the pre-registration (§2: 0.58, 0.65 and
    0.66 for l = 20, 50 and 100, labelled "worst corner of the window") are not
    each corner's own float λ, which is 0.648, 0.661 and 0.667.
    - They are the minimum of E over the corners for one common ground state,
      designed at the low-left corner. Re-running the scratch study gives 0.584,
      0.650 and 0.664.
    - The certificate uses the same construction (one ground state for the whole
      ω² hull), and so loses more than the per-corner λ suggests (§3).
- **A4.** `python -m p2.m3 verify` replays every verdict and μ_cert identically,
  and replays the l = 100 P2-2 criterion (61 s; re-run at `85fad9f1` after the
  regeneration).
- **P2-2 overlap.** a_l ≤ window top at every l, with margins 4.296, 2.994 and
  1.008κ₋. For l = 100 the threshold was computed in S2 with the S1 code
  (`certificates/m3/lg_PB_l100.json`, criterion 1.881 < 2). The margin there is
  thin but certified: an exact dyadic against a ball.

So at l = 20, 50 and 100 the whole set {Re ω ≥ ω_l − 8γ₀, −(½+δ)κ₋ ≤ Im ω ≤ 0}
is QNM-free at P_B: the window covers the part below a_l, and P2-2 the rest.

## 3. Error constant and projected l₀ (pre-registered fit)

| fit of m(l) = μ_cert/μ_HO | m_∞ | C | max residual | projected l₀ (formal) |
|---|---|---|---|---|
| m_∞ − C/(l+½) | 0.577 | 0.94 | 0.010 | 2 |
| m_∞ − C₂/(l+½)² | 0.567 | 16.0 | 0.007 | 5 |

**The projected l₀ = 2 is the formal output of the pre-registered fit.** The
fit is mis-specified, and the projection has no predictive content.

Where the margin goes (float; review §4, plus a scan of the node ends; units of
μ_HO):

| l | float λ at the low corner | × 0.95 (design slack) | min E at interval midpoints (one ground state for the ω² hull) | min E at node ends (spline sawtooth) | μ_cert |
|---|---|---|---|---|---|
| 20 | 0.661 | 0.628 | 0.551 (−0.077) | 0.528 (−0.023) | 0.528 |
| 50 | 0.681 | 0.647 | 0.616 (−0.031) | 0.568 (−0.048) | 0.568 |
| 100 | 0.679 | 0.645 | 0.630 (−0.015) | 0.560 (−0.070) | 0.560 |

- The loss has four parts:
  - **the design slack**: the target is 0.95 λ_float, about 0.033 μ_HO;
  - **one ground state for the whole ω² hull**: about 1/l (0.077, 0.031,
    0.015). Part of it is the hull rectangle being larger than the window's ω²
    image (review W1);
  - **the spline sawtooth**: cG″ jumps at the nodes of the piecewise-linear G′,
    and at the fixed node spacing 0.01 the loss grows about like √l (0.023,
    0.048, 0.070);
  - **the mean-value cells**: < 1e-4 μ_HO.
- The first version attributed 0.11–0.13 μ_HO to "the spline ground state and
  mean-value cells". That was half right: the cells cost nothing measurable.
- The 1/l and √l parts have opposite l-dependence, which explains the
  non-monotone m(l) = 0.528, 0.568, 0.560. A curve m_∞ − C/(l+½) through three
  such points predicts nothing. With the design fixed, m(l) eventually falls
  like −√l instead of tending to 0.577.

What the three certified points do and do not imply:
- **They do imply** that the barrier-top window is QNM-free at exactly l = 20,
  50 and 100, and, with P2-2 and P2-0, that {Re ω ≥ ω_l − 8γ₀} ∩ S_δ is QNM-free
  at these l.
- **They imply nothing at any other l**: not between 20 and 100, not above 100,
  and not below 20, where nothing was computed. The first version's "the
  barrier-top regime poses no obstruction at any l ≥ 20, and on this fit at any
  l ≥ 2" is withdrawn.
- **The method is nearly exact; the design is not.** At the l-th PS mode, a true
  QNM where the ideal margin is 0, the float margin is −0.020, −0.004 and
  −0.001 μ_HO for l = 20, 50 and 100, i.e. O(l⁻²). The losses above are design
  losses, not losses of the method.
- **What the √l term means for S4.** A fixed design cannot be uniform in l. At
  node spacing 0.01 the sawtooth alone would use up the ≈ 0.63 μ_HO left after
  the other losses at l of order 10⁴ (a rough extrapolation of three points).
  S4 must scale the design with l:
  - node spacing ∝ (l+½)^{−1/2}, i.e. fixed in σ = s(l+½)^{1/2};
  - window width ∝ the oscillator quantum.

  With that scaling the hull and sawtooth losses become l-independent at
  leading order, and the closing box at h = 1/(l+½) = 0 is the
  harmonic-oscillator ground-state bound. Whether the O(h) corrections stay
  below the margin on a whole h-box is S4's question; S2 does not answer it.
- The low-l range of S3 is therefore not chosen from l₀ = 2
  (`s3-preregistration.md`).

## 4. What S2 does not establish (carried to S3/S4)

- **Uniformity in l.** Only l = 20, 50 and 100 are certified. Part (i) needs the
  inequality for all l ≥ l₀ via the scaled variable σ = s(l+½)^{1/2} and closing
  boxes at h = 0, with the l-scaled design of §3. There the condition reduces to
  the harmonic-oscillator ground-state bound with margin 0.673 minus design
  losses; its h-structure matches the P4 "closing box" pattern.
- **Sub-barrier regime** (Re ω < ω_l − 8γ₀, down to the low-frequency region)
  and the **low-frequency regime** near the horizon throat. The probe did not
  address them. The float design study suggests the same contour keeps a
  positive margin at least to 8 quanta below the top. Lower energies may need a
  different phase profile, or the forbidden-region Riccati argument of s0-gate
  §4.
- **P_A** has not been probed. Its strip reaches 0.667γ₀, so the ideal M3 margin
  shrinks to 1 − 0.667 = 0.333 μ_HO. The design losses seen here (0.11–0.13 μ_HO
  at these l, with the fixed design) would leave about 0.2. That is still
  positive, but tight.

## 5. Cost

All runs used at most 2 cores at nice 10; one P4 probe shared the machine.

| item | CPU |
|---|---|
| geometry certificate | 0.6 s |
| window certificates | 59 s (l=20, deep bisection), 0.9 s (l=50), 0.8 s (l=100) |
| A2 negative controls | 55 / 44 / 43 s |
| A2′ near-mode controls (after the review) | 203 / 129 / 102 s |
| a_100 | 0.9 s |
| replay | 61 s |
| tests | 11.5 s fast (51 pass, 6 skipped); 316 s with `P2_SLOW=1` (57 pass) |
| float design study (scratch), and the re-run of §2 A3 | ≈ 3 min, ≈ 1 min |
| **S2 production total** (first run plus regeneration) | ≈ 15 min |

## 6. Ledger text (from the review, §10; for the coordinator's bundle, not applied)

This replaces the first version's suggested text, which put two levels in one
row.

**P2-1w (new row):**

| P2-1w | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20) and for l = 20, 50, 100: W_l(ω) ≠ 0 for every ω with (l+½)Ω_P − 8γ₀ ≤ Re ω ≤ (l+½)Ω_P + 6γ₀ and −(½+δ)κ₋ ≤ Im ω ≤ 0 (the barrier-top window; Re ω/κ₋ ∈ [87.37, 110.94], [234.94, 258.50], [480.88, 504.44]). Certified coercivity margins μ/μ_HO = 0.528, 0.568, 0.560, with μ_HO = 2(l+½)Ω_Pγ₀. With P2-2 (a_l below the window top: a_20, a_50 from the S1 certificates; a_100 = 503.43κ₋ from certificates/m3/lg_PB_l100.json) and P2-0, no zero of W_l lies in {Re ω ≥ (l+½)Ω_P − 8γ₀} ∩ S_δ at these l. Three discrete l only, so not part (i): there is no uniformity in l, and the sub-barrier and low-frequency regimes are open | WRITTEN_PROOF. Machine-checked (Arb, 128-bit): the ground-state inequality E ≥ μ on s-cells covering [−64, 32] and on both tails, for the whole ω² hull of each window; c ≥ 0.119; the contour enclosure; x_P < Re r*(0); the two special vertical segments; the tail discs; the decay condition Re ω + Im ω > 0. Written (notes/s2-m3.md): the sector analyticity of r(x) (Lemma 1); the continuation and decay of the Jost solutions on the rotated line (Lemma 3); the coercivity and ground-state integrations by parts (§3, Lemma 4). Independent review 2026-10-17: proven as stated. The margins were reproduced by an independent Arb implementation, the identities were checked in sympy, and the complex-scaled float spectrum reproduces the PS modes and has no eigenvalue in the windows. The Im ω > 0 part is at P2-0's level | problems/P2/src/p2/m3.py, m3_contour.py; problems/P2/certificates/m3/; problems/P2/notes/s2-m3.md, s2-preregistration.md, s2-report.md; problems/P2/tests/test_m3.py; docs/review-2026-10-17-p2-s2-adversarial.md; problems/P2/certificates/review_s2/ | 2026-10-17 |

**P2-2 amendment.** Two changes to the existing row; the level is unchanged
(CERTIFIED on the lower strip):
- **Statement:** after "for every 0 ≤ l ≤ 50", insert "and, at P_B, also l = 100
  (a_100 = 503.43κ₋)".
- **Evidence:** append "problems/P2/certificates/m3/lg_PB_l100.json (S1 code,
  unchanged since a78a3057; run at c4159e7f; criterion 1.881 < 2; replayed by
  `python -m p2.m3 verify`, not by `p2.lg_bounds verify`); independent review
  2026-10-17: replayed; actual deviation 0.012 at l = 100".

Since the review, `test_m3.py::test_lg_l100_replay` also replays this
certificate under pytest.

## 7. Commits

| commit | subject |
|---|---|
| `3aa4af27` | P2: record the approval of problems-document revision 1.2 |
| `0d115788` | P2: pre-register the S2 part-(i) feasibility probe |
| `c4159e7f` | P2: implement the S2 part-(i) route M3 at the barrier top |
| `f950b34a` | P2: keep the KAT-6 re-run test from overwriting its certificate |
| `6b933896` | P2: S2 barrier-top certificates, proof note and report |
| `03f949f5` | P2: independent adversarial review of S2 (route M3 barrier-top windows) |
| `85fad9f1` | P2: apply the S2 review's code fixes (tail bound, controls, remnants) |
| (this commit) | P2: apply the S2 review's should-fixes to the notes, regenerate the m3 certificates, add tests |

## 8. Reproduction (cwd `problems/P2/src`)

```
python -m p2.m3 run 20 50 100   # geometry, windows, A2/A2′ controls, overlap, fit (~10 min)
python -m p2.m3 verify          # A4 replay (~1 min)
python -m pytest problems/P2/tests -q                 # from the repo root (~12 s)
P2_SLOW=1 python -m pytest problems/P2/tests -q       # ~5.5 min
```
