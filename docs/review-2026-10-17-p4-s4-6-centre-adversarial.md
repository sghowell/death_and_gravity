# Independent adversarial review: P4 S4-6, the centre of the Theorem C strip

Date: 2026-10-17. Reviewer: an independent agent with no part in S4-6.
Reviewed: branch `p4/theorem-c-centre` at `36c0c802`:
- notes: `problems/P4/notes/s4-6-centre.md`, and `s4-5-assembly.md` §9 with the
  passages S4-6 replaced;
- source: `src/p4/validated/tc_layer_centre.py` (sha256 `ec32ccd6…`) and
  `tc_layer_centre_mv.py`;
- scripts: `s46_centre.py` (`735336cc…`), `s46_sample.py`, and the merge rule
  in `s43_parts.py`; tests: `test_tc_layer_centre.py`, `test_s43_merge.py`;
- data: the 60 parts (`parts/centre_t14/`, `parts/centre_t224/`),
  `centre_nb_identity.json`, `centre_sampling.json`, `centre_probe/` and the
  merge output `tc_commutator.json`.

`P4/` means `problems/P4/`.

How the review ran:
- I worked in a git worktree on branch `p4/theorem-c-centre-review`, made from
  `36c0c802`. The gitignored tube caches were copied in; their content hashes
  equal the cache records'.
- No tracked file was edited. The recompute and its merge did overwrite the
  tracked parts and `tc_commutator.json` in that worktree. I compared them
  (§8), then restored them with `git checkout --`.
- Everything ran at nice 10, with at most three of my processes at once:
  1.6 core-h in all (§9).
- My scripts and outputs are in `P4/results/theorem_c/review_s46/`. Its
  README gives each command.

**Overall verdict.**
- **I found no blocking flaw, and no should-fix flaw in the certificate.**
  - Each committed `centre_t14` part certifies the following on its piece.
    For every y > 0 in the piece, every ε ∈ [0, 6] and every τ ≥ 14:
    - dT, κdB, dB and nB are all nonzero;
    - Σ_c = S/(1 + (yτ)²) ≤ b < 0, hence S < 0.
  - That is exactly what assembly §4.1 (S3, S9) and I8 ask of the centre.
  - Together with the 546 `layer_t14` pieces and the exact cover, the strip's
    pointwise condition holds on all of x ∈ (−∞, 0).
- **The identity behind the closing piece is right, and stronger than the
  note says.**
  - I re-derived it exactly with independent code: U6 rebuilt as U1/Dp⁶ from
    the unreduced pair, and symbolic series substituted directly.
  - The y⁸ and y⁹ coefficients of nB's κ¹…κ⁷ coefficients vanish for every
    background with A = 1 + O(y²). The evenness of the background is not
    needed.
  - On the certified data, the dropped coefficients lie inside the TM's
    exact-coefficient range (ve ≥ 25) and contain 0 in all three jet slots.
    The y⁹ ones are exactly 0.
- **The cover is exact.**
  - Replays of c00 and c29 in both τ modes, and of c27 and c28 at τ ≥ 14,
    reproduce the committed parts bit for bit.
  - Their logged box trees prove, in rational arithmetic, that the accepted
    boxes cover piece × [0, 6] × [τ_min, ∞). On c00 this starts from y = 0
    and s = 0, and on every piece it runs through τ = ∞.
  - Independent direct values at 13818 points inside all 6873 accepted boxes
    of six logged covers, and at 231 adversarial points, all lie in the
    boxes' enclosures.
- **The recompute agrees.** All 60 parts, recomputed in the review worktree,
  are identical in every certificate field. The merge on them reproduces
  `tc_commutator.json` except for its git head and merge time.
- **The merge rule is strict and opt-in.**
  - 50 end-to-end attacks on the real merge behave as the rule requires.
  - The default merge and `--substitute-layer-t14` are unchanged: I compared
    main's script with the branch's on the same parts.
- **Nits (§11), none affecting the certificate:**
  - the `Z_LARGE` and `S_SPLIT` comments are stale, and so is the `s_split`
    value every part records;
  - the evenness argument and a few other phrases need rewording;
  - two pre-existing items: an ungated slow test, and a `ve` bookkeeping
    looseness in pinned `tc_tm`.
- **Ledger (§12).** P4-9 and P4-10 merit CERTIFIED.
  - The decisive steps are machine-checked.
  - The written lemmas are the routine ones the S4-5 review accepted, plus two
    routine centre steps that now have exact checks.
  - The README requires the certificate to be committed under
    `certificates/`, and these artifacts are still under `results/`. Promote
    the rows together with the move, not before.

## 1. Verdicts

Basis: **X** exact identity checked by machine; **M** ball/interval
certificate or computation; **W** written argument.

| # | Item | Verdict | Basis |
|---|---|---|---|
| A1a | What the centre must show (§4.1 S3, S7, S9; I8) | re-derived (§2); the parts' statement implies it | W |
| A1b | Σ_c has the sign of S on the whole domain | PROVEN AS STATED: 1 + s² > 0; both regimes' normalizations checked in sympy | X |
| A1c | dT, κdB, dB, nB ≠ 0 on every accepted box, not only S < 0 | PROVEN AS STATED: every regime evaluation inverts all four pivot forms and refuses a hull containing 0 (§2) | code + M |
| A2a | Layer 8 of U6 = −nB vanishes for κ¹…κ⁷; no odd layer | PROVEN AS STATED, independently, in exact arithmetic | X |
| A2b | The y⁸ and y⁹ coefficients of nB's κ¹…κ⁷ coefficients vanish on the background | PROVEN, in a stronger form: for every background with A = 1 + O(y²) | X + W |
| A2c | The certified background is even; A = 1 + y²Â | PROVEN for the computed coefficients (odd ones exactly 0 up to K = 30); Â agrees with the system's own series | M |
| A2d | `drop_lead` drops only coefficients that vanish in truth, within `ve`, jet slots included | PROVEN AS STATED: ve ≥ 25, balls ∋ 0, the y⁹ ones exactly 0; `raise_nb` drops exactly two | M |
| A2e | nY's offsets; the sonic (tt·dT_k)′ raise | nY at offset 10; the raise never fires on centre data | M |
| A3a | GV folding, the Neumann remainder, negative exponents | PROVEN AS STATED (sympy; 20000 randomized containment cases with wide hulls) | X + M |
| A3b | Enclosures against direct evaluation at adversarial points | no failure: y down to 1e−40, τ up to 1e40, ε ∈ {0, 0.03, 6}, grid edges | M |
| A4 | The regime closed forms; the intersection valid on the whole box | PROVEN AS STATED (sympy; each regime's hull contains every box point) | X + code |
| A5a | c00 down to y = 0, from s = 0 | PROVEN (exact tree check, both modes) | X |
| A5b | τ → ∞ with closed boxes | PROVEN (pseudo-boxes; τ ranges tile [τ_min, ∞]) | X |
| A5c | All of ε ∈ [0, 6] | PROVEN (closed cells tiling [0, 6]) | X |
| A5d | Im κ ≤ −14 | by S13 (conjugation), accepted by the S4-5 review; nothing centre-specific | W |
| A5e | The 576-piece x-cover; c29 against the deep end | PROVEN with independent exact code; junction overlap 0.0600 | X |
| A6 | The certified (a*, μ*) balls, never midpoints | PROVEN AS STATED | code + M |
| B7a | `--substitute-centre` is strict and opt-in | PROVEN AS STATED (50 scenarios) | M |
| B7b | The default merge and `--substitute-layer-t14` are unchanged | PROVEN (main's code against the branch's: core fields identical) | M |
| B8 | Recompute | all 60 parts identical; the merge on them identical except head and time | M |
| B9 | Tests | 401 passed, 12 skipped by their gates, 1 deselected (pre-existing, known slow) | — |
| B10 | Pinned modules; the 546 layer substitutions | PROVEN: byte-identical to main; all 1092 layer parts current; 546 substituted | — |
| C11 | The notes | accurate; wording nits in §10 | — |
| C12 | The `Z_LARGE` comment | stale; harmless, since every form is a pure w-form (§10) | — |

## 2. What the centre has to show (A1)

**From assembly §4.1.**
- Fix κ = ε + iτ with ε ∈ [0, 6] and τ ≥ 14.
- S7 needs S = K + ν²/4 < 0 at every x ∈ (−∞, 0), strictly and pointwise.
  Then ∫g|S||u|² ≤ BT(x₂) − BT(x₁) → 0 forces u ≡ 0.
- S3 to S5 and S12 need the following at every x:
  - dT ≠ 0. By the exact factor structure (S4-5 review §5) this gives
    D2, Bc, D2r (so dn, e and g_N), SΔ̃ and Bcr all nonzero;
  - nB, dB ≠ 0, so that B0 is finite and nonzero and (log B0)′, (log B0)″
    exist.
- The centre's share of the line is y = eˣ ∈ (0, e^{−4.56}] at least, where
  the deep pieces stop. It must reach y → 0 (x → −∞), and it must hold for
  the closed ranges ε ∈ [0, 6] and τ ∈ [14, ∞): P4-9's claim includes
  |Im κ| = 14.
- I8 says the same.

**The parts' statement implies it.** A passing part certifies, for every
y ∈ [lo, hi] with y > 0, every ε ∈ [0, 6] and every τ ≥ 14:
- dT, κdB, dB, nB ≠ 0;
- Σ_c ≤ b < 0.

Σ_c = S/(1 + s²) with s = yτ real, so 1 + s² ≥ 1 > 0 and Σ_c has the sign of
S. Both normalizations reduce to it in sympy (`sympy_forms.json`):
- regime L: Σ = S/s², and Σ_c = Σ·s²/(1 + s²);
- regime S: Σ = S, and Σ_c = Σ/(1 + s²).

**Nonvanishing is enforced, not inferred.**
- `sigma_graded` (`tc_layer_centre.py:613-617`) forms, in every regime
  evaluation and unconditionally:
  - V and V′, which divide by dT;
  - Y, which divides by κdB;
  - (log dB)′ and (log dB)″, which divide by dB;
  - (log nB)′ and (log nB)″, which divide by nB.
- Each division is `GV.inv` of the pivot form's bracket Q, with D = yᵐzᵖQ.
  `GV.inv` raises if the leading coefficient's ball q₀ contains 0, or if the
  hull of 1 + u contains 0. So Q = q₀(1 + u) ≠ 0 at every point of the hull.
- With y > 0 and z = yκ ≠ 0, each of the four denominators is nonzero at
  every point of the box.
- `sigma_box_centre` returns an enclosure only if some regime succeeded. A
  failure raises, `_upper` turns it into +∞, and the box is split, never
  accepted.
- Empty boxes have no points, and §5 checks that exactly.

## 3. The identity and the raise (A2)

**Independent exact check** (`identity_independent.py` →
`identity_independent.json`, 7 s). The script never calls `tc_layer_centre`.

1. **U6 built a second way.** U1 = Bc·UN3c − Ac·dBc·D2 is built from the
   unreduced reduction pair, the numerator the S4-5 review checked in sympy
   (I3). Exact division gives U1 = Dp⁶·U6 with zero remainder. This U6 equals
   the formula `build_pairs` evaluates (Bcr·UN3r − Acr·EBc·D2r) and
   `check_nb_identity` uses. It has 34291 terms.
2. **Symbolic series substitution, with no layer bookkeeping.** Every monomial
   of U6 gets the truncated series:
   - n = n₀ + n₁y + n₂y², and likewise w and v;
   - Â = a₀ + a₁y + a₂y², with A = 1 + y²Â, N = n/y, W = y²w, V = yv.

   Odd coefficients are included. The script collects the exact coefficient of
   y^k κ^j for k ≤ 10, as a polynomial in all twelve background symbols.
3. **Results.**
   - For j = 1…7, the lowest monomial layer is 8, and the y⁸ and y⁹
     coefficients are the zero polynomial. The y¹⁰ coefficient is nonzero and
     involves only n₀, w₀, v₀ and a₀.
   - For j = 0, 8 and 9, nothing appears below y¹⁰.
   - U6 has no odd-parity monomial.
   - For dB, the y⁸ coefficient is nonzero for κ¹…κ⁷. So dB rightly keeps
     offset 8.
4. **Every atom read by `build_pairs`** (value, d- and dd- slots) involves only
   A, N, W, V and κ. `tc_energy.compile_mp` keeps only those exponents, so
   nothing is silently dropped.
5. **Controls.** Each of these is reported:
   - A = 1 + αy + y²Â makes y⁹ nonzero for all of κ¹…κ⁷;
   - adding N²W⁵κ³ makes y⁸ nonzero at κ³;
   - adding N·W⁵κ² gives one odd-parity term, and y⁹ becomes nonzero at κ².

**Why y⁹ needs no evenness.**
- The y⁹ coefficient gets contributions from two places:
  - layer-9 monomials at order 0. There are none, since there is no odd
    parity;
  - layer-8 monomials at order 1. Since A = 1 + O(y²), this is the derivative
    of the layer-8 polynomial L₈(n, w, v) in the direction (n₁, w₁, v₁).
- L₈ is identically 0, so that derivative vanishes.
- The note's route instead uses the evenness of the background. That is also
  valid, but only the order-0 and order-1 coefficients matter, and those are
  exactly 0 in the certified series.
- The A-structure, unlike evenness, is load-bearing (control 5a). It holds by
  construction: `centre_components` builds A = 1 + y²Â from the momentum
  constraint.

**On the certified balls** (`raise_inputs_check.json`; pieces c00, c01, c15,
c28 and c29):
- **The (a, μ) balls** are rebuilt here from `a3_midpoint.json` with my own
  parsing.
  - They contain the box midpoint + `K_minus_m`.
  - They exclude the bare a3 midpoint, the bug fixed before 2026-10-16.
- **The certified series.**
  - All 45 odd coefficients of n, w, v up to K = 30 are exactly zero balls.
  - The tail radius ν = 0.08 ≥ YCUT.
- **The components.**
  - Offsets are (0, −1, 2, 1).
  - A's constant term is exactly 1 and its y¹ coefficient exactly 0.
  - Â's 29 TM coefficients overlap the system's own Â series
    (`CentreExpansion.Ahat_series`, from `systems.centre_Ahat_fraction`).
- **nB's κ¹…κ⁷ coefficients, in all three jet slots and on every piece checked.**
  - The trimmed offset is 8, with `ve = 25`. So `drop_lead(2)` drops exact
    Taylor coefficients, and its `k ≤ ve` guard is far from binding.
  - The y⁸ coefficient is a ball around 0 with radius 295 to 7.4e6.
  - The y⁹ coefficient is exactly 0, radius 0. Only the y⁸ drop relies on the
    identity.
  - The first kept coefficient (y¹⁰) is 3e15 to 4e18. Across slots 0, 1, 2 it
    scales 1 : 10 : 100, as y d/dy requires.
  - `raise_nb`'s output equals dropping exactly those two coefficients.
- **nB's κ⁰, κ⁸ and κ⁹ coefficients** are at offset 10. The pivot κ⁹ at y¹⁰ is
  9.899e13, the note's "9.9e13".
- **The sonic derivative raise** `LayerData` applies to (tt·dT_k)′ never fires
  on the centre data: der.off = val_off ≥ target. It would be sound if it did,
  because d/dx = y d/dy keeps Laurent offsets.

**Jet slots.**
- Slots 1 and 2 are the true x-derivatives, by `AtomKP.jet_kp`:
  X′ = D1(X)/SΔ̃ on the true background.
- The y⁸ and y⁹ coefficients of y d/dy f are 8f₈ and 9f₉, both 0.

**`drop_lead` semantics** (`tc_tm.py:34-111`).
- A TM is f = P + e: P's coefficients below `ve` enclose the exact Taylor
  coefficients, and the remainder e is pointwise with valuation ≥ ve.
- Under that reading, dividing by ℓᵏ for k ≤ ve, with the k lowest true
  coefficients 0, gives a valid TM.
- `ve` is propagated correctly by add, mul, shift, trim and the spill into
  CAP − 1.

## 4. The graded arithmetic and the regimes (A3, A4)

**The identity used to assemble S.**
- `sigma_graded` transcribes `tc_layer.sigma_box`'s formula T(λ).
- I re-derived T(λ) = λ²S/τ² in sympy for a free real λ. S comes from the
  S4–S6 definitions, and the blocks are what `LayerData`'s ratios mean
  (V = λ(A0 + κtt), Y = λ²(B0 − u2κ²)/κ, …).
- The residual is 0. Flipping one sign gives a nonzero residual.
- So λ = g⁻¹ is legitimate in both regimes. The S4-5 review's I6, which the
  note cites, also used a free λ.

**The closed forms** (sympy, all true):
- **Regime L:** 1/τ = gσ, z = is + gε, ζ = σ(−i + yεσ)/(1 + (yεσ)²) and
  w = yζ.
- **Regime S:** y = rs, z = is + rsε and w = r(rε − i)/(1 + r²ε²).
- **The fold** g^a = g^{K+1}g^{a−K−1} is correct.
- **The Neumann identity** 1/(1 + u) = Σ_{k≤K+1}(−u)^k + (−u)^{K+2}/(1 + u)
  is correct.

**The `GV` class** (`gv_stress.json`: 20000 random cases at g-hulls up to 0.5,
some touching 0).
- `value()`, products, powers, inverses and inverses of products all contain
  the exact value at a random point of the hull and of every coefficient ball.
- Exactly-zero negative-exponent coefficients are accepted. A 1e−40 ball at
  g⁻¹ is refused.
- 73 inverses were refused, where 1 + u reached 0.
- Reading the code adds three points.
  - **Exact zeros are rigorous.** A coefficient ball that is exactly zero is
    a rigorous statement, because ball arithmetic is inclusion-isotone.
  - **g-dependent coefficients are fine.** Some coefficients depend on g: the
    TM values over the strip, D, and Dr. Every operation keeps the identity
    value(point) = Σ g^a c_a(point) with c_a(point) in the ball, and folding
    preserves it.
  - **The remainder uses the hull of u.** The remainder (−u)^{K+2}·(1/(1 + u))
    multiplies by 1/(1 + hull(u)), which contains 1/(1 + u(point)) at every
    point.

**The intersection of regimes is sound.**
- Each regime encloses the same function, Σ_c, over a set that contains every
  point of the box:
  - Regime L reads (y ∈ L, σ ∈ [1/s_hi, 1/s_lo]) and ignores τ, whose value
    is s/y.
  - Regime S reads (r ∈ [1/τ_hi, 1/τ_lo], s ∈ [s_lo, s_hi]), with the TM
    values over the y-strip, which contains the actual y = rs.
- If the two enclosures were disjoint the evaluator would raise, and the box
  would be split, not accepted.
- Only `regimes_of` decides applicability: L needs s > 0, S needs bounded s.
  `S_SPLIT` and `Z_LARGE` play no part in production (§10).

**Adversarial points** (`adversarial_points.json`).
- The points:
  - on c00, y ∈ {1e−40, 1e−30, 1e−20, 1e−13, hi(1 − 1e−12), hi},
    τ ∈ {14, 14(1 + 1e−12), 224, 2e4, 1e8, 1e15, 1e40} and
    ε ∈ {0, 0.03, 6};
  - on c29, y at both ends and in between, ε ∈ {0, 3.4, 6}, and τ putting s
    at s_lo, 1e−3, 0.5, 8 and 1e3, or τ = 14, 224, 1e12, 1e30.
- For each point, every accepted box of the logged covers that contains it is
  found exactly.
- The independent direct value lies in every such box's enclosure.
  - 231 points give 396 point-in-box memberships, since a boundary point lies
    in up to 4 boxes.
  - 0 failures, 0 uncovered points, 0 undecided values.
  - The largest value is −6, at the s → 0 corner.
- The direct value is my own chain (`direct_chain.py`), written from the
  assembly's definitions. It runs at the midpoint background of the certified
  balls, with precision raised until tight. On the same pairs it gives the
  same balls as `tl.point_S_direct`, to all printed digits, radii included.

**The optional mean-value form** (`tc_layer_centre_mv.py`).
- No part uses it, and the merge rejects any part that records it.
- I only read it. Its centres are exact ball midpoints: there is no float
  ι-centre, the S4-4 3d family.
- δ is frozen with ρ bounded over the true y-range.

## 5. The cover (A5)

**Replays** (`cover_replay.py`).
- `certify_logged` is `tc_layer.certify_piece` copied statement for
  statement, with only logging added.
- It ran on the parts' own inputs (`centre_layer_data(..., raised=True)`,
  `sigma_box_centre`, target 0, 150000 evaluations, precision 320).
- Worst bound, boxes, evaluations, empty boxes, the depth histogram and the
  cover parameters all equal the committed parts':
  - c00 at τ ≥ 14 (505 boxes / 680 evaluations) and at τ ≥ 224 (433 / 608);
  - c29 at τ ≥ 224 (137 / 172) and at τ ≥ 14 (3898 / 29924, 640 s);
  - c27 at τ ≥ 14 (533 / 2920, 61 s) and c28 at τ ≥ 14 (1367 / 9592, 200 s),
    the other two thin-margin pieces (bounds −1.6e−3).

**Exact verification of the logged trees** (`cover_verify.py`, no p4 code,
Fractions, ∞ kept as ∞). A box means P(b) = {y ∈ [la, lb], y > 0,
ε ∈ [e0, e1], τ ∈ [t0, t1] finite, yτ ∈ [s0, s1]}.
1. **The roots.**
   - Each has y-range [lo, hi].
   - Their ε-cells tile [0, 6].
   - Per cell, the s-cells tile [s_first, ∞) with s_first ≤ lo·τ_min.
   - Per s-cell, the τ ranges tile [τ_min, ∞).
2. **Every evaluated node.**
   - The exact projections of P(raw) onto y, τ and s lie inside the
     normalized box. So P(raw) ⊆ P(nb) ⊆ the hulls the enclosure was computed
     on.
   - A node called empty has empty P(raw).
3. **Every split** cuts nb at one point of one coordinate, so the children's
   union is P(nb).
4. **Every leaf** is accepted with upper bound < 0, or empty. No node failed,
   and the tree is consistent.

All six trees pass with 0 errors. Five corruptions of the c29 tree are
each reported:
- nb shrunk below raw;
- a sliver cut out of a split;
- a nonnegative accepted box;
- a false empty;
- a missing root cell.

By induction the accepted boxes cover D = {y ∈ [lo, hi], y > 0} × [0, 6] ×
[τ_min, ∞). That includes y → 0 on c00 (s from 0, with the pre-split τ slabs)
and τ = ∞ through the pseudo-boxes on every piece. The enclosures over each
accepted box are tested next.

**Points inside the accepted boxes** (`point_check.json`).
- Every accepted box of the logged trees got points: y from a per-tree grid
  inside the box's y-projection, τ log-uniform in its τ-range at that y (up to
  1e18), and ε at the cell ends or inside.
- The result covers 13818 points in the 6873 accepted boxes of the six trees
  (c00 and c29 in both modes, c27 and c28 at τ ≥ 14): **0 failures and 0
  undecided**. Every value lies in its box's logged enclosure [lower, upper].
- The values nearest a box's upper bound show tightness, not near misses:
  - in the c00 corner boxes (s ≤ 1e−3) the enclosure is
    [−6.00000014, −5.99999986], around the value −6.0;
  - in the large-s boxes the enclosure's upper end, −595.4, is nearly
    attained at the cell edge ε = 4.6, where k₂(ε)s²/(1 + s²) is largest
    (value −596.7).
- The largest sampled Σ_c is −6 on c00 and −6.131 on c29, as the note's
  samples say.
- The thinnest margins sit in boxes whose enclosures are wide. The c29 box
  with upper bound −5.9e−4 has enclosure [−13.4, −5.9e−4] at
  y ∈ [0.01002, 0.01017], ε ∈ [0.12, 0.2], s ∈ [0.140, 0.164] and
  τ ∈ [14, 16.4].
- The cover stops refining as soon as the upper bound is negative. With
  rigorous arithmetic that is sound, and the in-box points confirm the
  enclosures.

**The x-cover and the junction** (`cover_x_full.py`, no s43_parts code).
- The 546 non-centre pieces sit at their exact cache steps; their union is
  exactly [−5134103575202365/2⁵⁰, 0] = [−4.56, 0], with no gap.
- The centre pieces touch exactly (hiᵢ = loᵢ₊₁ as Fractions), with
  c00.lo = 0, so their union is [−∞, log YCUT].
- log hi(c29) minus the deep end = 0.059999999999999571913 ± 4.3e−22 > 0.
- The merge's own exact cover agrees: removing c29 opens exactly
  (−4.970, −4.56) (§7).

**Im κ ≤ −14.**
- This is S13 (conjugation), a written lemma the S4-5 review accepted (§6).
- It acts on the eigenvalue problem as a whole, so nothing about the centre
  pieces needs separate treatment.

## 6. The certified balls (A6)

- `s46_centre.Ctx` builds the series with
  `ec_centre.certified_centre_expansion(256, K=30)`.
- `certified_centre_balls` returns the a3 midpoint plus `K_minus_m`, read
  from `a3_midpoint.json`. My own parsing confirms that these balls contain
  the box and exclude the bare midpoint (§3).
- `centre_layer_data` passes that `ce` to `centre_components` and then
  `build_pairs`. No code path in the parts touches a midpoint.
- Midpoints appear in two places, both clearly marked as checks and not
  certificates:
  - `s46_sample.py`'s `direct_mid`, which uses the midpoint of the certified
    balls, a point inside them;
  - my `direct_chain.py`.
- `centre_checks` tests the record's `inputs.centre_series` string. The real
  guarantee is the pinned driver, which hard-codes the certified expansion
  (§11, nit 3).

## 7. The merge rule (B7)

**Code read** (`s43_parts.py:1995-2237` and the `cmd_merge` changes). The rule
is opt-in through `--substitute-centre` and always strict.
- **Per part:**
  - current hashes, with accepted variants not counting;
  - no late module, a passing write check, a git head and a clean tree;
  - the required modules recorded, and the mean-value modules not;
  - the driver pinned to the current `s46_centre.py`;
  - status pass, ok, no failing box, and one finite negative bound;
  - τ_min of the record, the cover and the scope equal and ≤ b₁, with τ up
    to ∞;
  - ε-cells tiling [0, 6];
  - the cover's y-range equal to the plan's [lo, hi];
  - s_lo ≤ lo·τ_min, and s_lo = 0 on c00;
  - laurent on and combine off;
  - the evaluator, the identity block, and the TM domain [0, hi].
- **Global:**
  - checks.json valid for the current code;
  - b₁, c and EPS_MAX unchanged;
  - the identity artifact valid.
- A failed or errored S4-3 part goes to manual review.

**End to end** (`merge_tamper.py`).
- The script runs the real `cmd_merge` on the committed parts. It serves
  edited copies of chosen parts, of the identity artifact and of checks.json,
  and writes its output to a temporary directory.
- All 50 scenarios behave as the rule requires.
  - **The baseline** is complete and certified, with 576 substituted.
  - **Each single-field tamper** leaves the merge partial and uncertified,
    with the gap at that piece: bound, status, ok, fails, a stale hash, a
    missing required module, the MV module recorded, a dirty tree, no head, a
    late module, a failed write check, the driver, an id swap, the spec, the
    TM domain and the series string. So do these cases:
    - a τ ≥ 224 part served as the τ ≥ 14 part;
    - its scope edited to τ ≥ 14;
    - its τ_min forged in all three fields. `s_lo = lo·224` still gives it
      away;
    - a τ end of 1e6;
    - an ε-cell removed, or cells ending at 5;
    - the y-range halved;
    - s_lo ≠ 0 on c00, or s_lo too large;
    - the evaluator, laurent, combine, raised or the identity block changed.
  - **A missing c00 or c15 part** leaves the matching gap, (−∞, −18.13) or the
    c15 interval.
  - **Each bad identity artifact blocks all 30 centre pieces**: missing,
    ok false, a layer-8 entry, odd terms, a wrong κ list, a stale hash, a
    missing module, or a wrong driver. So does a checks.json that is invalid
    or stale.
- **The inherent limit.** A part forged consistently in every field, bounds
  or τ_min and s_lo together, is accepted. That is a merge that doesn't
  recompute, the same as for `layer_t14`. Integrity rests on git history and
  on the aggregate hash the ledger cites (`4ac8a7b0…`), which §8 re-derives.

**Unchanged paths** (`merge_unchanged.py`).
- I ran main's `s43_parts.py`, from `git show main:` and executed in-process
  with its real path, against the branch's, on the same parts and in both
  modes.
- In both modes the core fields are identical: ok, status, cover, detail,
  missing, failed, substitution, `strip_exclusion_certified`, checks,
  conflicts and the layer summaries.
  - The default merge is partial, with 576 missing.
  - `--substitute-layer-t14` is partial, with 546 substituted, 30 missing and
    the gap (−∞, −4.56), as in `merge_substitute_layer_t14.json`.
- The only differences are:
  - the S4-5 review's wording fix in `scope`: `full_strip_tau` and
    `xi_open_max` replaced, and `tau_target`, `tau_by_carrier` and
    `open_at_centre_end` added;
  - the new `centre_substitution` block;
  - the `phase_settings` note;
  - provenance.

## 8. Recompute, tests, provenance (B8–B10)

**Recompute** (`recompute_compare.json`).
- All 60 parts were recomputed with `s46_centre.py work --force`, in three
  disjoint detached groups, in the review worktree.
- Each recomputed part was compared with `git show HEAD:` of the same path.
  All 60 are identical in every certificate field: status, bound, worst,
  boxes, evaluations, empty, depth histogram, cover parameters, ladder
  attempts (target, ok, evaluations), forms, engine, identity, inputs, scope,
  recorded source hashes and driver.
- Only timing, env, dates and the git head (`36c0c802`) differ. The CPU was
  0.71 core-h, as the note says.
- Every new part records a clean tree, no late module, a passing write check
  and the same source hashes.
- The merge CLI on the recomputed parts (`merge_recomputed.json`) gives
  status complete and `strip_exclusion_certified = true`, with 576
  substituted (546 + 30) and nothing missing, failed, rejected or blocked.
  That is identical to the committed `tc_commutator.json` (`6c228a73…`)
  except for `provenance.git_head` and `provenance.merged`.
- The tracked files were then restored with `git checkout --`.

**Tests.**
- Run: `PYTHONPATH=problems/P4/src python -m pytest problems/P4/tests` on all
  414 tests (`logs/pytest_p4_full.txt`).
- Result: 401 passed, 12 skipped and 1 deselected, in 901 s.
  - The 12 skips are the suite's own gates (`P4_SLOW=1`, `P4_TUBE_CACHE`),
    reported by `-rs`.
  - The deselected test is
    `problems/P4/tests/test_tc_commutator.py::test_centre_first_piece`.
- About the deselected test:
  - It runs the old S4-3 centre path on one piece.
  - It is not gated, and it is documented not to finish in 15 minutes
    (`s4-3-recovery.md:83`, and the test file's docstring).
  - My first attempt, which I stopped after it had spent about five minutes
    in that test, is kept as `logs/pytest_p4_full_attempt1_stopped.txt`.
  - The deselected test predates S4-6 and is not S4-6 code.
- `test_tc_layer_centre.py` (22 tests) and `test_s43_merge.py` passed inside
  the full run.

**Provenance** (`provenance_check.json`).
- Every recorded source hash is current. That covers the 60 centre parts (all
  at git head `e9cc5086`, with a clean tree, no late module, a passing write
  check and driver `735336cc`), the identity artifact (at `da7a7f28`, clean)
  and all 1092 `layer_t14`/`layer` parts.
- The 26 sources the layer parts pin are byte-identical on main and at HEAD,
  and so is `s44_layer.py` (`3a86e50f`).
- `git diff --stat main HEAD -- problems/P4/src` lists only the two new files.
  The same holds from `fda58897`.
- Between `e9cc5086` and `da7a7f28` and HEAD, neither `src` nor
  `s46_centre.py` changed.
- Every sha256 and aggregate the notes quote matches: `centre_t14`
  `4ac8a7b0…`, `centre_t224` `ba424c3f…`, `layer_t14` `e9835a41…`, the
  identity `0fad868b…`, `tc_commutator.json` `6c228a73…` and the sampling
  `a740dd19…`.

## 9. Computations

All ran in the review worktree at nice 10.

| script (`P4/results/theorem_c/review_s46/`) | result | CPU |
|---|---|---|
| `identity_independent.py` | U6 two ways equal; y⁸, y⁹ ≡ 0 for κ¹…κ⁷ on every A = 1 + O(y²) background; controls reported | 7 s |
| `raise_inputs_check.py` (5 pieces) | balls, series parity, A-structure, ve ≥ 25, dropped balls ∋ 0, raise = drop of two | 174 s |
| `sympy_forms.py` | T(λ) = λ²S/τ² for free λ; regime forms; Neumann, fold | < 1 s |
| `gv_stress.py` (2 × 20000) | no containment failure | 7 s |
| `cover_replay.py` | c00 ×2, c29 ×2, c27, c28: bit-identical to the parts | about 1100 s |
| `cover_verify.py` | six trees exact, 0 errors; 5 controls reported | < 2 s |
| `point_check.py` | 13818 points in 6873 accepted boxes (six trees): 0 failures, 0 undecided | 368 s (+ 144 s first passes) |
| `adversarial_points.py` | 231 points / 396 memberships: 0 failures, 0 uncovered | about 20 s |
| `cover_x_full.py` | 576 pieces, no gap, junction 0.0600 | 1 s |
| `merge_tamper.py` | 50/50 scenarios as required | 20 s |
| `merge_unchanged.py` | default and layer_t14 merges unchanged in their core | 2 s |
| `provenance_check.py` | all current; pins identical to main; quoted hashes match | 3 s |
| `notes_check.py` | the §5 table, CPU sum, merge and sampling claims match | 1 s |
| `s46_centre.py work --force` (60 parts) + `recompute_compare.py` | identical; merge identical | 2562 s (0.71 core-h) |
| `s46_centre.py identity --out` → `identity_replay.json` | identical in content (identity block, recorded source hashes, driver); only `created`, `env`, `setup_s` and the code record's head and time differ | 6 s |
| pytest, P4 | 401 passed / 12 skipped / 1 deselected | 901 s (+ 480 s, stopped attempt) |

Total: about 1.6 core-h, within the 3 core-h budget.

## 10. Documentation (C11, C12)

**`s4-6-centre.md`** is accurate, checked number by number (`notes_check.json`):
- the §5.2 table (30 rows of bound, boxes, evaluations and CPU);
- §5.1's CPU, 0.707 core-h;
- §5.3's merge claims, including the worst bound per region;
- §5.4's sampling claims: 14400 points, 0 failures, 13795 finite direct,
  1.14e−13, 3.07e−6, and −6.112 / −6.051 / −6.017 on c29 / c28 / c27.
- The §4 table's test names exist, and the file holds 22 tests.

Wording nits:
- §3.2, "Why the y⁸ and y⁹ coefficients vanish."
  - The argument asserts that the background is even. The evidence covers the
    coefficients up to K, and the argument needs only orders 0 and 1.
  - Better: "the y⁹ coefficient is the derivative of L₈ ≡ 0 in the direction
    (n₁, w₁, v₁) (A = 1 + O(y²)), so it vanishes on any background; on the
    certified one n₁ = w₁ = v₁ = 0 anyway."
  - Also say that the y⁹ TM coefficients are exact zeros, so only the y⁸ drop
    rests on the identity.
- §3.2, "Statement" says the identity uses "the compiled atoms of
  `tc_energy.Atoms`". It uses the exact fmpq_mpoly atoms (`at.mp`). The
  compiled term lists (`at.cmp`) are what evaluation reads. They are the same
  polynomials, and `compile_mp` drops no variable these atoms use (§3).
- §3.2, "against ~1e14 for the true leading terms": that is nB's κ⁹ pivot,
  9.9e13. The κ¹…κ⁷ y¹⁰ coefficients are 3e15 to 4e16.
- §1 writes y ∈ (0, e^{−4.5}]. The pieces cover (0, YCUT], with
  YCUT = e^{−4.5}(1 − 3.7e−17), as the assembly's I8 says. The difference is
  harmless, since the deep pieces overlap by 0.06.
- §5.4 calls the samples independent evidence. Their "graded" column
  evaluates point boxes only, and the box enclosures were tested only on the
  4 boxes of the tests. §4 and §5 here close that gap.

**`s4-5-assembly.md`.**
- Every S4-6 hunk replaces lines one for one (equal old and new counts), and
  §9 is appended after the footer. So no cross-referenced line number moved.
- The replaced passages are correct: I8's status, S10, S11, A2, the §4.4 row,
  §5 items 1 and 2, the P4-9 "level now", and §7.
- §9 matches the S4-6 note.

**The `Z_LARGE` comment (C12).**
- The comment: `Z_LARGE = 8.0  # regime L: small-z forms for s <= Z_LARGE,
  large-z (zeta) forms above`.
- What the code does: regime L always uses the large-z (ζ) forms
  (`use_large = rg == "L"`). `Z_LARGE` only decides whether `Grade` builds a
  `Z`, which nothing then uses.
- Does it mislead? Only about which form regime L evaluates. A reader
  checking regime L's small-z algebra for s ≤ 8 checks code that never runs.
- No soundness consequence: on the raised data all 19 forms are pure w-forms
  (pivot = top), so the two forms coincide, and the large form is valid for
  every s > 0 anyway.
- **`S_SPLIT` has the same problem, and the note does not mention it.**
  - The comment says "regime S below, regime L above". In production both
    regimes are evaluated wherever they apply, and `S_SPLIT` only steers
    `regime_of`, the default when `sigma_graded` is called without a regime.
    Only the tests do that.
  - Every part records `engine.s_split = 0.5`, which suggests a split that was
    not used.
- Neither can be fixed in place: the 60 parts record the module's hash. Fix
  them in the note, or at the next deliberate change of the module.

## 11. Findings, ranked

**Blocking: none found.** Checked:
- the identity, independently and exactly;
- the raise on the certified data;
- the graded algebra, in sympy and by randomized containment;
- the cover, by exact tree verification of the replayed parts;
- the enclosures, at points inside every accepted box and at adversarial
  points;
- the x-cover with independent code;
- the merge, under 50 attacks;
- the recompute of all 60 parts;
- provenance.

**Should-fix: none.**

**Nits:**
1. **The `Z_LARGE` and `S_SPLIT` comments are stale, and so is the recorded
   `engine.s_split`** (§10).
   - *Failure scenario:* a reader concludes that regime L uses the small-z
     forms below s = 8, or that a box with s < 0.5 was evaluated in regime S
     only, and audits the wrong code path.
   - *Fix:* a sentence in the note now; edit the comments at the next
     deliberate module change.
2. **The evenness wording in §3.2** (§10).
   - *Failure scenario:* a reader takes "the background is even" as a claim
     about the series tail beyond K, which is not shown and not needed.
   - *Fix:* the derivative argument, or "orders 0 and 1".
3. **`engine.raised`, `inputs.*` and `engine.s_split` are constants the
   driver writes, not measurements.**
   - The merge checks them, but the guarantee that the data were raised and
     built from the certified balls is the pinned `s46_centre.py`, which
     hard-codes both.
   - *Failure scenario:* none, while the driver pin holds. With a
     differently-written driver, the record would not prove those facts.
   - *Fix:* say so in §6 of the note.
4. **`test_centre_first_piece` is not gated** (pre-existing, not S4-6).
   - *Failure scenario:* `pytest problems/P4/tests` never finishes.
   - *Fix:* gate it behind `P4_SLOW`, like the other long runs. The test file
     isn't pinned.
5. **`tc_tm.TM.inv` records `ve = m`, ignoring `self.ve`** (`tc_tm.py:157`,
   `:182`; pinned, pre-existing).
   - *The issue:* the recorded exact-coefficient range of an inverse can
     exceed the input's true one.
   - *Failure scenario:* a future `drop_lead(k)` applied directly to an
     inverse, with k between the input's true range and the recorded one,
     would pass the guard wrongly.
   - *Not triggered:* every TM on the centre path is also capped by the
     component TMs (ve 31), and `drop_lead(2)` runs at ve = 25.
   - *Fix:* none possible without invalidating 1152 parts. Record it.
6. **Small wording points in the notes** (§10): "compiled atoms";
   "~1e14"; (0, e^{−4.5}] instead of (0, YCUT]; the description of the
   samples.

## 12. Ledger

**Level.**
- Under README "Ledger levels", a row is CERTIFIED when three things hold:
  - its certificate is committed under `certificates/`;
  - the decisive steps are machine-checked;
  - any written lemmas are routine.
- For P4-9 the decisive steps are all machine-checked:
  - the pointwise condition, S < 0 with dT, dB, nB ≠ 0, on every piece ×
    ε ∈ [0, 6] × τ ∈ [14, ∞): 546 layer parts and 30 centre parts, all
    interval certificates;
  - the exact cover of x ∈ (−∞, 0);
  - the exact identities: E1/E2/E3/E5, the raise identity, S6 and the
    layer normalization;
  - the boundary-term data: the residue identities and the centre layers.
- The written lemmas are S5 to S8, S12 and S13, which the S4-5 review accepted
  as routine given its checks. The centre adds two steps, both routine and
  both now backed by exact checks:
  - the grading argument for the raise;
  - "slots 1 and 2 are y d/dy of slot 0".
- So the proof is not WRITTEN_PROOF. No analytic step is prose-only.
- P4-10 adds P4-6 and P4-8 (both CERTIFIED) and Lemma G1 (reviewed), so it
  merits CERTIFIED with P4-9.
- **The one unmet condition is location.** The artifacts are under
  `results/theorem_c/`, and `certificates/README.md` ("Not here yet") keeps
  them there until promotion.
- **Recommendation:** promote both rows to CERTIFIED in the same commit that
  moves the artifacts below into `certificates/`. Until then, keep them
  CONJECTURED, not WRITTEN_PROOF.

**Artifacts to move into `problems/P4/certificates/theorem_c/` on promotion**
(`git mv`, unchanged, as for the 2026-10-17 move; add them to `SHA256SUMS`
and to the README's old-path table):

| artifact (now under `P4/results/theorem_c/`) | proposed path | regenerated by (repository root) |
|---|---|---|
| `parts/plan.json` | `strip/plan.json` | `python problems/P4/scripts/s43_parts.py plan` |
| `parts/caches/{full,gap,deep,band}.json` (the cache records; content hashes of the gitignored caches) | `strip/caches/` | `python problems/P4/scripts/s43_parts.py build-caches` |
| `parts/layer_t14/` (546 parts, aggregate `e9835a41…`) | `strip/layer_t14/` | `python problems/P4/scripts/s44_recompute.py launch --workers 4 --nice 10` (the six queues of `s4-4-layer-engine.md` §14; about 3 h 52 min on 4 workers) |
| `parts/centre_t14/` (30 parts, aggregate `4ac8a7b0…`) | `strip/centre_t14/` | `python problems/P4/scripts/s46_centre.py launch --workers 4 --nice 10 --tau-mins 224,14 --max-evals 150000` (`--force` in a scratch checkout; about 0.7 core-h) |
| `centre_nb_identity.json` (`0fad868b…`) | `strip/centre_nb_identity.json` | `python problems/P4/scripts/s46_centre.py identity` |
| `tc_commutator.json` (the complete merge, `6c228a73…`) | `strip/tc_commutator.json` | `python problems/P4/scripts/s43_parts.py merge --substitute-layer-t14 --substitute-centre` |
| `review_s46/` (this review's checks) | `review_s46/` | its README |

- Already in `certificates/theorem_c/`: `checks.json` (P4-8, P4-9) and
  `review_s45/`.
- Supporting, optional, not cited for the theorem:
  - `parts/centre_t224/` and `parts/layer/` (the second τ mode, same
    commands);
  - `centre_sampling.json` (`python problems/P4/scripts/s46_sample.py`);
  - `layer_recompute/merge_substitute_layer_t14.json` (historical).

**How to move without breaking the tooling.**
- The drivers find parts through `s43_parts.PARTS`:
  - `s44_layer.py` and `s46_centre.py` are hash-pinned, by the 1092 and 60
    parts respectively;
  - `s43_parts.py` is pinned by nothing.
- So point `s43_parts.py`'s paths for `layer_t14`, `centre_t14`, the plan, the
  caches, the identity artifact and `tc_commutator.json` at the new location.
- Re-run the merge, and check that only `provenance` changed.
- The S4-3 part directories and the second-mode parts can stay where they
  are.

**Paste-ready rows** (CLAIMS.md columns: id | problem | statement | level |
evidence | updated; evidence paths assume the move above):

| P4-9 | P4 | Strip exclusion (Theorem C, part 2): in the class C of P4-8 there is no eigenvalue with 0 ≤ Re κ ≤ 6 and \|Im κ\| ≥ 14. Method: a once-commuted (m = 1) energy identity for χ₁ = χ′ with weight g = exp∫(Re A1 − c·tt), c = −1/32; the pointwise condition S = K + ν²/4 < 0, together with dT, dB, nB ≠ 0, certified on piece × ε ∈ [0, 6] × τ ∈ [14, ∞) for every piece of an exact cover of x ∈ (−∞, 0) by 576 pieces: 546 S4-4 layer-engine pieces on x ∈ [−4.56, 0] (Σ = Sξ²/τ²) and 30 S4-6 centre pieces on y = eˣ ∈ (0, YCUT] (Σ_c = S/(1 + (yτ)²), graded pivot forms in y on the certified (a*, μ*) balls; the closing piece reaches y = 0 through the exact centre layer identity of U6 = −nB, layer 8 ≡ 0 for κ¹…κ⁷ and no odd layer, which raises nB's κ¹…κ⁷ Laurent offsets from 8 to 10); boundary terms γ = 1 + σ₀ + (1 − σ₀)(ε + c) ≥ 1.0617 at the sonic point and BT = O(t^{3+k}) → 0 at the centre (B0 = O(t^k), k ≥ 2); Im κ ≤ −14 by conjugation | CERTIFIED — 2026-10-17: merge --substitute-layer-t14 --substitute-centre reports status complete and strip_exclusion_certified = true (576 substituted, exact cover with no gap; worst bounds: centre −5.9e−4 (c29), deep −0.0217, tube −0.209, sonic −0.452, gap −0.603, band −3.95). The layer-engine review and its re-check cover the 546 pieces; the S4-5 review accepts the written lemmas S5–S8, S12, S13 as routine given its committed checks; the S4-6 centre review (docs/review-2026-10-17-p4-s4-6-centre-adversarial.md) found no flaw: the raise identity re-derived exactly by independent code, all 60 centre parts recomputed identically, the covers of c00, c27, c28 and c29 verified exactly, the merge rule attacked end to end | problems/P4/notes/s4-commutator.md §2–§3; notes/s4-4-layer-engine.md §1–§5, §12–§14; notes/s4-5-assembly.md §4.1, §9; notes/s4-6-centre.md; certificates/theorem_c/strip/{plan.json,caches/,layer_t14/ (aggregate e9835a41…),centre_t14/ (aggregate 4ac8a7b0…),centre_nb_identity.json,tc_commutator.json}; certificates/theorem_c/checks.json; certificates/theorem_c/review_s45/; certificates/theorem_c/review_s46/; src/p4/validated/{tc_layer.py,tc_layer_centre.py}; scripts/{s44_layer.py,s46_centre.py,s43_parts.py}; tests/{test_tc_layer.py,test_tc_layer_centre.py,test_s43_merge.py}; docs/review-2026-10-17-p4-layer-engine-adversarial.md; docs/review-2026-10-17-p4-s4-5-adversarial.md; docs/review-2026-10-17-p4-s4-6-centre-adversarial.md | 2026-10-17 |

| P4-10 | P4 | Theorem C (eigenvalue form) and the mode count: in the class C of P4-8 there is no eigenvalue in Ω = {Re κ ≥ 0} \ [0, 15] × [−14, 14] (P4-8 covers Re κ ≥ 6; P4-9 and conjugation cover the rest). With P4-6, the eigenvalues with Re κ ≥ 0 are exactly κ̄ ∈ [0.35569919907, 0.35569920835] (pure gauge) and κ₁ ∈ [2.81055254439, 2.81055255326], both simple zeros of E; so κ₁ is the unique non-gauge mode with Re κ ≥ 0 and γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231]. Corollary (Lemma G1): E(κ) ≠ 0 on Ω. Not claimed: kink modes, non-spherical modes, nonlinear stability, Re κ < 0 | CERTIFIED — 2026-10-17, with P4-9 (P4-6 and P4-8 CERTIFIED; Lemma G1 checked by the S4-5 review) | the P4-6, P4-8 and P4-9 evidence; problems/P4/notes/s4-5-assembly.md; docs/review-2026-10-17-p4-s4-5-adversarial.md; docs/review-2026-10-17-p4-s4-6-centre-adversarial.md | 2026-10-17 |

If the owner promotes without the move, the README's definition is not met,
and the rows should stay CONJECTURED with the text "certified pending the move
to certificates/".

**Related edits for the owner:**
- assembly §6 "Related edits": widen P4-2's scope note to Re κ ≥ 0 once P4-10
  is CERTIFIED, and record FORMULATION v1.5;
- `certificates/README.md`: move P4-9 and P4-10 out of "Not here yet";
- the nits of §11.
