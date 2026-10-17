# S4-6: the centre of the Theorem C strip (x ≤ −4.5)

> **Status 2026-10-17: reviewed (no blocking or should-fix finding) and
> promoted: P4-9 and P4-10 are CERTIFIED, artifacts in `certificates/` (§10).**
> - All 30 centre pieces pass, in both τ modes (τ ≥ 14 and τ ≥ 224): 60
>   parts in `certificates/theorem_c/strip/centre_t14/` and `results/theorem_c/parts/centre_t224/`
>   (`c4d86e83`).  Worst bound at τ ≥ 14: Σ_c ≤ −5.91e−4 (c29).
> - The merge `s43_parts.py merge --substitute-layer-t14 --substitute-centre`
>   reports `status = complete` and **`strip_exclusion_certified = true`**:
>   576 pieces (546 `layer_t14` + 30 `centre_t14`), exact cover of x ∈ (−∞, 0)
>   with no gap, nothing missing, failed, rejected or blocked
>   (`certificates/theorem_c/strip/tc_commutator.json`, `a038640c`).
> - The closing piece c00 rests on an exact identity (§3.2), committed as an
>   artifact with provenance (`certificates/theorem_c/strip/centre_nb_identity.json`).
> - Neither checkpoint of the build plan triggered: c29 certified at τ ≥ 14,
>   and the exact identity for c00 holds.  No prose corner lemma is used.
> - Compute: 0.71 core-h for the 60 parts (4 workers at nice 10, 11.7 min
>   wall), well inside the 5 core-h budget.  Besides: the probe 0.9 core-h,
>   the point-sample check (§5.4) 0.33 core-h per run (run twice).
> - Promoted 2026-10-17 with the owner's approval (§10).  The 26
>   hash-pinned modules of `src/p4` and `scripts/s44_layer.py` (`3a86e50f`)
>   are unchanged, so the 1092 `layer_t14`/`layer` parts stay valid
>   (`git diff --stat fda58897..HEAD -- problems/P4/src` lists only the two
>   new files of §3).

Paths are relative to `problems/P4/`.  §8 is the reviewer's checklist.

## 1. The problem

The strip exclusion (ledger P4-9; `s4-5-assembly.md` §4.1, I8) needs, on an
exact cover of x ∈ (−∞, 0):

    S < 0 and dT, dB, nB ≠ 0  for every piece × ε ∈ [0, 6] × τ ∈ [14, ∞).

The S4-4 layer engine covers x ∈ [−4.56, 0] with 546 pieces.  The rest,
y = eˣ ∈ (0, YCUT], YCUT = e^{−4.5}(1 − 3.7e−17), is 30 plan pieces in y: `centre/c00` = [0, 1.337e−8]
(closes at y = 0) to `centre/c29` = [0.006943, 0.011109] (overlapping the
deep tube end at −4.56 by 0.06 in x).  The committed S4-3 path never
finished there (≥ 25 CPU-h for the direct phase; its τ-quadratic phase fails
on c29 by +77, `s4-3-recovery.md` 11.5), and the layer engine's sonic
normalization (Σ = S ξ²/τ², λ = ξ) does not fit the centre.

## 2. What S looks like at the centre

From the direct chain at points (`tc_layer.point_S_direct`):

    S = G(s, ε) + O(y),   s = y τ,   G(0+, ε) = −6,   G(s, ε) ~ k2(ε) s² (s → ∞),

k2(0) = −19.0, k2(6) = −772 (about −19 − 125.6 ε).  The blocks tend to a
constant-coefficient limit: V = nV/dT → 3, W → 1, ρ → 1, ν → 0,
tt ~ 2.616 y², u2 ~ −1.962 y², B0 ~ u2 κ² → 0.  Laurent offsets of the
pairs: nT, dT, TN [9, 7, 7, 7, 7, 7, 9]; nB [10, 8 × 7, 10, 10] as computed
(10 throughout after the raise of §3.2); dB [10, 8 × 7].  The certified
series have zero odd coefficients up to K = 30 (§3.2).

## 3. Design

Two new source files; nothing else in `src/` changed.

| file | role | sha256 (as recorded by the parts) |
|---|---|---|
| `src/p4/validated/tc_layer_centre.py` | the centre σ-function, the identity, the raise | `ec32ccd6…64ea` |
| `src/p4/validated/tc_layer_centre_mv.py` | optional mean-value form (§3.7), not used by the parts | not loaded |
| `scripts/s46_centre.py` | driver: identity, work, launch, status, grid | `735336cc…0db0` |

The module docstrings give the full derivations; this section is the
argument a reviewer needs.

### 3.1 Inputs and layer data

* The certified centre series at the certified (a*, μ*) balls:
  `ec_centre.certified_centre_expansion(256, K=30)` (a3 midpoint plus the
  Krawczyk offset; as `tc_commutator._centre_expansion`).
* `tc_tm.centre_components(ce, hi)`: Taylor models in y on [0, hi] with the
  Banach tail at order K + 1, A = 1 + y² Â, N = n/y, W = y² w, V = y v.
* `tc_commutator.build_pairs(AtomKP(at, comps, sd_shift=0))`, exactly as
  the S4-3 `centre_piece` builds its pairs: nT = −TN1r, dT = denT,
  nB = −U6, dB = denU, as κ-polynomial jets (value, x-derivative, second).
* `tc_layer.LayerData(pr, laurent=True)` (unchanged) builds the blocks nV,
  dT, nVp, dT1, nY, kdB, dB, dB1, dB2, nB, nB1, nB2, tt, tt′, u2 and checks
  that every denominator has its pivot and top on the layer line.
  `CentreLayerData` adds the centred numerators of §3.5.

### 3.2 The exact centre layer identity (the closing piece)

**The problem it solves.**  `build_pairs` computes nB's κ¹…κ⁷ coefficients
at Laurent offset 8, with constant Taylor coefficients that are balls
around 0 (radius 3e2–7e6 from the width of the (a*, μ*) balls) against
3e15–4e16 for the true y¹⁰ coefficients (nB's κ⁹ pivot is 9.9e13).  With
them no cover closes c00 at y τ → 0 (the probe, §9).

**Statement** (`tc_layer_centre.check_nb_identity`, exact `fmpq_mpoly`
arithmetic on the exact atoms `tc_energy.Atoms.mp`; about 10 s).  Grade a
monomial A^a N^b W^c V^e κ^j of U6 = Bcr·UN3r − Acr·EBc·D2r by the centre
layer p = −b + 2c + e (A = 1 + O(y²), N ~ y⁻¹, W ~ y², V ~ y), as E2 does
(`tc_commutator.check_layer_identity`).  Then:
* for κ¹…κ⁷ the lowest layer of U6 is 8, and the layer-8 polynomial, summed
  at A = 1, is identically zero in (n, w, v); the lowest nonzero layer is 10;
* κ⁰, κ⁸, κ⁹ start at layer ≥ 10;
* U6 has no odd-parity term (b + e odd), hence no odd layer.

**Why the y⁸ and y⁹ coefficients vanish.**  Write A = 1 + O(y²) and expand
n, w, v in y.  A monomial of layer p contributes y^p × (a series in y), so the
y⁸ coefficient of the κ^j coefficient comes from layer 8 at y = 0 only: the
layer-8 polynomial L₈ at A = 1, which is identically zero.  The y⁹ coefficient
is the derivative of L₈ ≡ 0 in the direction (n₁, w₁, v₁), plus odd layers, of
which there are none; so it vanishes on any background with A = 1 + O(y²).  On
the certified one n₁ = w₁ = v₁ = 0 anyway (zero odd coefficients up to K = 30,
`test_centre_series_is_even`), and the y⁹ TM coefficients are exact zeros, so
only the y⁸ drop rests on the identity: nB's κ¹…κ⁷ offset is truly 10.  The
x-derivative slots are y d/dy of slot 0: their y⁸, y⁹ coefficients (8c₀, 9c₁ and 64c₀, 81c₁) vanish too.

**Use** (`raise_nb`, then `centre_layer_data(..., raised=True)`).  Each of
nB's κ¹…κ⁷ coefficients, in all three jet slots, is re-expressed at offset
10 by `tc_tm.lc_raise(t, 2)` → `TM.drop_lead(2)`.  `drop_lead` refuses unless
2 ≤ ve (the dropped coefficients are exact Taylor coefficients, not
remainder-bearing) and the dropped balls contain 0; the identity supplies
their exact vanishing.  `raise_nb` refuses any other offset.  nY follows
(nY = nB − u2 κ² dB with u2·dB_{j−2} of offset ≥ 10).  On the raised data
every pivot form is a pure w-form (pivot = top), nB's pivot coefficient is
TM₉ (9.9e13 at y = 0), and c00 evaluates cleanly to −6 at the corner.

**Artifact.**  `certificates/theorem_c/strip/centre_nb_identity.json`, sha256
`0fad868ba62fdecb836e6f8a78de23bd3a87f096e50bd117b3cb0475c29bd505`, written
by `s46_centre.py identity` at `da7a7f28` (clean tree, driver `735336cc`,
the strict code record of `s43_parts.code_record`).  Every part's process
also re-runs `check_nb_identity` before building data, and refuses to work
if it fails; each part records the result (`identity` block).  The merge
checks the artifact (§6).  Consistent with the S4-5 review's
`certificates/theorem_c/review_s45/centre_layers.json` (val U6 = 10,
val denU = 8, no odd layers).  dB keeps offset 8: its layer 8 is nonzero
for κ¹…κ⁷ (tested).

### 3.3 Pivot forms and the normalization

The pivot forms are `tc_layer.PForm`, unchanged, in the Laurent variable y,
with z = y κ and w = 1/κ:

    P = y^m z^p [ Σ_{j≥p} TM_j y^{e_j} z^{j−p} + Σ_{j<p} TM_j y^{g_j} w^{p−j} ].

The certified quantity is

    Σ_c = S / (1 + s²),   s = y τ,

O(1) and bounded away from 0 as s → 0 (→ −6) and as s → ∞ (→ k2(ε)), with
the sign of S.  A passing part states: for all y ∈ [lo, hi] with y > 0,
ε ∈ [0, 6] and τ ≥ τ_min, dT, dB, nB ≠ 0 and Σ_c ≤ certified_bound < 0.

### 3.4 Graded evaluation

No single power λ = y^a normalizes every block at the centre (V is O(1),
tt and u2 are O(y²), and the τ-linear terms of S are O(1) only through
τ y = s).  So every quantity is a polynomial in one grading variable g with
exactly tracked integer exponents and ball coefficients (class `GV`):

| | regime L (s > 0) | regime S (s bounded) |
|---|---|---|
| g | y | r = 1/τ |
| 1/τ | g σ, σ = 1/s | g |
| y | g | g s |
| z = y κ | i s + g ε | i s + g s ε |
| w = 1/κ | g ζ, ζ = 1/z = σ(−i + y ε σ)/(1 + (y ε σ)²) | g (g ε − i)/(1 + r² ε²) |
| λ | g⁻¹ | g⁻¹ |
| Σ_c | (S/s²) · s²/(1 + s²) | S · 1/(1 + s²) |

S is assembled by tc_layer's identity S λ² τ⁻² = T(λ) (the formula of
`sigma_box`, valid for every λ > 0) in this arithmetic.  The identities used
are z w = y, τ y = s and the closed forms of w and ζ.  `sigma_graded` uses
the large-z forms in regime L and the small-z forms in regime S.  On the
raised data all 19 pivot forms (the blocks and the centred numerators) are
pure w-forms, pivot = top (checked on c00 and c29), so the two kinds
coincide and z never enters; κ enters only through w (regime L:
w = y ζ).  The comment on `Z_LARGE` in the module ("small-z forms for
s ≤ Z_LARGE") is stale: the constant only decides whether regime L builds
a z that nothing uses.  The comment cannot be fixed without changing the
module's hash, which the 60 parts record.

* **Exact cancellation.**  λ = g⁻¹ creates negative exponents.  The
  exponent-0 coefficient of a pure w-form is the real ball TM_p(y), and
  every imaginary part carries w = O(g), so the imaginary parts that meet
  τ or λ sit at exponent ≥ 1.  `GV.value` refuses a box unless every
  negative-exponent coefficient of the final Σ is exactly zero.  The
  structure is needed only for completeness: if it failed, boxes would be
  refused, never mis-certified.
* **Folding.**  Exponents above KMAX = 4 fold into KMAX + 1 with the g-hull:
  Σ_{a>K} g^a c_a = g^{K+1} Σ g^{a−K−1} c_a ⊂ g^{K+1} Σ gh^{a−K−1} c_a.
* **Division.**  1/Q for Q = q₀(1 + u), u = O(g): the Neumann series to
  order KMAX + 1 plus the exact remainder (−u)^{KMAX+2}/(1 + u), bounded with
  the hull of u.  A box whose q₀ or 1 + u hull contains 0 is refused.  Every
  regime evaluation divides by dT, κ dB, dB and nB, so a box that passes in
  any regime has all four nonzero (y > 0, τ < ∞).
* **Both regimes.**  `sigma_box_centre` evaluates every regime that applies
  (L if s_lo > 0, S if s_hi < ∞) and intersects the enclosures; disjoint
  enclosures raise (a bug).  Each regime encloses the same function over
  the same box, since its (g, s, ε) hull contains every box point.

### 3.5 Centred ratios

At the centre N ≈ c D for the log-derivative pairs, so N/D is evaluated as
c + (N − c D)/D, with N − c D formed coefficientwise on the Taylor models
before any hull (`CentreLayerData.CENTRED`): (nV, dT) 3, (nVp, dT) 21,
(dT1, dT) 7, (dB1, dB) 8, (dB2, dB) 64, (nB1, nB) 10, (nB2, nB) 100.  This
is an exact identity for any constant c; the constant only makes N − c D
small.

### 3.6 The cover

`tc_layer.run_layer_piece`, unchanged, with `sigma_fn = sigma_box_centre`,
no margin ladder (target 0) and at most 150000 evaluations.
`tc_layer.box_geometry` with laurent=True reads ℓ = ξ = y, so the boxes are
(y-strip) × (ε-cell) × (s = y τ) × (τ), over y ∈ [lo, hi], the 12 ε-cells of
[0, 6], s ≥ lo·τ_min (from s = 0 on c00) and τ ∈ [τ_min, ∞).  Unbounded
ranges enter through pseudo-boxes: τ = ∞ as g = 1/τ ∋ 0 (regime S) and
s = ∞ as σ = 1/s ∋ 0 (regime L).

### 3.7 The mean-value form (built, not used)

The plan made a mean-value form the first milestone for c29 at τ ≥ 14.  It
is in `tc_layer_centre_mv.py`: forward-mode derivatives (class DV) through
the graded arithmetic (class GVD), midpoint polynomial plus frozen
remainder δ as in `tc_layer_mv`'s mv3, intersected with the natural
enclosure.  Tests show it contains point values.  It was not needed: on
the identity-raised data, with both regimes intersected, the natural
enclosure certified c29 at τ ≥ 14 first (29924 evaluations), and the MV
form was tighter only sometimes, at about six times the cost.  The parts do
not load it, and the merge rejects a part that records it (§6).

## 4. Identities and where each is checked

| identity | kind | where |
|---|---|---|
| layer 8 of U6 = −nB vanishes for κ¹…κ⁷; no odd layers | X (exact fmpq) | `check_nb_identity`; artifact §3.2; `test_nb_layer_identity` |
| certified centre series even in y (odd coefficients exactly 0) | M | `test_centre_series_is_even` |
| nB = −U6, dB = denU, dT = Dp·SDr·D2r·Bcr | code | `tc_commutator.build_pairs` (pinned, S4-3/S4-4 reviewed) |
| pivot-form identities, small and large z, both regimes | algebra; tests | module docstring; `test_pivot_forms_graded` |
| regimes L and S give the same Σ_c | algebra; tests | `test_regimes_agree_and_normalization` |
| graded Σ_c = direct chain at points (c29, c00) | M (tests) | `test_agrees_with_direct_chain_c29/_c00` |
| negative exponents cancel exactly on boxes | M (tests) | `test_negative_exponents_cancel_exactly_on_boxes` |
| box enclosures contain point values | M (tests) | `test_box_enclosure_contains_points`, `test_mean_value_form_contains_points` |
| c + (N − cD)/D = N/D | exact | `CentreLayerData` |
| S λ² τ⁻² = T(λ) | X (sympy I6, S4-5 review) | `tc_layer.sigma_box` |

`tests/test_tc_layer_centre.py`: 22 tests, about 75 s (two piece fixtures
build the c29 and c00 data).

## 5. Results

### 5.1 The run

Launch `s46-centre-20261017` (`logs/s46/launch-s46-centre-20261017.json`,
machine-local): 4 detached workers at nice 10, `work --tau-mins 224,14
--max-evals 150000`, 19:09:50Z–19:21:30Z, at `e9cc5086` (clean tree).  Each
worker claims a piece (O_EXCL lock), builds its data once and runs both τ
modes.  All 60 parts record git head `e9cc5086`, `git_dirty_src = false`, a
passing write-time check, no late module and driver `735336cc`.  CPU: 1008 s
of builds (30 × 34 s) + 288 s (τ ≥ 224) + 1250 s (τ ≥ 14) = 0.71 core-h.

Aggregates (paths relative to the repository root, C-locale sort,
`find … -name '*.json' | LC_ALL=C sort | xargs shasum -a 256 | shasum -a 256`):
* `parts/centre_t14/` (now `strip/centre_t14/`, §10): `4ac8a7b00f9904d56dd5a26e7b656fd6387f1457614be7467dc6cad61e1f17d7`
* `parts/centre_t224/`: `ba424c3f78ea8d4a001072dd07025b7c6cd56acbfce12c5304ad14cf7f9fc836`

### 5.2 Per piece

Bounds are the worst upper bound of Σ_c over the accepted boxes.  The cover
stops refining a box once its upper bound is below the target 0, so a bound
is far above the function's supremum, which is about −6 on every piece
(§5.4).

| piece | y-range | τ ≥ 14: bound | boxes | evals | cpu s | τ ≥ 224: bound | boxes | cpu s |
|---|---|---|---|---|---|---|---|---|
| c00 | [0, 1.337e-08] | -0.184 | 505 | 680 | 12 | -0.184 | 433 | 11 |
| c01 | [1.337e-08, 2.14e-08] | -1.19 | 440 | 664 | 15 | -1.19 | 392 | 13 |
| c02 | [2.14e-08, 3.423e-08] | -1.29 | 430 | 668 | 14 | -1.29 | 382 | 13 |
| c03 | [3.423e-08, 5.477e-08] | -0.228 | 410 | 592 | 12 | -0.228 | 362 | 11 |
| c04 | [5.477e-08, 8.763e-08] | -1.21 | 416 | 640 | 14 | -1.21 | 368 | 13 |
| c05 | [8.763e-08, 1.402e-07] | -1.31 | 406 | 644 | 13 | -1.31 | 358 | 12 |
| c06 | [1.402e-07, 2.243e-07] | -0.27 | 385 | 560 | 12 | -0.27 | 337 | 11 |
| c07 | [2.243e-07, 3.59e-07] | -1.22 | 392 | 616 | 13 | -1.22 | 344 | 12 |
| c08 | [3.59e-07, 5.743e-07] | -0.12 | 381 | 612 | 13 | -0.12 | 333 | 12 |
| c09 | [5.743e-07, 9.189e-07] | -0.672 | 361 | 536 | 11 | -0.672 | 313 | 10 |
| c10 | [9.189e-07, 1.47e-06] | -1.24 | 368 | 592 | 13 | -1.24 | 320 | 12 |
| c11 | [1.47e-06, 2.352e-06] | -1.25 | 357 | 588 | 13 | -1.25 | 309 | 11 |
| c12 | [2.352e-06, 3.764e-06] | -0.241 | 336 | 504 | 11 | -0.241 | 288 | 9 |
| c13 | [3.764e-06, 6.022e-06] | -0.486 | 343 | 560 | 12 | -0.486 | 295 | 11 |
| c14 | [6.022e-06, 9.636e-06] | -1.33 | 333 | 564 | 12 | -1.33 | 285 | 11 |
| c15 | [9.636e-06, 1.542e-05] | -0.0248 | 311 | 472 | 10 | -0.0248 | 263 | 9 |
| c16 | [1.542e-05, 2.467e-05] | -1.29 | 319 | 536 | 12 | -1.29 | 271 | 10 |
| c17 | [2.467e-05, 3.947e-05] | -1.12 | 309 | 540 | 12 | -1.12 | 261 | 11 |
| c18 | [3.947e-05, 6.315e-05] | -0.35 | 287 | 448 | 9 | -0.35 | 239 | 8 |
| c19 | [6.315e-05, 0.000101] | -1.3 | 295 | 512 | 11 | -1.3 | 247 | 10 |
| c20 | [0.000101, 0.0001617] | -1.11 | 285 | 516 | 11 | -1.11 | 237 | 10 |
| c21 | [0.0001617, 0.0002587] | -0.0491 | 262 | 416 | 9 | -0.0491 | 214 | 7 |
| c22 | [0.0002587, 0.0004138] | -0.0108 | 270 | 480 | 10 | -0.0108 | 222 | 9 |
| c23 | [0.0004138, 0.0006621] | -0.702 | 261 | 492 | 11 | -0.702 | 213 | 9 |
| c24 | [0.0006621, 0.001059] | -0.109 | 252 | 504 | 10 | -0.385 | 192 | 7 |
| c25 | [0.001059, 0.001695] | -0.0754 | 261 | 576 | 13 | -0.0359 | 188 | 7 |
| c26 | [0.001695, 0.002712] | -0.0126 | 293 | 916 | 20 | -0.231 | 170 | 6 |
| c27 | [0.002712, 0.004339] | -0.00163 | 533 | 2920 | 64 | -0.624 | 154 | 5 |
| c28 | [0.004339, 0.006943] | -0.00163 | 1367 | 9592 | 207 | -0.184 | 151 | 4 |
| c29 | [0.006943, 0.01111] | -0.000591 | 3898 | 29924 | 653 | -0.228 | 137 | 4 |

c29 at τ ≥ 14 reproduces the milestone run of the build session bit for
bit (bound −5.913728371814791e−4, 3898 boxes, 29924 evaluations).  The
τ ≥ 224 mode is a second, independent cover of part of the same statement;
the merge reports it and does not use it.

### 5.3 The merge

`python problems/P4/scripts/s43_parts.py merge --substitute-layer-t14
--substitute-centre` at `2560a06d` → `certificates/theorem_c/strip/tc_commutator.json`,
sha256 `6c228a73abb42c2a2a17ef0c59f4e6565e5c4381b19b7d1dafc15556150fdbab`:
* `status = complete`, `ok = true`, `strip_exclusion_certified = true`;
* substituted 576 = 546 (`layer_t14`) + 30 (`centre_t14`); missing 0,
  failed 0, manual review 0, rejected 0, blocked 0;
* exact cover: 576 pieces, no error, no gap for τ ≤ tstar or τ > tstar;
  `open_at_centre_end = false`;
* `checks_valid_for_code = true`, no code conflicts, no accepted variants;
* centre identity artifact valid; `centre_t224` second mode: 30 pass;
* worst bounds by region: centre −5.91e−4 (c29), deep −2.17e−2 (deep-128),
  tube −0.209 (full-244), sonic −0.452 (s29), gap −0.603 (gap-000),
  band −3.95 (band-000).  The layer values equal those of the S4-5 record.

### 5.4 Point samples

`scripts/s46_sample.py` (a check, not part of the certificate) →
`results/theorem_c/centre_sampling.json`, sha256
`a740dd19abea81c74193caf7f3f788cad11f57f1597c4ddbe476c0382f789c06`, written
at `90cbd611` with a clean tree.  At
14400 random points (seed 20261017; per piece 6 values of y, log-uniform in
the piece and from 1e−12 on c00; 8 of ε including 0 and 6; 10 of τ
including 14 and 224, log-uniform up to 1e9), it compares three enclosures
of Σ_c:
* *graded*: the certificate's evaluator on a point box;
* *direct*: the direct chain of `tc_layer.point_S_direct` (no pivot forms,
  no grading, no raise) at the certified balls, rigorous but wide at small s;
* *direct_mid*: the same chain at the midpoints of the certified (a*, μ*)
  balls, a background inside the certified box, so every enclosure over the
  box must contain it.

Results: **0 failures** of the three checks (graded meets direct, graded
meets direct_mid, direct_mid below the piece's `centre_t14` bound and, for
τ ≥ 224, its `centre_t224` bound).  The graded value is defined at all 14400
points; the direct chain at the certified balls is finite at 13795 (it
divides by balls containing 0 at the smallest s); direct_mid is tight
(radius < 1e−6) at all 14400.  The largest |graded − direct_mid| is 1.1e−13;
the largest graded radius is 3.1e−6.  The largest sampled Σ_c is −6.0
(s → 0).  Per piece it is −6.112 on c29, −6.051 on c28, −6.017 on c27 and
−6.000000 to six digits on c00–c15 (y ≤ 1.5e−5), against certified bounds
between −1.33 and −5.9e−4.  The output keeps, per piece, the 20 points
nearest the bound; `--rows all` regenerates all of them (the run is seeded
and deterministic).  20 min wall, 0.33 core-h.

### 5.5 The plan's checkpoints

* (a) c29 at τ ≥ 14: certified in the first build session, by the natural
  enclosure on the raised data with both regimes intersected (the probe had
  failed it at 20000 evaluations on unraised data).  The mean-value form was
  built as planned (§3.7) but is not needed.
* (b) the exact identity for c00's near-zero nB coefficients: true (§3.2).
  Both checkpoints passed; no prose corner lemma is used anywhere.

## 6. The merge rule (`s43_parts.py merge --substitute-centre`)

Opt-in and strict, in the pattern of `--substitute-layer-t14` (`2560a06d`;
`centre_checks`, `centre_identity_checks`, `substitute_centre`).

* **Which pieces.**  A centre piece whose S4-3 part is missing, stale or
  invalid is carried for every τ ≥ b₁ by its passing `centre_t14` part.  A
  failed or errored S4-3 part paired with a passing centre part goes to
  `manual_review`, never to substitution; a usable S4-3 part is kept.
* **Per part** (`centre_checks`, always strict): schema `s46-centre/1`, id,
  region and the plan spec; a non-empty code record with every hash current
  (accepted variants do not count), no module loaded after the import-time
  snapshot, a passing write-time check, a recorded git head and
  `git_dirty_src` exactly false; the modules tc_layer, tc_commutator, tc_tm,
  tc_energy, tc_reduce, tc_layer_centre and ec_centre recorded, and the
  mean-value modules not; `driver_sha256` equal to the current
  `s46_centre.py`; status pass, ok, no failing boxes, a finite negative bound
  equal to the result's bound and worst; τ_min of the record, cover and scope
  equal and ≤ b₁, τ up to ∞; ε-cells tiling [0, EPS_MAX] and the scope's ε
  range; the cover's y-range [lo, hi] and the scope's; s_lo ≤ lo·τ_min and
  s_lo = 0 on the closing piece; laurent; combine off; evaluator
  `sigma_box_centre` in the cover and the engine record, by
  `tc_layer.run_layer_piece`, on raised data; the identity block equal to
  the raise of the current `tc_layer_centre` (2, κ¹…κ⁷); the certified
  centre series on the TM domain [0, hi].
* **Global** (each blocks every centre substitution): `checks.json` valid
  for the current code; b₁, c, EPS_MAX unchanged; the identity artifact
  valid (schema; ok; raise_by and κ list equal to the current raise; at every
  raised κ power no odd terms and lowest nonzero layer ≥ 10; the strict code
  conditions with tc_layer_centre, tc_energy and tc_commutator recorded and
  the current `s46_centre.py` as driver).
* **Cover.**  The substituted centre pieces join the exact cover with the
  other used and substituted pieces: [log lo, log hi] in exact arithmetic,
  consecutive centre pieces touching exactly, c00 from −∞.
* **Output.**  `centre_substitution` (substituted, rejected with reasons,
  blocked, manual review, no part, S4-3 parts kept, the identity artifact's
  hash and validity, the second τ mode's states) and a detail entry per
  piece (bound, boxes, evaluations, engine, code fingerprint, the
  nonvanishing statement, the second mode).
* **Wording fixes** (S4-5 review §8e): `scope.full_strip_tau` is replaced by
  `tau_by_carrier` (S4-3 pass-both [b₁, ∞), pass-direct [b₁, tstar], layer
  join (tstar, ∞), substituted `layer_t14` and `centre_t14` [b₁, ∞)) and
  `tau_target = [b₁, ∞)`; `xi_open_max` is given only for an open set
  bounded in x, and `open_at_centre_end` says whether a gap reaches x = −∞.
  The ignored `ncells` argument of `run_centre` stays: `tc_commutator.py` is
  hash-pinned.
* **Tests** (`tests/test_s43_merge.py`, 214 tests, about 11 s): one per
  rejection mode of a centre part (58, each pinned to its one reason) and of
  the identity artifact (20); the global blocks; positive cases (missing,
  stale and invalid S4-3 parts substituted, manual review, S4-3 part kept,
  both rules together, a two-piece centre with and without the closing
  piece, the second mode reported but not used, accepted code ignored); the
  scope wording; the constants shared with the driver; a record written by
  `s46_centre.write_part` passes `centre_checks`.

## 7. Reproducing

From the repository (worktree) root, Python 3.11 with python-flint 0.9.0:

    # the identity artifact (about 10 s)
    python problems/P4/scripts/s46_centre.py identity --out /tmp/centre_nb_identity.json
    # the 60 parts (about 0.7 core-h); resumable: skips parts valid for the
    # current code, so add --force to recompute (in a scratch checkout, below)
    python problems/P4/scripts/s46_centre.py launch --workers 4 --nice 10 \
        --tau-mins 224,14 --max-evals 150000 --tag replay
    python problems/P4/scripts/s46_centre.py status -v
    # the merge (seconds); rewrites certificates/theorem_c/strip/tc_commutator.json
    # (only the provenance's head and time change; check with git diff)
    python problems/P4/scripts/s43_parts.py merge --substitute-layer-t14 --substitute-centre
    # the point samples (about 20 min)
    python problems/P4/scripts/s46_sample.py --out /tmp/centre_sampling.json
    # tests
    PYTHONPATH=problems/P4/src python -m pytest problems/P4/tests/test_tc_layer_centre.py \
        problems/P4/tests/test_s43_merge.py

To recompute committed parts without touching them, work in a scratch
checkout (`git worktree add /tmp/s46 <head>`) and run `s46_centre.py work
--pieces c29,c00 --force` there; compare `certified_bound`, `boxes` and
`evals` (timing, env and dates differ).  The computation is deterministic:
the production c29 τ ≥ 14 part equals the earlier milestone run.

## 8. What an adversarial reviewer must check

1. **The identity and the raise** (§3.2): the grading argument (layer 8 at
   A = 1 and no odd layers ⇒ the y⁸ and y⁹ coefficients vanish); the
   evenness of the certified background; that nB = −U6 in `build_pairs`;
   that the jet slots are y d/dy images; that `drop_lead(2)` drops exact
   Taylor coefficients (k ≤ ve); that nY's offsets follow.
2. **The graded arithmetic** (`GV`): folding, the Neumann inverse and its
   remainder bound with the hull of u, the refusal of nonzero
   negative-exponent coefficients.
3. **The two regimes**: the closed forms of z, w, ζ, 1/τ and y, the
   normalization factors s²/(1 + s²) and 1/(1 + s²), the applicability
   conditions (L needs s > 0, S bounded s), and the soundness of the
   intersection.
4. **The centred ratios**: that the `C_N_D` forms are exactly N − c D.
5. **The cover**: that `run_layer_piece` with laurent=True covers
   [lo, hi] × [0, 6] × [τ_min, ∞) through s = y τ, from s = 0 on c00, and
   that every accepted box divides by dT, κ dB, dB and nB.
6. **Nonvanishing**: that dT, dB, nB ≠ 0 on y > 0 meets I8's requirement
   (with Bcr, D2r ≠ 0 from dT = Dp·SDr·D2r·Bcr ≠ 0).
7. **Inputs**: the certified (a*, μ*) balls, the TM domain [0, hi] with its
   Banach tail, and pairs built exactly as the S4-3 `centre_piece`.
8. **The merge rule and cover** (§6): the strict checks, the identity
   artifact check, the exact cover with c00 from −∞ and c29 overlapping the
   deep end.
9. **Provenance**: the 60 parts at `e9cc5086` with a clean tree and driver
   `735336cc`; the identity at `da7a7f28`; every recorded hash current; the
   26 pinned modules and `s44_layer.py` unchanged since `fda58897`.
10. **Independent evidence**: the point samples (§5.4) and a spot
    recompute of c29 (τ ≥ 14) and c00 in a scratch checkout.

The S8 centre boundary term (BT = O(t^{3+k}) → 0, `centre_layers.json`) is
already reviewed (S4-5 review §4) and unchanged; the centre parts supply
only the pointwise condition for y > 0.

## 9. The probe (history)

The probe (`8dc4e203`; outputs in `results/theorem_c/centre_probe/`,
superseded) ran the natural enclosure on unraised data, with one regime
chosen per box:
* c29 at τ ≥ 224 passed (Σ_c ≤ −1.53, 138 boxes); c29 at τ ≥ 14 failed in
  20000 evaluations: overestimation at s ∈ [0.1, 0.8] with ε ≥ 3.4, where the
  true Σ_c is −6 to −20;
* c00 failed at the corner y τ → 0: the junk balls of nB's leading
  coefficients dominated the pivot forms for y τ ≲ 1e−3, while the midpoint
  evaluation gave −6.0000 there;
* the float grid gave max Σ_c = −6.10 on c29 (τ ≥ 14) and −6.0000 on c00.
It recommended a mean-value form and an exact identity (or a corner lemma)
for the closing piece.  The build found that the identity, plus evaluating
and intersecting both regimes, closed both pieces; the probe's 0.9 core-h is
not included in §5.1.

## 10. Review and promotion (2026-10-17)

**Review.**  `docs/review-2026-10-17-p4-s4-6-centre-adversarial.md`
(`e17ff004`, `2bda78ba`): no blocking and no should-fix finding.  It
re-derived the raise identity exactly with independent code, recomputed all
60 parts identically, verified the covers of c00, c27, c28 and c29 exactly,
checked the x-cover of all 576 pieces independently, and attacked the merge
rule 50 ways.  Its checks are in `certificates/theorem_c/review_s46/`.

**Promotion.**  With the owner's approval, P4-9 and P4-10 are CERTIFIED
(CLAIMS.md), and the strip certificate moved unchanged (git renames at 100%
similarity) from `results/theorem_c/` to `certificates/theorem_c/strip/`:
`plan.json`, `caches/`, `layer_t14/`, `centre_t14/`, `centre_nb_identity.json`
and `tc_commutator.json`.  The review's checks moved to
`certificates/theorem_c/review_s46/`.  `centre_t224/`, `centre_sampling.json`
and the probe stay in `results/` (supporting, not cited).
* `s43_parts.py` (not pinned) reads the moved files at the new paths
  (`STRIP`, `layer_dir_path`).  The pinned `s46_centre.py` still writes new
  parts under `results/theorem_c/parts/centre_t<T>/` and the identity under
  `results/theorem_c/` (or `identity --out`), so a replay never overwrites
  the certificate.
* The merge re-run at the new paths reproduces `tc_commutator.json` except
  `provenance` (git_head, merged, driver_sha256 and the new strip_dir) and
  the identity's recorded path.  The committed file is the reviewed one
  (`6c228a73…`).
* Aggregates (the §5.1 command): `strip/centre_t14/` gives `7b087f9a…` at the
  new paths.  The reviewed `4ac8a7b0…` is the same content at the old paths
  (map the paths back to `results/theorem_c/parts/centre_t14` before the outer
  shasum and it reproduces).  `strip/layer_t14/`: `f21d32bd…` new,
  `e9835a41…` old.

**Review nits recorded here.**  Fixing them in code would change hashes that
the 60 parts (or all 1152 parts) record.
* `tc_layer_centre.Z_LARGE`'s comment is stale: regime L always uses the
  large-z (ζ) forms; `Z_LARGE` only decides whether `Grade` builds a `Z` that
  nothing uses.  `S_SPLIT`'s comment ("regime S below, regime L above") is
  stale too: production evaluates both regimes wherever they apply and
  intersects them, and `S_SPLIT` only steers `regime_of`, the default when
  `sigma_graded` is called without a regime (tests only).  So every part's
  `engine.s_split = 0.5` records no split that was used.
* `engine.raised`, `inputs.*` and `engine.s_split` are constants the driver
  writes, not measurements.  The merge checks them, but the guarantee that
  the data were raised and built from the certified balls is the pinned
  `s46_centre.py`, which hard-codes both (§6).
* `tc_tm.TM.inv` records `ve = m`, ignoring `self.ve` (pinned, pre-existing),
  so the recorded exact-coefficient range of an inverse can exceed the
  input's.  Not triggered here: every TM on the centre path is also capped by
  the component TMs (ve 31), and `drop_lead(2)` runs at ve = 25.
* §3.2's identity uses the exact `fmpq_mpoly` atoms (`at.mp`); evaluation
  reads the compiled term lists (`at.cmp`), the same polynomials.
* §5.4's samples evaluate the graded form on point boxes only.  The box
  enclosures are covered by the review's containment checks (its §4–§5).
* `tests/test_tc_commutator.py::test_centre_first_piece` is now gated behind
  `P4_SLOW` (it did not finish within 15 min).
