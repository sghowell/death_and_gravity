# S4-6: the centre variant of the layer engine (probe, 2026-10-17)

> **Status (branch `p4/theorem-c-centre`, from `d6ccd565`).**  Probe done.
> - New files only: `src/p4/validated/tc_layer_centre.py` (the centre
>   sigma-function), `scripts/s46_centre.py` (driver), `tests/test_tc_layer_centre.py`
>   (16 tests), outputs in `results/theorem_c/centre_probe/`.  No existing
>   file changed: the 1092 `layer_t14` parts keep their hashes and the
>   `s44_layer.py` pin.
> - Float grid: S < 0 with margin on c29 and c00, both tau ranges, as
>   y -> 0 and at eps = 0: Sigma_c = S/(1 + (y tau)^2) <= -6.10 on c29 at
>   tau >= 14 (-15.2 at tau >= 224), -6.0000 on c00.
> - Natural enclosure: **c29, tau >= 224: certified** (Sigma_c <= -1.53, 138
>   boxes, 2 s).  c29, tau >= 14: not closed in 20000 evaluations
>   (overestimation at small s = y tau with large eps; true margin -6 to
>   -20 there).  c00: not closed (near-zero leading Taylor coefficients of
>   nB and its derivatives dominate the pivot forms for y tau below ~1e-3).
> - Recommendation (section 6): go, conditional on two items: a mean-value
>   form for the graded evaluation, and an exact identity (or a lemma) for
>   the closing piece.  About 3-4 sessions and 2-5 core-h for the 30 pieces.

Paths are relative to `problems/P4/`.

## 1. The problem

The centre region is y = e^x in (0, e^{-4.5}], split into 30 plan pieces
(`centre/c00` = [0, 1.34e-8] closes at y = 0; `centre/c29` = [0.00694,
0.0111]).  The committed S4-3 path never finished there: the direct phase on
tau in [14, 224] measured >= 25 CPU-h without completing (45-90 core-h
estimated), and the tau-quadratic phase from tau = 14 fails on c29 by +77,
because its tau^2 coefficient K2 is O(y^2) against an O(1) tau-linear slack
(`s4-3-recovery.md` 11.5).  The layer engine (`s4-4-layer-engine.md`)
certified the 546 other pieces; its sonic normalization (Sigma = S xi^2/tau^2,
lam = xi) does not fit the centre.

## 2. What S looks like at the centre (float data)

Point evaluation of S by the direct chain (`tc_layer.point_S_direct`) on
point pairs from the certified centre series shows a layer in s = y tau:

    S = G(s, eps) + O(y),   G(0+, eps) = -6,   G(s, eps) ~ k2(eps) s^2 (s large),

with k2(0) = -19.0, k2(1) = -144.5, k2(6) = -772 (roughly -19 - 125.6 eps)
and G(1, 0) = -24.96 = -6 - 19.0.  The blocks tend to a constant-coefficient
limit: V = nV/dT -> 3, W = V - (log B0)' -> 1, rho -> 1, nu -> 0,
tt ~ 2.616 y^2, u2 ~ -1.962 y^2, Y = (B0 - u2 kappa^2)/kappa = O(y^2), so
B0 ~ u2 kappa^2 -> 0 as y -> 0 at fixed kappa, and Re H = Re(V' - V Lh) -> -6
is the s -> 0 limit of S.  The pairs have the Laurent offsets nT, dT, TN:
[9, 7, 7, 7, 7, 7, 9]; nB: [10, 8 x 7, 10, 10]; dB: [10, 8 x 7]; every Taylor
model is even in y; at y = 0, dT ~ kappa (kappa - 1)^4 and dB ~ kappa
(kappa - 1)^6 (binomial coefficients).

## 3. Design (`tc_layer_centre.py`)

**Layer data.**  `tc_layer.LayerData(pr, laurent=True)` with the TM variable
ell read as y builds the blocks (nV, dT, nVp, dT1, nY, kdB, dB, dB1, dB2, nB,
nB1, nB2, tt, tt', u2) from `tc_commutator.build_pairs` over
`tc_tm.centre_components(ce, hi)`, with `ce` the certified centre series at
the certified (a*, mu*) balls (`ec_centre.certified_centre_expansion(256,
K=30)`, as `tc_commutator._centre_expansion`).  Its checks pass: every
denominator (dT, kappa dB, dB, nB) has its pivot and its top on the layer
line.  Its derivative raise for (tt dT_k)' does not fire (the computed
offsets are at the target), and would be sound here too, since
d/dx = y d/dy keeps Laurent offsets.

**Pivot forms** (`tc_layer.PForm`, unchanged): P = sum_j y^{o_j} TM_j(y)
kappa^j = y^m z^p Q, z = y kappa, w = 1/kappa, with the small-z form
Q = sum_{j>=p} TM_j y^{e_j} z^{j-p} + sum_{j<p} TM_j y^{g_j} w^{p-j} and the
large-z form Q/z^{d-p} in zeta = 1/z.  A ratio N/D is y^{mN-mD} z^{pN-pD}
QN/QD; a negative z-power is rewritten exactly as y^{-1} w per power.

**Normalization.**  The certified quantity is Sigma_c = S/(1 + s^2),
s = y tau: O(1) and bounded away from 0 both as s -> 0 (-> G(0+) = -6) and as
s -> infinity (-> k2(eps)).  It has the sign of S.

**Graded evaluation.**  No single power lam = y^a normalizes every block
(V is O(1), tt and u2 are O(y^2), and the tau-linear terms of S are O(1) only
through tau y = s), so every quantity is a polynomial in one grading
variable g with exactly tracked integer exponents and ball coefficients
(class `GV`):

| | regime L (s >= 0.5) | regime S (s <= 0.5) |
|---|---|---|
| g | y | r = 1/tau |
| 1/tau | g/s | g |
| y | g | g s |
| z = y kappa | i s + g eps | i s + g s eps |
| w = 1/kappa | g zeta, zeta = (g eps - i s)/(s^2 + y^2 eps^2) | g (g eps - i)/(1 + r^2 eps^2) |
| lam | g^{-1} | g^{-1} |
| Sigma = S lam^2 tau^{-2} | S/s^2 | S |
| Sigma_c | Sigma * s^2/(1 + s^2) | Sigma/(1 + s^2) |

The assembly of S is tc_layer's identity S lam^2 tau^{-2} = T(lam) (the
formula of `sigma_box`, valid for every lam > 0) evaluated in this graded
arithmetic.  The identities used: z w = y, tau y = s, the two closed forms of
w and zeta above, and the pivot-form identities.  Regime L uses the large-z
forms (for nB, nB1, nB2, where z^2 dominates, the s-dependence then cancels
in QN^L/QD^L).

**Exact cancellation.**  lam = g^{-1} creates negative exponents.  They
cancel exactly, not by interval arithmetic: the exponent-0 coefficient of
every pivot form is real (only even powers of z appear without a y factor,
e.g. TM_9 z^2 in nB but TM_8 y z, and z = i s at g = 0), so the imaginary
parts that meet tau or lam are carried at exponent >= 1.  `GV.value` refuses
a box unless every negative-exponent coefficient of the final Sigma is
exactly zero (tested on boxes).

**Division.**  1/Q for Q = q0 (1 + u), u = O(g): the Neumann series to order
KMAX + 1 = 5 plus the exact remainder (-u)^{KMAX+2}/(1 + u), with 1 + u
bounded by the hull of u; a box whose q0 or 1 + u contains 0 is refused, so
dT, dB, nB != 0 on every certified box (y > 0, kappa != 0).  Exponents above
KMAX fold into KMAX + 1 by containment (g^a c = g^{K+1} g^{a-K-1} c).

**Centred ratios.**  At the centre N ~ c D for the log-derivative pairs
(nB ~ y^10, dB ~ y^8, dT ~ y^7, V -> 3), so N/D is evaluated as
c + (N - c D)/D with N - c D formed on the Taylor coefficients before any
hull: (nV, dT) 3, (nVp, dT) 21, (dT1, dT) 7, (dB1, dB) 8, (dB2, dB) 64,
(nB1, nB) 10, (nB2, nB) 100.  An exact identity for any constant c; it
removes the dependence on (y, s) that the separate hulls of N and D cannot
cancel.

**Cover.**  `tc_layer.run_layer_piece`, unchanged, with
`sigma_fn = sigma_box_centre`.  `tc_layer.box_geometry` with laurent=True
reads ell = xi = y, so its boxes are (y-strip) x (eps-cell) x (s = y tau) x
(tau).  Statement of a passing part: for all y in the piece (y > 0), eps in
[0, 6] and tau >= tau_min, dT, dB, nB != 0 and Sigma_c <= certified_bound < 0.

## 4. Probe results

All on `1c247e25` (clean tree), nice 10.

**Float grid** (`results/theorem_c/centre_probe/grid.json`, `s46_centre.py
grid`): 336 points per tau range on c29 (y in {lo, mid, hi}, 8 eps values,
tau = tau_min 2^k up to 2^39 tau_min), 560 on c00 (y from hi down to 1e-9 hi).

| piece | tau range | max Sigma_c | at | eps <= 0.03 |
|---|---|---|---|---|
| c29 | >= 224 | -15.18 | y = lo, eps = 0, tau = 224 | -15.18 |
| c29 | >= 14 | -6.10 | y = lo, eps = 0, tau = 14 | -6.10 |
| c00 | >= 224 | -6.0000000007 | y = 1.3e-14, s = 6e-6 | same |
| c00 | >= 14 | -6.0000000007 | y = 1.3e-11, s = 6e-6 | same |

So S < 0 with margin at least 6 (1 + s^2) in this design check, including
y -> 0 and eps = 0; the worst cell is eps = 0 at the smallest tau.  On c29
the graded sigma with the certified-ball Taylor models equals the direct
chain to 6e-14.  On c00 the certified-ball evaluation is undefined at 176
(tau >= 224) and 240 (tau >= 14) of the 560 points, all at s <= ~1e-3 (next
paragraph); the c00 values in the table use the Taylor-model midpoints (a
float check only).

**Natural enclosure** (`s46_centre.py run`, launch
`s46-centre-probe-20261017b`, parts in `results/theorem_c/centre_probe/t224/`
and `t14/`, ladder target -1 then 0, at most 20000 evaluations):

| piece | tau >= 224 | tau >= 14 |
|---|---|---|
| c29 | pass, Sigma_c <= -1.53, 138 boxes, 180 evals, 2.1 s | fail, 2520 boxes certified in 20021 evals, 283 s |
| c00 | fail at the corner, 274 boxes, 54 s | fail at the corner, 274 boxes, 53 s |

Build of the pairs: 33-54 s per piece.  Evaluation: about 3 ms per point,
10-25 ms per box.

**Why c29 fails at tau >= 14: overestimation.**  The failing boxes sit at
s = y tau in [0.097, 0.78] with eps >= 3.4 (tau in [14, ~110]), where the true
Sigma_c is -6 to -20.  The hulls there are far wider than the
function's variation over the box (15 times on the smallest failing boxes,
much more on the larger ones): the graded coefficients carry the box's
(y, s, eps) ranges independently, and the expansion parameter eps/tau is up
to 0.43 at tau = 14.  Finer root grids (s ratio 1.1, eps cells of 0.35) move
the failure but do not remove it within 25000 evaluations.  This is the
dependency problem a mean-value form solves (as `tc_layer_mv` did for the deep
tube steps); it is not a sign problem.

**Why c00 fails: the near-zero leading coefficients.**  The constant Taylor
coefficients of nB's kappa^1..kappa^7 coefficients (and nY's) are balls
around 0, radius 3e2-8e3 for nB, up to 2.3e5 (nB1) and 7.4e6 (nB2), against
~1e14-1e15 for the true leading terms; their true offset is 10, not 8
(every coefficient of nB is y^10 times an O(1) Taylor model, e.g. nB's
kappa^7 coefficient is -2.96e15 y^10 + ...).  The radii come from the widths
of the certified (a*, mu*) balls through cancellations, not from rounding.
With them, nB's pivot coefficient is TM_7 = [+/- 296] - 2.96e15 y^2 + ... and
the pivot form's leading coefficient TM_7 - TM_9 s^2 is dominated by the junk
for s <= ~1e-3 (y <= 1.34e-8): the hulls of fp and Im B1 grow like 1/s^2
(+/- 5 at s = 4.7e-4, +/- 115 at s = 1e-4) while the midpoint evaluation
gives Sigma_c = -6.0000 there.  Away from that corner c00 evaluates cleanly
(s in [0.01, 0.02]: [-6.18, -5.83]).

**Compute used:** about 0.9 core-h, all at nice 10 (exploration ~0.7 core-h,
of which three capped tau >= 14 attempts on c29 took 0.4; the committed runs
and grids 0.2).

## 5. What the closing piece needs

Yes, a lemma or an exact identity.  With the junk present, no cover closes
c00 at y tau -> 0: every box touching s = 0 has a pivot coefficient ball
containing 0.  Two routes:
- **An exact identity** (preferred, cheap): the y^8 layers of nB's
  kappa^1..kappa^7 coefficients and of nY's vanish identically, i.e. their
  true offset is 10, checked by exact computation on the centre series in the
  manner of E2 (`tc_commutator.check_layer_identity`).  With it the offsets
  can be raised (`lc_raise`), nB becomes a pure w-form y^10 kappa^9 [TM_9 +
  ...] with TM_9(0) = 9.9e13, and the pivot forms close at s -> 0.  The
  midpoint grid (-6.0000 down to y = 1e-17) suggests the identity holds.
- **A corner lemma**: S <= -6 + C (y tau)^2 + C' y for y tau <= s_c, from
  the constant-coefficient limit at x -> -infinity with explicit error
  bounds.  More prose, more review.

On the larger pieces the junk is a small relative perturbation: 1e-5 of
nB's pivot coefficient at y = 1e-4 (c20-c29 have y >= 1e-4), more for nB2,
whose centred difference nB2 - 100 nB is itself small.  Where it stops being
negligible has to be measured piece by piece; c29 is unaffected.

## 6. Recommendation and revised estimate

**Go, conditionally**, against the committed S4-3 default (45-90 core-h
estimated, never finished, and its tau-quadratic phase fails on c29 even
from tau = 14, so its tau >= 224 half is uncertain too).  The probe shows
the centre has a clean layer structure (S ~ -6 - 19 (y tau)^2 at eps = 0)
that the graded pivot forms capture exactly, and the asymptotic part
(c29 at tau >= 224) certifies in seconds.

Build:
1. Mean-value form for the graded sigma in (y, s, eps): forward-mode
   derivatives through the GV arithmetic (a DV class over GV coefficients),
   intersected with the natural enclosure, as `tc_layer_mv` does for tube
   steps.  1-2 sessions.  First milestone: c29 at tau >= 14.
2. The junk identity for nB, nB1, nB2, nY at the centre (exact check, then
   `lc_raise` in a new centre-data constructor), or the corner lemma.
   1 session.
3. A merge rule for centre parts (new code in `s43_parts.py`, outside this
   probe), the 30-piece run in both tau modes, and the review package.
   1 session.

Total: 3-4 sessions.  Compute: tau >= 224 is seconds per piece; tau >= 14
with a centred form is estimated at 1-5 min per piece, so 30 pieces x 2
modes = 1-3 core-h, budget 5 core-h.  The src side stays in new files (a
new module beside tc_layer_centre), so the 1092 parts' recorded hashes and
the s44_layer.py pin stay valid; the merge rule for centre parts is new code
in s43_parts.py, which no part pins.  If tc_layer.py itself had to change,
all 1092 parts would go stale: about 15 core-h to recompute, plus a
re-review.

Risks:
- The mean-value form may not tame the small-s, large-eps corner at tau
  near 14 (eps/tau up to 0.43).  Fallback: shrink the eps cells there, or a
  direct (non-graded) evaluation for tau in [14, ~50] at y >= 1e-4.
- The junk identity may not hold exactly.  Then pieces with y below ~1e-5
  need the corner lemma, a prose argument for review.
- Smaller-y pieces put more of tau >= 14 into regime S (s = y tau small), so
  per-piece cost may grow toward c00.
- Review: the graded arithmetic (exact cancellation, folding, remainder),
  the two regimes, and the centred ratios are new machinery, and need an
  adversarial review like the one of the layer engine.

## 7. Reproducing

    python problems/P4/scripts/s46_centre.py grid
    python problems/P4/scripts/s46_centre.py launch --pieces c29,c00 \
        --tau-mins 224,14 --max-evals 20000 --ladder-evals 5000
    python problems/P4/scripts/s46_centre.py status
    PYTHONPATH=problems/P4/src python -m pytest problems/P4/tests/test_tc_layer_centre.py
