# S4-4: the near-sonic layer engine (tau >= 224 on the band)

> **Status, 2026-10-17 morning (branch `p4/theorem-c-s4-3`; sections 9-13
> are new).**
> - **Band run** (section 9): all 140 band pieces certified in both
>   modes, 280 parts, no failure (commit `ac91f9fa`).  tau >= 224: worst
>   bound -0.750 on the sonic pieces, -3.51 on the gap steps (each at least
>   half of the step's own scale).  tau >= 14 (an independent second
>   certificate of the direct range): worst -0.452 (sonic/s29, after a
>   `--norm piece` re-run of s29-s50, whose default-engine bounds had gone
>   down to -8.5e-6), gap worst -0.604.
> - **Merge hook** (section 10): `s43_parts.py merge` joins layer parts to the
>   S4-3 parts by piece id; a valid passing layer part carries tau > 224 for a
>   pass-direct piece.  19 tests.
> - **Deep probe** (section 11): with a per-piece margin scale, derivative
>   blocks combined at the polynomial level and a centred (mean-value)
>   enclosure in (sigma, eps, 1/tau), the engine certifies all six probe steps
>   (tube/full-170, -200, -240, deep/deep-000, -064, -128) for tau >= 224 and
>   for tau >= 14, every certified bound at 50-54 % of the piece's own float
>   scale, in 1-4.6 min per piece (both modes, one pair build); also
>   full-138, a band-A step (full-175) and the widest step (full-244).
>   Without the centred form, tau >= 14 fails on full-240, deep-000 and
>   deep-064 (section 11.3: the obstruction).  Projected cost of the 236 tube/deep pieces of the
>   unlaunched September-style runs (tube 138-244, deep 0-128, 470-1070
>   core-h): about 8 core-h (section 11.5).
> - Not done: no layer part replaces an S4-3 part in the merge (that is a
>   semantic change for review, section 13), the centre (30 pieces) is
>   untouched, and CLAIMS.md is not edited.
>
> Earlier status (2026-10-17 night): engine built, 3 of the 140 band pieces
> certified (`sonic/s00`, `sonic/s25`, `gap/gap-088`, tau >= 224 bounds
> -0.76, -0.75, -0.85; the same three for tau >= 14).  Those six parts
> became stale on the S4-3 branch (their spec lacks the `tq` key added in
> `ed4cab14`, and the code hashes changed) and were recomputed by the band
> run.

Code: `src/p4/validated/tc_layer.py` (module), `src/p4/validated/tc_layer_mv.py`
(centred enclosures, section 11), `scripts/s44_layer.py` (driver),
`scripts/s43_parts.py` (merge, section 10), `tests/test_tc_layer.py`,
`tests/test_s43_merge.py`, `tests/test_s44_driver.py`.  Paths are relative to
`problems/P4/`.

## 1. What is bounded

The S4-3 condition (`s4-commutator.md` section 3): with kappa = eps + i tau,
c = -1/32,

    S = Re B1 - (tau/c) Im B1 - (f' + rho f)/2 + nu^2/4,
    A1 = A0 - (log B0)',  B1 = B0 + A0' - A0 (log B0)',
    nu = Im A1 + tau tt,  rho = Re A1 - c tt,  f = (c + tau^2/c) tt - (tau/c) nu,

and S < 0 pointwise on (-inf, 0) excludes an eigenvalue at kappa.  The open
S4-3 item is tau > 224 on the band xi = -x <= 0.0357: the 51 sonic pieces and
the gap-tube steps 0-88, 140 pieces.

A layer part certifies, for one piece P (an xi-interval), the statement

    for all xi in P, eps in [0, 6], tau >= 224:
        dT, dB, nB != 0   and   Sigma := S xi^2 / tau^2 <= b_P < 0.

Since xi > 0 and tau > 0, Sigma < 0 is S < 0.  The nonvanishing of dT, dB and
nB on the piece is part of the certificate (a denominator hull containing 0
rejects the box).  It makes the commutation legal (B0 = nB/dB != 0), the
eta-elimination legal (Bcr is a factor of dT), and keeps the coefficient
poles away, so the separate Rbar and real-pole checks of `analyze_piece` are
not needed for tau >= 224.  The boundary terms and the function-space argument
are kappa-independent and unchanged (`s4-commutator.md` section 3).  Combined
with the S4-3 direct phase of the same piece (tau in [14, 224]) the piece is
closed for tau >= 14; alternatively `--tau-min 14` runs the engine on
tau >= 14 directly (section 6).

## 2. The blocks of S

All quantities are ratios of kappa-polynomials whose coefficients are
Laurent Taylor models in the piece's TM variable (tc_tm `LC` = ell^off TM),
taken from `tc_commutator.build_pairs` exactly as S4-3 builds them (on sonic
pieces raised by `sonic_raise_inputs`, which rests on the exact identities
E1/E3/E5 of `parts/checks.json`).  With tt = TN_top/dT_top and
u2 = nB_top/dB_top:

    V  := A0 + kappa tt = nV/dT,       nV  = nT + kappa tt dT   (kappa^6 slot = 0),
    V' = nV'/dT - V dT'/dT,            nV' = nT' + kappa (tt dT)'  (kappa^6 slot = 0),
    Y  := (B0 - u2 kappa^2)/kappa = nY/(kappa dB),  nY = nB - u2 kappa^2 dB  (kappa^9 slot = 0),
    (log B0)'  = nB'/nB - dB'/dB,
    (log B0)'' = nB''/nB - (nB'/nB)^2 - dB''/dB + (dB'/dB)^2,
    W := A1 + kappa tt = V - (log B0)',   nu = Im W,  nu' = Im W',
    B1 = u2 kappa^2 + kappa (Y + X) + (V' - V (log B0)'),   X = tt (log B0)' - tt'.

The dropped top slots vanish identically by the definitions of tt and u2.
These two symbolic subtractions are the cancellations that defeat any
hull evaluation as tau -> infinity: A0 ~ -kappa tt and B0 ~ u2 kappa^2, while
nu, rho and Im B1 / tau stay bounded.  Everything else is value arithmetic.

One Laurent raise is new: on sonic pieces the Leibniz value of (tt dT_k)' has
offset -1, and it is raised to offset 0.  Justification: tt dT_k =
TN_top (dT_k/xi)/(dT_top/xi) is analytic at xi = 0 (every dT coefficient
carries the exact factor xi from Dp, E3, and the raised dT_top has a constant
term bounded away from 0, checked in `LayerData`), so its derivative has no
xi^{-1} Taylor coefficient.  `drop_lead` checks 0-containment and the
exact-coefficient range ve, as for the S4-3 raises.

## 3. The normalization (an identity, not an approximation)

Write lam = xi and multiply each block by lam to its natural order:
V^ = lam V, V'^ = lam^2 V', Y^ = lam^2 Y, L^ = lam (log B0)',
L'^ = lam^2 (log B0)'', t = lam tt, t_p = lam^2 tt', u2^ = lam u2, and
G = Y^ + t L^ - t_p,  H = V'^ - V^ L^,  W^ = V^ - L^,  W'^ = V'^ - L'^,
nu^ = Im W^,  nu'^ = Im W'^,  rho^ = Re W^ - (eps + c) t,  iota = 1/tau.  Then

    Sigma = u2^ lam (eps^2 iota^2 - 1) + eps iota^2 Re G - iota Im G + iota^2 Re H
          - (1/c) [2 eps u2^ lam + Re G + eps iota Im G + iota Im H]
          - (c iota^2 + 1/c) t_p/2 + (iota/c) nu'^/2
          - rho^ [(c iota^2 + 1/c) t - (iota/c) nu^]/2 + nu^2 iota^2/4.

Every kappa-power has been split exactly into its eps and tau parts.  No
term carries a positive power of tau.  So Sigma is bounded whenever the
normalized blocks are bounded, and these are bounded uniformly in z = xi kappa
in [0, infinity] (section 4).  As iota -> 0 at fixed xi, Sigma -> xi^2 K2 (the
S4-3 tau-quadratic hat).  As xi -> 0 at fixed s = xi tau it tends to the
layer function F(s, eps) of `s4-commutator.md` section 8.  The formula is
checked against the direct S4-3 chain at points (section 5).

## 4. Pivot forms: uniformity up to the corner xi -> 0, tau -> infinity

For P = sum_j a_j kappa^j with a_j = ell^{o_j} TM_j(ell), let
m = min_j (o_j - j) (the layer weight), p = max{j : o_j = min_k o_k} (the
fixed-kappa leading index), e_j = o_j - j - m >= 0 and g_j = o_j - m - p >= 0.
With z = ell kappa and w = 1/kappa = ell/z, exactly

    P = ell^m z^p Q,   Q = sum_{j >= p} TM_j ell^{e_j} z^{j-p} + sum_{j < p} TM_j ell^{g_j} w^{p-j},
    Q / z^{d-p} = sum_{j >= p} TM_j ell^{e_j} zeta^{d-j} + zeta^{d-p} sum_{j < p} TM_j ell^{g_j} w^{p-j},

zeta = 1/z.  Both brackets are polynomials in two small variables with TM
coefficients.  If the pivot lies on the layer line (e_p = 0), Q -> TM_p(0) at
(z, w, ell) -> 0.  If the top does (e_d = 0), Q/z^{d-p} -> TM_d(0) at
(zeta, w, ell) -> 0.  `LayerData` requires both for every denominator (dT,
kappa dB, dB, nB).  A ratio is N/D = ell^{mN-mD} z^{pN-pD} QN/QD (small-z
form) = ell^{mN-mD} z^{dN-dD} QN^L/QD^L (large-z form, |z| >= 4).  For all
eight blocks the ell-power cancels against lam^k exactly (k + mN - mD = 0 on
sonic pieces).  The only z-power left is z^1 in lam^2 nB''/nB (small form).
A negative ell-power on a strip touching 0 raises an error, and the box is
rejected.

Offsets on the closing piece (`result.forms` in the part):

| poly | offsets (kappa^0 .. top) | d | p | m |
|---|---|---|---|---|
| nT, nV, nV', dT' | all 0 | 6/5/5/5 | top | -d |
| dT | all 1 | 5 | 5 | -4 |
| dB, dB', dB'' | all 3 / 2 / 1 | 7 | 7 | -4/-5/-6 |
| nY | all 1 | 8 | 8 | -7 |
| **nB** | 1 (kappa^0..8), **2** (kappa^9) | 9 | **8** | -7 |
| nB' | 0 (kappa^0..8), 1 (kappa^9) | 9 | 8 | -8 |
| nB'' | all 0 | 9 | 9 | -9 |

All polynomials except nB and nB' are top-dominated: their kappa-expansion
at fixed xi is uniform in xi, so the w-form is uniform.  **The layer lives
only in nB.**  Its kappa^8 and kappa^9 terms share the layer weight, so
Q_nB = TM_9 z + TM_8 + sum_{j<8} TM_j w^{8-j}, whose layer symbol
6.85e4 z + 1.58e5 vanishes only at z = -2.30 = -u1^/u2^ (u1^ = -0.9099,
u2^ = -0.3951).  Since Re z = xi eps >= 0, that zero is never reached.  The
tau-quadratic phase expands 1/nB in powers of w at fixed xi, i.e. in powers
of 2.30/z.  That series diverges for |z| < 2.30, which is the layer, and this
is the precise reason it cannot close on the band.  The pivot form keeps
1/Q_nB exact.

## 5. Boxes, pseudo-boxes and the cover

A box (ell, eps, s, tau) in [la, lb] x [e0, e1] x [s0, s1] x [t0, t1] (s1, t1
may be infinite) stands for the points with s = xi tau.  `box_geometry`
intersects the ranges using that relation, with 1e-12 outward slack on the
float endpoints.  Every hull is then built by ball arithmetic: z = xi eps + i s,
w = 1/kappa, iota = 1/tau, lam = xi (on tube steps xi = -x_k + sigma, with
x_k an arb).
- tau = infinity (t1 = inf): w in {Re w in [0, e1/t0^2], Im w in [-1/t0, 0]},
  iota in [0, 1/t0].
- s = infinity (s1 = inf): zeta = u/(r + i), u in [0, 1/s0], r in [0, xi e1/s0].
- s -> 0 (the closing piece only): the z-box contains 0; the small form is
  valid there (Q_nB -> TM_8 != 0).
So the non-compact domain [lo, hi] x [0, 6] x [224, inf) is covered by
finitely many boxes.  Roots: the 12 S4-3 eps-cells; s geometric (ratio 2)
from 224 xi_lo (from 0 on the closing piece: [0, 1e-3] first) to
s_big = max(1e3, 1e4 xi_hi), then [s_big, inf]; tau in [224, inf].  On a piece
touching xi = 0, roots with s0 < 2e4 xi_hi are pre-split in tau (ratio 1.5
to 2e4, then [2e4, inf]).  The 1/tau terms of Sigma cancel in truth (there
Sigma is tau-independent to 1e-4) but not in interval arithmetic.  All
splits tile, so the cover is unchanged.  Refinement is greedy: each failing
box is bisected in whichever of s, eps, tau, ell gives the lowest worse
child.  A margin ladder certifies the first of Sigma <= -0.75, -0.5, -0.25
that succeeds within 20000 evaluations, else Sigma < 0 (400000 evaluations).

Soundness ledger: (i) the S4-3 inputs and their raises (E1/E3/E5); (ii) the
exact top-slot drops of section 2; (iii) the (tt dT_k)' raise of section 2;
(iv) exponents e_j, g_j >= 0 always, and ell-prefactor exponents checked;
(v) outward-rounded hulls; (vi) any denominator hull containing 0 gives a
non-finite Sigma, and the box is rejected.

## 6. Validation

**Float check.**  `point_S_direct` evaluates S by the direct S4-3 chain at a
point (constant TMs from the certified sonic series).  It reproduces the
September `corner3.py` values (`s4-commutator.md` section 8):

| s = xi tau (eps = 0) | 0.03 | 0.1 | 0.3 | 1 | 3 | 10 | 100 | 1e4 |
|---|---|---|---|---|---|---|---|---|
| September, xi = 1e-4 | -0.904 | -0.929 | -1.14 | -3.21 | -10.1 | -14.7 | -15.4 | -15.5 |
| now, xi = 1e-4 | -0.9036 | -0.9285 | -1.1440 | -3.2098 | -10.059 | -14.724 | -15.449 | -15.457 |
| now, xi = 1e-5 | -0.9006 | -0.9256 | -1.1411 | -3.2078 | -10.059 | -14.723 | -15.449 | -15.456 |

At eps = 3, s = 0.3, 1, 3: -40.891, -42.956, -49.803 (September -40.9,
-43.0, -49.8).  At fixed small xi (1e-6, 1e-7) Sigma is tau-independent
for tau in [14, 1e4]: -0.8978 (eps = 0), -1.2952 (0.03), -14.144 (1),
-80.377 (6), a slope of about -13.25 per unit eps.  So the band's smallest
true margin is 0.898, at eps = 0 and s -> 0.

**Layer formula against the direct chain.**  On tiny boxes on the closing
piece the engine gives, for example, s = 1: -3.2076228 +/- 7e-8 (direct chain
at xi = 1e-6: -3 +/- 0.47, at 1e-5: -3.2078).  It also gives s = 1e4:
-15.456150 and eps = 3, s = 1: -42.947235.  Note the radii: the direct chain
loses about 1/xi to cancellation near the sonic point, and the layer form does
not.  Random boxes with random interior points: 0 disjoint of about 40 on
`s25`.  On `gap-088`, 0 disjoint of 30, with the point values taken from the
sonic series, a background independent of the tube step.

**Tests** (`tests/test_tc_layer.py`).  The fast set (about 10 s) covers
pivot-form identities on synthetic Laurent polynomials, small/large-form
agreement, box-geometry containment, the September values, the layer algebra
against the direct chain at five (eps, tau) up to tau = 1e9, an end-to-end
cover at one xi for tau >= 224 and tau >= 14, and the part-file consistency.
With P4_SLOW=1: the three pieces end to end (174 s; the gap step needs
P4_GAP_TUBE or the cache).  [2026-10-17: the file now has 38 tests (about
20 s fast), adding the combined blocks, the piece scale and the centred forms
(section 12); `tests/test_s43_merge.py` (19) covers the merge hook and
`tests/test_s44_driver.py` (5) the driver's selection and skip rule.]

## 7. Results, certificate format and cost

| piece | xi range | tau >= 224: bound | boxes / evals | cover | tau >= 14: bound | boxes / evals | cover | build |
|---|---|---|---|---|---|---|---|---|
| sonic/s00 | (0, 1e-6] | -0.7633 | 862 / 2168 | 0.6 s | -0.5003 | 2048 / 11896 | 2.7 s | 53 s |
| sonic/s25 | [7.95e-5, 9.54e-5] | -0.7520 | 352 / 1388 | 0.2 s | -0.5001 | 2130 / 15276 | 16.8 s | 53 s |
| gap/gap-088 | [0.03475, 0.03569] | -0.8466 | 286 / 1532 | 0.3 s | -1.45e-4 | 29243 / 232852 | 19.4 s | 55 s |

"bound" is the largest certified upper bound of Sigma over the cover (the
certified margin is its negative).  All six parts come from commit
`f65c536c` (`git_dirty_src` false).  Without the tau pre-split of section 5,
`s00` only reached -0.257 at tau >= 224.  At tau >= 14 the ladder settles
at -0.5 on the sonic pieces.  The thin bound on `gap-088` comes from tau near
14, where Sigma is far from its layer limit.
Two more pieces were measured at tau >= 224 during design but not stored:
`s50` (xi in [0.00758, 0.0091], the widest Laurent domain) gave -0.783 with
351 boxes, and `gap-000` (xi = 0.0085) gave -0.760 with 343 boxes; both took
52-54 s to build and under 0.5 s to cover.

Part schema `s44-layer/1` (`parts/layer/<id>.json`, about 8 KB):
- `scope` (tau range, eps range, xi range, the statement) and `certified_bound`;
- `result`: ok, worst, boxes, evals, the ladder attempts, the depth histogram,
  up to 8 failing boxes, `params` (s-grid, ZSPLIT, tau_min, eps cells), and
  `forms` (the offset/pivot table of section 4);
- `inputs` (sonic series, or gap cache with content hash `a5428c96...`),
  `region_call`, the `code` hashes of every loaded p4 module, the git head,
  `env`, and `timing` (build/cover wall and cpu);
- `spec` equal to the S4-3 plan entry, so a merge can join the two parts by id.

**Cost.**  Per piece: one pair build, about 55 s on one M4 Pro core.  This is
the same `build_pairs` call that the S4-3 run of the piece makes.  The
cover takes 0.2-0.6 s for tau >= 224 (up to 3.5 s when the margin ladder
falls back) and 3-19 s for tau >= 14, plus about 8 s of context (Atoms) per
worker.  **All 140 band pieces: about 2.2 core-h**
(140 x 55 s plus cover), i.e. about 17 min on 8 workers or 65 min on 2 cores.
With `--tau-min 14`, about 2.5 core-h, and this replaces the S4-3 direct
phase on the band as well.  The September direct-only band runs took 16062 s
of wall time for the 51 sonic pieces and about 19000 s of direct phase for the
gap steps.  Fused into the S4-3 piece run, which has the pairs already, the
engine costs only the cover time.

## 8. Open items and risks

- **The other 137 band pieces** have not been run (instructed: no band run
  before review).  [2026-10-17 morning: all 140 run, section 9.]  Command: `python problems/P4/scripts/s44_layer.py work`
  (add `--tau-min 14` for the single-engine variant).  Five pieces spanning
  the band (both ends of the sonic ladder, its middle, both ends of the gap
  band) certify with margins of 0.75 or more at tau >= 224.  That is
  evidence, not a certificate, for the rest.
- **Merge.**  `s43_parts.py merge` counts a pass-direct piece only for
  tau <= 224.  It must learn to treat pass-direct plus a passing
  `parts/layer/<id>.json` (same `spec`, code hashes current) as covering
  tau >= 14.  Alternatively it can accept `parts/layer_t14` alone for the
  band.  Not done here, to keep `s43_parts.py` untouched for the merge with
  the S4-3 branch.  [2026-10-17 morning: done, section 10.]
- **Inputs.**  The certificate inherits the S4-3 inputs: `build_pairs`,
  `sonic_raise_inputs` (E1/E3/E5) and the tube steps.  The one new
  analytic input is the (tt dT_k)' raise of section 2.
- **Thin margins away from the sonic point.**  The normalization xi^2/tau^2 is
  the natural scale only in the layer.  On `gap-088` at tau >= 14, and in a
  scratch run on the deep step `deep/deep-064` (x = -3.75) at tau >= 224,
  the engine certifies, but with bounds of -1.5e-4 and -1.6e-6.  The
  deep-step cover took 32795 boxes and 22 s.  These are valid certificates.
  Using the engine to retire the tau-quadratic phase outside the band (as
  `s4-3-recovery.md` section 8 (ii) suggests) is plausible but untested
  beyond that one step.  A region-appropriate scale (e.g. dividing by
  |xi^2 K2|) would make its margins readable.  [2026-10-17 morning: the
  -1.6e-6 was a ladder artifact, not a thin truth (the float sup of Sigma on
  deep-064 is -0.148); with the piece scale, the combined blocks and the
  centred form the six probe steps certify at 50 % of their scale for
  tau >= 224 and tau >= 14; section 11.]
- **TSTAR coupling.**  The tau >= 224 parts assume the S4-3 direct phase uses
  tstar = 16 b1 = 224 (`TSTAR_FAC`).  If that changes, rerun with
  `--tau-min` equal to the new tstar.

## 9. The band run (2026-10-17)

All 140 band pieces (`sonic/s00`-`s50`, the 51 Laurent pieces; `gap/gap-000`-
`gap-088`, the 89 gap-tube steps with xi_hi <= 0.0357; membership is
`s44_layer.band_pieces`, the same set as the plan's direct-only pieces), in
both modes, one pair build per piece (`work --tau-mins 224,14`):

| pieces | engine | parts |
|---|---|---|
| sonic, 51 | the reviewed default: natural layer form, `norm=layer` (ladder -0.75 / -0.5 / -0.25 in Sigma units, then Sigma < 0 at 400000 evaluations) | `parts/layer/sonic/`, `parts/layer_t14/sonic/` |
| gap, 89 | `--combine --norm piece --mv mv3 --piece-ladder-evals 200000 --piece-final-evals 600000` (sections 11.2-11.4) | `parts/layer/gap/`, `parts/layer_t14/gap/` |

Why the gap steps use the extended engine.  They are tube steps (no Laurent
structure), so `combine` and the centred form apply.  With the default engine
`gap-088` certified tau >= 14 only at the final Sigma < 0 attempt
(-1.45e-4, 29243 boxes), and with `combine` alone `gap-000` gave -1.1e-5
after 110 s; with mv3 `gap-000` gives -0.604 (50 % of its scale 1.206) in
7.4 s.  At tau >= 224 the default ladder stops at its first rung
(-0.85 on gap-088) although the true margin there is ~14.5 (Sigma -> K2 xi^2
as s = xi tau grows); the piece ladder certifies half of it.  The
certified statement is the same in both engines.

**Results (2026-10-17 07:14Z; 280 parts, all pass, no error; commit
`ac91f9fa`).**  "bound" is the certified upper bound of Sigma; for the gap
steps the fraction of the step's own float scale is given too.

| pieces | tau >= 224 | tau >= 14 |
|---|---|---|
| sonic s00-s50 (default engine) | 51 pass; bounds -0.750 ... -0.783, all at the first rung -0.75; 351-862 boxes, cover 0.1-1.1 s | 51 pass; s00-s28 at -0.50, s29-s34 at -0.25, s35-s50 only at the final Sigma < 0 attempt (bounds -8.5e-6 (s40) ... -3.6e-4); 1598-15442 boxes, cover 2.4-27 s |
| gap 000-088 (combine, piece, mv3) | 89 pass; bounds -3.508 ... -9.513, every one >= 0.500 of its scale (scales 6.89 ... 14.52); 114-170 boxes, cover 0.4-0.7 s | 89 pass; bounds -0.604 (gap-000) ... -1.231, every one >= 0.500 of its scale (1.21 ... 2.46); 818-1419 boxes, cover 6.3-12.5 s |

Cost: 140 pair builds of 51-54 s plus about 0.4 core-h of covers, about
2.5 core-h in all, 83 min wall on two processes at nice 10 next to the six
S4-3 processes.  `sonic/s00` reproduces the reviewed parts bit for bit
(-0.7633164782585093 and -0.5003413981238266); `sonic/s25` gives -0.75197
and -0.50009 as before.  The six old parts were stale (spec without `tq`,
old code hashes) and are replaced.

The thin tau >= 14 bounds on s35-s50 are valid certificates: the default
ladder gives each rung only 20000 evaluations, and on the outer sonic pieces
the -0.25 rung needs more at tau near 14.  A re-run of s29-s50 with
`--tau-min 14 --norm piece` (same Laurent engine, ladder relative to the
piece scale, larger budgets) replaces those parts; results in the table
below.

Re-run of sonic s29-s50, `--tau-min 14 --norm piece` (07:23-07:41Z, two
processes; logs `logs/s44/s44-sonic-t14-piece-w0.log`, `-b-w0.log`): all 22
pass at the first rung, bounds -0.452 (s29) ... -0.585 (s50), each 0.50 of
the piece's scale (0.903 ... 1.170), 2515-51899 boxes, cover 5.7-62 s.  These
parts replace the default-engine ones in `parts/layer_t14/sonic/` (those are
kept in commit `ac91f9fa`).  So the band's tau >= 14 parts now certify
Sigma <= -0.45 or better everywhere: sonic s00-s28 at -0.50 (default
engine), s29-s50 at -0.452 ... -0.585, gap at -0.604 ... -1.231.

**Processes, checking, resuming.**  Launched 2026-10-17 05:51Z at nice 10,
two processes (the cap); both finished by 07:14Z (the s29-s50 re-run below
ran 07:23-07:41Z):
- `s44-band-gap-mv3` worker, PID 14909, log
  `logs/s44/s44-band-gap-mv3-w0.log` (launch record
  `logs/s44/launch-s44-band-gap-mv3.json`);
- a detached chain, `sh` PID 14911 (`logs/s44/launch-s44-chain.json`): the deep
  probe (done 06:03Z, `logs/s44/s44-probe-mv3.log`), then the sonic band
  (`logs/s44/s44-band-sonic.log`), then a second gap worker on the same lock
  queue (`logs/s44/s44-band-gap-mv3-w1.log`).

Check: `$PY problems/P4/scripts/s44_layer.py status` (tau >= 224) and
`status --tau-min 14`: both report 140 pass, 0 fail, 0 error, 0 stale,
0 missing.  Nothing is left to resume.  To recompute after a code change
(valid parts are skipped; a dead worker's lock is reclaimed only under the
same host name, which changed from `Seans-Mac-mini` to `Mac` at 06:03Z, so
older locks of dead workers go by hand from `logs/s44/locks/`):

    $PY problems/P4/scripts/s44_layer.py launch --workers 1 --nice 10 --tag s44-band-sonic --regions sonic --tau-min 224
    $PY problems/P4/scripts/s44_layer.py launch --workers 1 --nice 10 --tag s44-band-sonic-t14 --regions sonic --ids 'sonic/s0?,sonic/s1?,sonic/s2[0-8]' --tau-min 14
    $PY problems/P4/scripts/s44_layer.py launch --workers 1 --nice 10 --tag s44-band-sonic-t14p --regions sonic --ids 'sonic/s29,sonic/s3?,sonic/s4?,sonic/s50' --tau-min 14 --norm piece
    $PY problems/P4/scripts/s44_layer.py launch --workers 1 --nice 10 --tag s44-band-gap-mv3 --regions gap --tau-mins 224,14 --combine --norm piece --mv mv3 --piece-ladder-evals 200000 --piece-final-evals 600000

(the engine options are part of the skip rule, so a sonic tau >= 14 run with
the default options would replace the s29-s50 piece-ladder parts).

Commit the parts with `git add problems/P4/results/theorem_c/parts/`.

## 10. The merge hook (`s43_parts.py merge`)

`merge` now reads every `parts/layer*/<id>.json` and joins it to the plan
piece of the same id (`layer_join`).  `layer_state(p, d, tstar, eps_max)`
returns pass only if all of the following hold, and otherwise missing,
invalid, stale, fail or error:
- schema `s44-layer/1`, `id` and `spec` equal to the plan entry (the spec
  carries `tq`, the cache hash and the step);
- for tube pieces, the plan's cache content hash in the cache record and in
  the part's `inputs`;
- a non-empty `code.files` record, every hash current (or accepted with
  `--accept-code`);
- `status` pass, `result.ok`, a finite negative `certified_bound` equal to
  the result's;
- `scope.eps = [0, eps_max]`, `scope.tau = [tau_min, "inf"]` with tau_min equal
  to the cover's `params.tau_min` and tau_min <= tstar of the piece (224);
- `scope.xi` equal to the piece's xi-range.

Coverage: a piece counts for tau > tstar if it is pass-both, or pass-direct
with a passing layer part (`parts/layer` preferred, else the passing part with
the lowest tau_min, e.g. `layer_t14`).  Layer parts never stand in for a
missing or failed S4-3 part, and they do not change what carries tau > tstar
on a pass-both piece; without layer parts the output is as before (tested).
The used layer parts' source hashes enter `code_conflicts` like any other
used part.  New output fields: per piece `tau_gt_tstar` ("tau-quadratic",
"layer:<dir>" or "open"), `layer` (state, bound, tau_min, engine and
relative margin per directory) and `direct_range_second_method` (a passing
layer part with tau_min <= b1, i.e. a second certificate of [14, 224]); the
summary `layer` (states per directory, the non-passing ids, worst bound per
directory, counts); `scope.tau_gt_tstar_by`; and `strip_exclusion_certified`
= complete and no tau > tstar gap.  `ok`/`complete` keep their meaning (every
planned S4-3 piece present and passing, checks valid, no unaccepted code).
Tests: `tests/test_s43_merge.py` (synthetic plan and parts; 19 cases,
including twelve rejection modes).

**On the real parts (2026-10-17 07:47Z; the partial certificate is not
committed, `merge` regenerates it):** status partial, because the S4-3 run is
still going (195 pieces used, 381 missing, 0 failed, no code conflicts).
The tau > tstar gaps equal the tau <= tstar gaps: every pass-direct band
piece present so far (77 gap steps; the S4-3 sonic pieces have not started)
has tau > 224 carried by its layer part, and none is open.
`tau_gt_tstar_by` = 118 tau-quadratic, 77 layer, 0 open;
`direct_range_second_method` = 77; all 149 joined layer parts pass in both
directories (140 band + 9 probe; the probe parts join tube/deep pieces whose
S4-3 parts are missing, so they carry nothing).  The remaining band pieces
are joined automatically as their S4-3 direct-only parts land.

## 11. Deep probe: away from the sonic point

### 11.1 The natural scale, and the -1.6e-6

Point values of Sigma = S xi^2/tau^2 at eps = 0 (float design data; the
mid-step of each piece):

| piece | x | tau = 14 | 30 | 100 | 224 | 1e3 .. inf |
|---|---|---|---|---|---|---|
| tube/full-170 | -1.564 | -4.170 | -4.211 | -4.240 | -4.242 | -4.243 |
| tube/full-200 | -2.134 | -1.651 | -1.621 | -1.623 | -1.623 | -1.623 |
| deep/deep-000 | -2.950 | -0.702 | -0.525 | -0.484 | -0.481 | -0.480 |
| deep/deep-064 | -3.750 | -0.565 | -0.238 | -0.156 | -0.150 | -0.148 |

Away from the sonic point Sigma is nearly tau-independent and tends to
K2(xi, eps) xi^2 (the tau-quadratic hat); S/tau^2 = Sigma/xi^2 decays like
e^{2x} (0.055 at x = -2.95, 0.0105 at -3.75).  So the fixed ladder
-0.75 / -0.5 / -0.25 of the band is infeasible on deep steps, every rung fails,
and the final Sigma < 0 attempt stops refining each box as soon as its upper
bound is below 0, so its recorded worst bound can be arbitrarily close to 0:
that is the -1.6e-6 of the scratch run, not a thin truth.  Option `norm="piece"`
(`tc_layer.run_layer_piece`, `piece_scale`): scale = minus the largest point
value of Sigma at eps = 0 over ell in {lo, mid, hi} and
tau in {tau_min 2^k} plus the tau -> infinity pseudo-point; the ladder is
-0.5 scale, -0.25 scale, then Sigma < 0.  The scale only sets targets; the
certificate is still the rigorous upper bound of Sigma.  Scales:

| piece | full-170 | full-200 | full-240 | deep-000 | deep-064 | deep-128 |
|---|---|---|---|---|---|---|
| scale (tau >= 224) | 4.168 | 1.598 | 0.481 | 0.476 | 0.1469 | 0.0433 |

### 11.2 tau >= 224 works

With the piece ladder the five steps run with the natural evaluation and
`combine` (`layer_probe/natural/layer/`; the run was stopped before deep-128)
certify at the 0.5 rung in 0.5-2.4 s of cover, and all six do with mv3
(section 11.4).
Without `combine`, deep-064 needed 32795 boxes (19 s) for Sigma < 0 and
90670 boxes (46 s) for -0.5 scale; with it, 1470 boxes (2.4 s).  So
`s4-3-recovery.md` 8 (ii) holds: the tau-quadratic phase is replaceable by
the layer engine on the deep band at negligible cost.  It is also cheap there
(about 4 s per piece in the bench log), so this alone saves nothing.

### 11.3 tau >= 14: the obstruction of the natural evaluation

The expensive part of the September-style runs is the direct phase on
[14, 224] (tube/full-200: 1928-3533 s per eps cell, 12 cells, in the running
bench; section 11.5).  The engine's tau >= 14 mode would replace both phases.
With the natural interval evaluation it fails:

| engine | full-170 | full-200 | full-240 | deep-000 | deep-064 |
|---|---|---|---|---|---|
| default (no combine) | not run | fail at 2M evals | not run | not run | fail at 2M evals |
| combine, natural | -2.048 (0.50), 12.9 s | -6.0e-6, 105 s | fail, 2M evals, 222 s | fail, 2M evals, 221 s | fail at 1M evals |

(certified bound, relative margin in parentheses, cover time; fails: no
cover within the evaluation budget, boxes still positive at depth 20-21.)
The failing boxes are tiny and the truth there is far from 0, so this is not
a sign problem but a dependency loss, first order in the box size:
- full-200, box sigma in [0.018125, 0.0184375], eps in [2.15, 2.2],
  tau in [14.80, 14.84]: the point value is -19.17, and the default
  evaluation's upper bound is +2.07 at depth 20.
  (log B0)'' = nB''/nB - (nB'/nB)^2 - dB''/dB + (dB'/dB)^2 is a difference
  of separately hulled ratios of size 80, 81, 47, 49 with an O(1) result, so
  Im W' = Im(V' - (log B0)'') gets the hull [-12.25, 12.69] (true 0.22),
  weighted by (1/tau)/(2c) = -16/tau in Sigma.  `combine` (section 12)
  removes this one.
- deep-064 with `combine`, box sigma in [0.00957, 0.00977], eps in
  [4.3, 4.45], tau in [17.53, 17.62]: point value -4.63, hull [-9.79, +0.53].
  H = V' - V (log B0)' is ~ -85 with Im H = 0.036; its hull has
  Im H in [-2.24, 2.31], weighted by -(1/tau)/c = 32/tau.  The ratios V and
  (log B0)' are ~11 and ~7.5 with ~1 % hulls on this box, and their
  product's imaginary part does not see that their true imaginary parts are
  O(1/tau).  Thinning one direction at a time: eps-thin [-7.36, -1.90],
  tau-thin [-6.96, -2.30], sigma-thin [-8.59, -0.67].
Bisection halves such widths, so the cover needs depth ~20 at tau ~ 15-35 and
eps >= 2 and runs out of budget.  (On the sonic pieces the layer
normalization keeps the blocks O(1) and the default engine certifies
tau >= 14 at -0.5; the gap steps are in between: certified, but only at the
final Sigma < 0 attempt, section 9.)

### 11.4 The remedy: a centred enclosure (`tc_layer_mv`), and the results

Since the true Sigma varies little across such boxes, a mean-value (centred)
form removes the first-order loss.  `tc_layer_mv.sigma_box_mv3` (section 12
for the soundness argument) encloses Sigma on a tube-step box by its value
at the box centre plus forward-mode gradient hulls in (sigma, eps,
iota = 1/tau) times the half-widths, intersected with the natural enclosure.
On 60 random deep-064 boxes the median hull radius is 8.29 (natural), 2.80
(centred in eps, iota only: `sigma_box_mv`) and 0.63 (mv3), with no point
value outside (240 points).  deep-064, tau >= 14, Sigma < 0: natural fails
at 2M evaluations; mv certifies with 93033 boxes in 403 s; mv3 with 9123
boxes in 90 s.  One mv3 evaluation costs ~1.2 ms (~14 natural ones).

Probe results with `--combine --norm piece --mv mv3` (parts in
`parts/layer/{tube,deep}/` and `parts/layer_t14/{tube,deep}/`; sources of
commit `9a99f51c`, recorded head `34d297f7` with the same sources,
`git_dirty_src` false); certified bound, its fraction of the piece scale,
boxes, cover wall time (one M4 Pro core, nice 10, 8-10 other processes on 12
cores); the pair build is 53-54 s per piece, shared by the two modes:

| piece | x range | tau >= 224 | tau >= 14 | piece total (both modes) |
|---|---|---|---|---|
| tube/full-170 | [-1.584, -1.564] | -2.148 (0.52), 99 boxes, 0.4 s | -2.051 (0.50), 669 boxes, 5.4 s | 59 s |
| tube/full-200 | [-2.154, -2.134] | -0.807 (0.51), 123, 0.8 s | -0.798 (0.50), 1626, 14.4 s | 69 s |
| tube/full-240 | [-2.954, -2.934] | -0.258 (0.54), 159, 1.4 s | -0.241 (0.50), 5137, 46.8 s | 102 s |
| deep/deep-000 | [-2.963, -2.950] | -0.239 (0.50), 116, 0.7 s | -0.238 (0.50), 3267, 30.2 s | 85 s |
| deep/deep-064 | [-3.763, -3.750] | -0.0737 (0.50), 157, 1.1 s | -0.0735 (0.50), 10498, 94.7 s | 149 s |
| deep/deep-128 | [-4.560, -4.550] | -0.0217 (0.50), 256, 2.3 s | -0.0217 (0.50), 24823, 222.7 s | 279 s |

Three more steps, run afterwards with the production budgets
(`--piece-ladder-evals 400000 --piece-final-evals 1000000`): the first
fine step, a band-A step (where the September direct phase failed at
tau ~ 14, eps ~ 0 under the old depth cap) and the widest step:

| piece | x range, h | tau >= 224 | tau >= 14 | piece total |
|---|---|---|---|---|
| tube/full-138 | [-1.354, -1.349], 0.005 | -5.932 (0.95), 84 boxes, 0.3 s | -3.379 (0.55), 221, 1.2 s | 58 s |
| tube/full-175 | [-1.674, -1.664], 0.01 | -2.530 (0.71), 84, 0.2 s | -1.750 (0.50), 433, 3.2 s | 58 s |
| tube/full-244 | [-3.05, -3.0], 0.05 | -0.210 (0.50), 398, 4.1 s | -0.209 (0.50), 14218, 135.4 s | 194 s |

Every probe part passed at the first rung (0.5 of the scale); deep-128 needed
197660 of the 200000 evaluations allowed for that rung, so a production run
should allow 400000.  The natural-evaluation probe parts (`combine`, no
centred form: five passes at tau >= 224; at tau >= 14 two passes, one of them
at -6.0e-6, and two failures) are kept in
`results/theorem_c/layer_probe/natural/` for the record; they are not read by
the merge or the status.

### 11.5 Cost of replacing the September-style runs

The tau >= 14 cover time grows roughly like e^{-1.6 x} on h = 0.02 tube steps
(5.4 s, 14.4 s, 46.8 s at x = -1.56, -2.13, -2.93) and like e^{-1.25 x} on
the h = 0.0125/0.01 deep steps (30, 95, 223 s at x = -2.95, -3.75, -4.55).
Integrated over the plan (per piece one 54 s build, tau >= 224 cover ~1-2 s):
- tube idx 138-244 (107 steps, x in [-3.05, -1.349]: 23 steps h = 0.005
  and 11 steps h = 0.01 in [-1.69, -1.349], 71 steps h = 0.02, one
  h = 0.0059 step and idx 244 with h = 0.05): about 107 x 54 s of builds plus
  ~1500 s of covers (cover time scaled with h) = ~2.1 core-h;
- deep idx 0-128 (129 steps): about 129 x 54 s + 129 x ~96 s = ~5.5 core-h;
so **~8 core-h for the 236 pieces (range 6-12 with the budget margin), about
1 h on 8 workers**, against 470-1070 core-h estimated for the September-style
both-phase runs (`s4-3-recovery.md` section 4: band A + band B + deep
207 pieces 470-1070, tube 138-166 1-5).  Nine of the 236 pieces are now
certified (the probe parts), including idx 244 (h = 0.05, the widest step,
135 s at tau >= 14) and a band-A step (idx 175, 3.2 s); the other 227 are
the estimate.  For comparison, the S4-3 bench of the committed both-phase
path on the same steps (`results/theorem_c/bench/`, still running): on
tube/full-200 the direct phase took 1928, 2658, 2125, 3533 and 2789 s for
the first five of its 12 eps cells (~8.7 h projected), and on deep/deep-064
the first eps cell had not finished after 3.8 h (log at 07:28Z).  The layer engine did both phases of these two
steps in 69 s and 149 s.

## 12. Soundness ledger of the additions

- **combine** (`LayerData(combine=True)`, tube steps only): the polynomials
  NVp = nV' dT - nV dT', NL = nB' dB - nB dB', Mn = nB'' nB - nB'^2,
  Md = dB'' dB - dB'^2 are formed by TM arithmetic (`tc_commutator.kp_mul`,
  `kp_sub`) from the same jet slots, and V' = NVp/dT^2,
  (log B0)' = NL/(nB dB), (log B0)'' = Mn/nB^2 - Md/dB^2 are identities of
  the true functions.  `ratio_hull_w` evaluates N/(D1...Dr) =
  w^{sum d_Di - d_N} Q_N/(Q_D1...Q_Dr) for pure w-forms (all offsets 0, pivot
  = top; anything else raises and rejects the box).  Denominators: the same
  dT, dB, nB.  Tests: the identity on synthetic forms; Sigma with combined
  blocks against the independent direct chain at six (eps, tau) up to
  tau = 1e9.
- **mv** (`sigma_box_mv`): for each fixed sigma in the strip and fixed
  coefficient values, Sigma is a real C^1 function of (eps, iota) on the
  closed convex box E x I, I = [1/tau_hi, 1/tau_lo] (iota = 0 included:
  w = iota/(eps iota + i)); the mean value theorem along the segment from the
  centre gives Sigma in Sigma(c) + G_eps (E - eps_c) + G_iota (I - iota_c),
  G the gradient hulls over the box by forward-mode differentiation in ball
  arithmetic (`DV`: every operation encloses value and derivatives for all
  points the inputs enclose), dw/deps = -w^2, dw/diota = i q^2,
  q = 1/(eps iota + i).
- **mv3** (`sigma_box_mv3`): a TM f(sigma) in sum [c_k] sigma^k on [0, R]
  equals P + delta with P the exact midpoint polynomial and
  |delta(sigma)| <= r = sum rad(c_k) R^k (sigma >= 0).  For fixed
  delta* = delta(sigma) the map (s, e, i) -> Phi(P(s) + delta*, e, i) is C^1
  on the box, hence Phi(f(sigma), eps, iota) lies in
  Phi(P(sigma_c) + [-r, r], eps_c, iota_c) + G . (box - centre), with the
  coefficient values P(S) + [-r, r] and their derivatives P'(S) in the
  gradient pass; lam = xi_k + sigma has derivative 1, w none.  Applied to
  every coefficient TM of the forms used and to tt, tt', u2.
- Both centred forms are intersected with the natural enclosure (never
  wider); a denominator hull containing 0 in the gradient pass makes it
  non-finite, and then the natural enclosure is returned (itself non-finite
  there, so the box is rejected as before).  Sonic (Laurent) pieces always
  get the natural enclosure.  Disjoint enclosures raise (never observed).
- Tests (`tests/test_tc_layer.py`): `DV` against finite differences of an
  exact rational function, the midpoint split against `TM.eval`, containment
  of point values on random boxes and no-wider for both centred forms, an
  end-to-end cover with mv3; P4_SLOW: mv3 containment on a real gap step.
  Scratch check on deep-064: 60 boxes x 4 points, 0 outside, for mv and mv3.
- `norm="piece"` changes targets only.  The default path (`norm="layer"`,
  `combine=False`, no `sigma_fn`) is unchanged: gap-088 reproduces -0.846557
  (tau >= 224) and -1.4514e-4 (tau >= 14).

## 13. Recommendation: finishing the strip exclusion at least cost

1. **Band, tau > 224:** done by the band run once every band part passes
   (section 9); `merge` then reports `tau_gt_tstar_by.layer` = 140 and no
   tau > tstar gap on the band.  The S4-3 direct-only runs of the same 140
   pieces (running) still certify [14, 224]; the layer_t14 parts are an
   independent second certificate of that range.
2. **Tube 138-244 and deep 0-128 (236 pieces; not launched):** do not
   launch the September-style both-phase runs (470-1070 core-h).  Run the
   layer engine instead:
   `s44_layer.py launch --workers 8 --nice 10 --beyond-band --regions tube,deep --ids 'tube/full-13[89],tube/full-1[4-9]?,tube/full-2[0-3]?,tube/full-24[0-4],deep/*' --tau-mins 224,14 --combine --norm piece --mv mv3 --piece-ladder-evals 400000 --piece-final-evals 1000000`
   (about 8 core-h; `list` with the same selection gives the 236 pieces:
   tube idx 138-244 and the 129 deep steps; the six probe parts used the
   200000/600000 budgets, so with these options they are recomputed).  Then
   extend `merge`, after review, with an opt-in rule: a piece whose S4-3 part is missing may be covered for
   all tau >= 14 by a passing `layer_t14` part (tau_min <= b1), with the same
   validity checks as section 10.  The mathematical justification is the one
   of section 1: the part certifies dT, dB, nB != 0 and S < 0 pointwise on
   [14, inf) x [0, 6], which is what `analyze_piece` certifies in two phases
   plus its Rbar and real-pole checks (Bcr and D2r are factors of dT).  The
   widest step (idx 244) and a band-A step were probed and pass (section
   11.4).
3. **Centre (30 pieces):** now the dominant open cost (the committed path
   never finished: >= 25 CPU-h measured, 45-90 estimated) and the only region
   the layer engine does not cover.  The obstruction there is the same kind
   (K2 = O(y^2) against an O(1) tau-linear slack in the tau-quadratic phase,
   `s4-3-recovery.md` 11.5; the direct phase is slow).  The next step is a
   centre variant: pivot forms in the centre variable y with its Laurent
   offsets, the y^2 normalization, and the centred form of section 11.4
   extended to Laurent coefficients.  Until that exists, the centre stays on
   the committed default.
4. Keep the running S4-3 workers (tube 0-137, band, gap, sonic): they are
   cheap and the merge needs them.
