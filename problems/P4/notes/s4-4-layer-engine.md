# S4-4: the near-sonic layer engine (tau >= 224 on the band)

> **Status, 2026-10-17 (branch `p4/s4-4-layer-engine`).**  Engine built and
> validated; **3 of the 140 band pieces certified**, nothing else.
> Certified, by stored parts produced from committed code:
> - `sonic/s00` (the closing piece xi in (0, 1e-6]), `sonic/s25`
>   (xi in [7.95e-5, 9.54e-5]) and `gap/gap-088` (xi in [0.03475, 0.03569],
>   the top band step): for every xi in the piece, eps in [0, 6] and
>   tau >= 224, the denominators dT, dB, nB are nonzero and
>   S xi^2/tau^2 <= -0.76, -0.75, -0.85 respectively.
>   Artifacts: `results/theorem_c/parts/layer/<id>.json`.
> - The same three pieces for tau >= 14 (both S4-3 phases at once):
>   `results/theorem_c/parts/layer_t14/<id>.json`.
>
> Not certified: the other 137 band pieces (cost estimate in section 7), and
> the merge of layer parts into `tc_commutator.json` (`s43_parts.py merge` does
> not read them yet; section 8).  CLAIMS.md is not touched.

Code: `src/p4/validated/tc_layer.py` (module), `scripts/s44_layer.py` (driver),
`tests/test_tc_layer.py`.  Paths are relative to `problems/P4/`.

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
P4_GAP_TUBE or the cache).

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
  before review).  Command: `python problems/P4/scripts/s44_layer.py work`
  (add `--tau-min 14` for the single-engine variant).  Five pieces spanning
  the band (both ends of the sonic ladder, its middle, both ends of the gap
  band) certify with margins of 0.75 or more at tau >= 224.  That is
  evidence, not a certificate, for the rest.
- **Merge.**  `s43_parts.py merge` counts a pass-direct piece only for
  tau <= 224.  It must learn to treat pass-direct plus a passing
  `parts/layer/<id>.json` (same `spec`, code hashes current) as covering
  tau >= 14.  Alternatively it can accept `parts/layer_t14` alone for the
  band.  Not done here, to keep `s43_parts.py` untouched for the merge with
  the S4-3 branch.
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
  |xi^2 K2|) would make its margins readable.
- **TSTAR coupling.**  The tau >= 224 parts assume the S4-3 direct phase uses
  tstar = 16 b1 = 224 (`TSTAR_FAC`).  If that changes, rerun with
  `--tau-min` equal to the new tstar.
