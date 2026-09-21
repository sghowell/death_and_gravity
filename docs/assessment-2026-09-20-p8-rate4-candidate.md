# RATE4 v1: a matter-calibrated, physically normalized EFT candidate

Date: 2026-09-20 (local). The user approved proceeding with the separately
named candidate proposed in the [parent-input comparison](assessment-2026-09-20-p8-match1-parent-input.md).
This document supplies its four defining rate conditions and constructs
their first-order finite matching. It does not change the frozen original
parent or claim to close original MATCH-1 or P8.

**Result:** the four conditions have a unique first-order solution in four
specified on-shell directions, for every fixed, finite choice of the
remaining inputs. The solution is independent of finite changes of the
hard subtraction reference. The complete retained Born-tree radiation
conversion is also bounded at this order, contributing less than
`10^-783*lambda` to the specified matching projection. Complete hard
matching, complementary curved data and higher-order physical error
bounds are not yet supplied.

## 1. Identity, retained data and new assumptions

The name is **QG2-H8A420-RATE4, v1 (matter-calibrated targets)**. Retain the
QG2-H8A420 classical action, spectrum, physical matter frame and original
parameters

```text
mu = 1,  kappa0 = kappa = 10^800,  lambda = 10^-600,
D = 10^200/512,  n = D+2,  g = 1/8192,
C_tree = -g^2*(3/D - 2/D^2).
```

The light pole mass/residue and heavy onepoint conditions remain as in
[S239](../problems/P8/s6/continuation/s6_239/FORMULATION.md).
Use the complete S239 H8A420-VAC-OS4 matter-loop amplitude `L_m`, including
its full-action jets and finite value condition, as the fixed reference.
[S297](../problems/P8/s6/continuation/s6_297/notes/source.md) identifies it
with the original finite-kappa one-loop contribution with no internal
graviton. The new order-hbar counterterms do not change those loop
diagrams at this order; inserting one into a loop costs another order.

This retains a **reference contribution**, not a fifth normalization
condition on the new total amplitude. In particular, the total
subthreshold symmetric value need not remain the old value after the new
finite gravitational matching. A scalar finite contact need not contain
an explicit graviton. Its disappearance in a quantum decoupling limit
must be demonstrated, not inferred from its bookkeeping label.

Likewise `kappa` is the classical/reference parameter here, not a newly
imposed measured Newton coupling. The finite `delta_kappa` direction is
one of the four unknowns. Adding an independent Newton normalization
would change the number of conditions and must not be done silently.
No new heavy pole/residue or cubic normalization is implicitly imposed.

The four rate targets below are **new model-defining physical conditions**,
not measurements, old-parent predictions, minimal-subtraction zeros or
naturalness priors. They specify a four-condition EFT candidate family;
they do not determine every complementary Wilson coefficient, quantum
state or curved/source response. All such remaining inputs stay explicit.

## 2. Exact targets, with the matter correction retained

At the four ordered nonforward points

```text
q1 = (8,-2,-2), q2 = (10,-3,-3),
q3 = (12,-4,-4), q4 = (16,-6,-6),
```

use the infrared-regulator-removed inclusive differential rate, normalized
by the full unexpanded reference Born rate at that point. Keep the
physical matter frame, total-unresolved-energy rule and detector
resolution `x_det=1/256` specified in the preceding comparison. The
regulator is removed at fixed resolution and fixed nonforward angle.
Pointwise theoretical normalizations are intended; finite bin averaging
would require its own correction. The inclusive rate is not itself a
crossing-analytic amplitude.

Let

```text
a_i = C_tree + g^2*sum_channels 1/(n-channel) - T(q_i)/kappa,
rho_m,i = 2 Re L_m(q_i)/a_i.
```

The defining targets are

```text
R_i^target = 1 + rho_m,i,       i=1,...,4.
```

Their formal perturbative specification is `1+hbar*rho_m,i`, with no
additional higher-order target coefficients. This is a proposed exact
normalization prescription; **only its first-order solution is constructed
here**. No existence or convergence of an all-order physical realization
is asserted. The defining `L_m` is the frozen complete amplitude, not a
rounded decimal or the interval bound below.

At these points the exact Born intervals are

```text
136 < a_1/lambda < 137,       228 < a_2/lambda < 229,
344 < a_3/lambda < 345,       648 < a_4/lambda < 649.
```

The original polynomial scalar tree is positive; both the exact heavy
exchange remainder and the gravitational Born contribution at these
points are positive. Thus `a_i` exceeds that polynomial tree. The complete
[S239 matter bound](../problems/P8/s6/continuation/s6_239/notes/matching.md)
then gives `|rho_m,i|<2*10^-199`. All four stipulated targets are positive.
This is a bound on the **chosen reference values**, not a bound on the
new candidate's missing quantum remainder.

Choosing `R_i=1` instead would cancel the retained matter correction at
the four centers as well. RATE4 v1 deliberately cancels only the additional
first-order contribution to these inclusive observables.

## 3. Constructive first-order matching

In mass-one units use the same four shapes as the preceding audit:

```text
S(q) = s^2+t^2+u^2,
H(q;n) = sum_channels (channel+2)/(n-channel),
T(q) = sum_cyclic (2-2a-bc)/a,
b(q) = [S(q), H(q;n), 1, T(q)],
M_ij = b_j(q_i),       theta = [A, B, C, D_N]^T.
```

Here the matching coordinate `C` is not the classical `C_tree`. Write the
complete first-order hard correction as

```text
L_m(q) + F(q) + b(q) theta.
```

`F` must include **all** remaining first-order hard contributions and any
fixed complementary matching data in one common reference. Let `D_i^rate`
be the complete finite conversion from that hard reference to the chosen
physical inclusive rate. It includes the hard/soft convention conversion,
the resolution change and the relevant real-emission/recoil terms; it is
not just the selected known soft term. The S302 analytic reference
`resolution/nu=1` must not be identified with `x_det=1/256`.

The first-order rate is

```text
R_i = 1 + hbar*[rho_m,i
       + 2*(Re F(q_i)+(M theta)_i)/a_i + D_i^rate] + O(hbar^2).
```

Finite insertions in loop or real-radiation corrections enter the next
order in this bookkeeping. A complete higher-order error estimate is not
provided by the notation `O(hbar^2)`.

Set `d_i=Re F(q_i)+a_i D_i^rate/2`. The solution is exactly

```text
theta_RATE4 = -M^(-1) d.
```

There is no zero default for either of those inputs. The determinant is

```text
det M = -2 P(n)/[5(n-16)(n-12)(n-10)(n-8)(n+3)(n+4)(n+6)],
P(n) = 3659n^4 + 309616n^3 + 815324n^2 + 30852016n - 93987840.
```

It is nonzero for every `n>=32`: after `n=32+z`, the relevant numerator
and denominator have strictly positive coefficients. The rate Jacobian
is `diag(2/a_i) M`, with determinant `16 det(M)/prod(a_i)`. Therefore these
four real normalization directions exist uniquely at first order, for
every fixed complete finite input vector. This is not a numerical
evaluation of `F`, `D^rate` or the four physical coefficients.

For eventual numerical matching, the diagnostic also supplies
`rate4_intervals`: it accepts four complete rational hard enclosures and
four conversion enclosures. Transcendental values are enclosed, not
rounded and declared exact. It returns the coordinate intervals and,
separately, the direct stable matching projection. Recombining the
individual coordinate intervals would discard their correlations and
can give a much looser bound. The conversion intervals may now use the
proved bound below; the complete hard intervals are still missing.

### First-order physical radiation conversion: an assembled bound

There is sufficient existing input to bound `D^rate` for the retained
classical tree, without assigning any hard finite coefficient. This
requires combining the earlier sectors, not reinterpreting the selected
S296 matter-only result as already complete.

The [S293](../problems/P8/s6/continuation/s6_293/FORMULATION.md)/S294 matter
virtual pole is `A_m Bsoft/(8pi^2 kappa epsilon)`. The full-D
[S302 metric pole](../problems/P8/s6/continuation/s6_302/notes/assembly.md)
is `-2 T_D Bsoft/(16pi^2 kappa^2 epsilon)`, hence
`A_G,D Bsoft/(8pi^2 kappa epsilon)`. Their sum has the **full dimensional
Born** `A_D` as coefficient. On the retained physical sheet
`K0=-2 Re Bsoft`.

Keep that same `A_D` in both the virtual interference and the real soft
subtraction. After division by the fixed four-dimensional Born rate, the
paired numerator, with the common `1/(8pi^2 kappa epsilon)` stripped, is

```text
(A_D/a)^2 * [exp(2 epsilon log x_det)*p(epsilon)*K_epsilon
            - exp(2 epsilon log x_ref)*K0],
x_ref = 1,  p(0)=1,  p'(0)=EulerGamma-2-log pi.
```

It vanishes at `epsilon=0`. Its derivative there is
`K1+p'(0)K0+2K0 log(x_det/x_ref)`; the derivative of `A_D` cancels.
This explicitly retains the full-D Born term whose omission would change
the S302 finite hard reference. The imaginary Coulomb pole is not
discarded from the amplitude; only its real-rate contribution is paired.

[S335](../problems/P8/s6/continuation/s6_335/notes/continuity.md) supplies
the fixed-resolution dimensional limit of the entire retained 47-graph
real-minus-soft remainder, with the same dimensional Born and reference
two-body phase. Therefore at first order

```text
D_i^rate = [K1+(EulerGamma-2-log pi)K0]/(8pi^2 kappa)
         + K0 log(1/256)/(4pi^2 kappa) + R_47(q_i;1/256).
```

The term `2 Re L_m/a` separately retains the gravity-Born/matter-loop
interference through the full denominator. It is not counted again as
part of `D^rate`.

All four points have `delta=1`, and `1/256<1/192`, so the sharper
[S335 low-window bound](../problems/P8/s6/continuation/s6_335/notes/bounds.md)
applies. Using the [S296 angular bounds](../problems/P8/s6/continuation/s6_296/notes/continuity.md),
`|K0|<=800`, `|Delta_soft|<112/kappa`, `pi^2>9` and `log 2<7/10`,
an explicit bound is

```text
|D_i^rate| < {112 + 800*(8*7/10)/36
 + (361/324)*[128 B x + B^2 x^2/2 + 12288 x]/36}/kappa
 < 10^13/kappa = 10^-787,
B=4*10^9, x=1/256.
```

The larger fixed-Born factor `361/324` is retained conservatively. Thus
the induced contribution to the four-direction matching projection obeys

```text
|sum_i w_i a_i D_i^rate/2|
 < 15*(650lambda)/2 * 10^-787 < 10^-783*lambda.
```

This resolves the **first-order conversion budget at these four points**
for the retained classical Born tree. It does not bound `F`, a one-loop
square, higher-order hard/radiative insertions such as the independent
S336 curvature term, or an added classical operator. No forward or
detector-zero limit is taken, and no all-order inclusive probability or
quantum decoupling theorem follows.

## 4. Explicit finite covariant representative

Use `Y=(nabla Phi)^2`, `X=Y/kappa0` and the original switch

```text
V(X) = (1-X)^1024/[X^1024+(1-X)^1024].
```

One explicit representative of the four finite directions is the
following addition to the Lagrangian scalar (with volume density restored):

```text
Delta L = hbar V(X) * [
   alpha Y^2/(32 pi^2 kappa^2)
 + (beta+4alpha) mu^2 Phi^4/(384 pi^2 kappa^2)
 + c_RH R_old H_heavy + (delta_kappa/2) R_old
].
```

`H_heavy` denotes the field, not the shape `H(q;n)`. In the selected
representative, at `mu=1`,

```text
A = alpha/(16 pi^2 kappa^2),     B = 2g c_RH/kappa,
C = beta/(16 pi^2 kappa^2),      D_N = delta_kappa/kappa^2.
```

The literal `Y^2` vertex is `2(S-4mu^2)` and the `Phi^4` vertex is `24`,
so the mass-dependent constant cancels correctly and the contact is
`[alpha S+beta mu^2]/(16pi^2 kappa^2)`. The `R H_heavy` exchange normalization
is the retained [S294 one](../problems/P8/s6/continuation/s6_294/notes/matching.md).
The Einstein term uses [S285's sign](../problems/P8/s6/continuation/s6_285/notes/volume.md):
`R_old=-R_P8`, with original action `+kappa R_old/2`. Expanding
`-T/(kappa+delta_kappa)` gives `+delta_kappa T/kappa^2`, as required.

This is an explicit new finite off-vacuum lift, not a reconstruction of
every curved counterterm from flat data. It introduces no light quadratic
term or flat heavy onepoint at this order. The switch satisfies
`V=1+O(X^1024)` at the vacuum and has a zero of order 1024 at the reference
clock. Its action jets through order 1023 on that clock are unchanged for
every finite choice of these coefficients. This exact statement does
**not** bound nearby clock-tube errors for arbitrarily large coefficients,
establish a healthy full higher-derivative theory, fix other curved
counterterms or prove a same-parent quantum bounce.

## 5. The actual matched object and its reference independence

For a real hard function define the four-sample subtraction

```text
(P f)(q) = f(q) - b(q) M^(-1) [f(q_1),...,f(q_4)]^T.
```

It obeys `P^2=P`, `P(b xi)=0` for every constant four-vector `xi`, and
`(P f)(q_i)=0`. The matched real non-matter hard correction is

```text
Re F_RATE4(q) = (P Re F)(q)
             - b(q) M^(-1) [a_i D_i^rate/2]_i.
```

The imaginary part of `F` is unchanged by real finite counterterms.
Matching four rates does not cancel absorptive cuts.

Two separate invariance statements matter:

- Under a finite hard-reference change `F -> F+b xi`, the solved
  coordinates change as `theta -> theta-xi`. The total hard amplitude
  is unchanged, including away from the four centers.
- Under a redistribution between the hard and soft references
  `Re F(q_i) -> Re F(q_i)+a_i z_i/2` and
  `D_i^rate -> D_i^rate-z_i`, `d_i` is unchanged. Changing only one side
  would be inconsistent. The physical rate, not a convention-dependent
  hard amplitude in isolation, is the invariant observable.

Thus the new calculation is a physical normalization construction, not
the instruction to set four minimal-subtraction coefficients to zero.

The next required bound is on the **subtracted remainder**, not on those
four arbitrary reference coordinates again. If a separately justified
linear analytic functional `L` has

```text
L b = [2, 2(n+2)/(n-2)^3, 0, 0],      w^T=(L b)M^(-1),
```

then its matched contribution is

```text
L F - sum_i w_i Re F(q_i) - sum_i w_i a_i D_i^rate/2,
```

with the appropriate real part understood. The four finite directions
cancel, and `sum |w_i|<15` controls the sample-error contribution. This is
an algebraic reduction **conditional on an admissible L**. An unqualified
forward derivative of a massless-gravity hard amplitude is not such a
functional: its pole, cuts, infrared convention and limiting operation
still require the V/G analysis. No dispersion relation is imposed on an
inclusive rate.

## 6. What four conditions do not determine

The curved coupling cannot generally be read individually from `B`. In
fact, with `Q(q)=sum 1/(n-channel)`, there is an exact identity

```text
H(q;n) = (n+2) Q(q) - 3.
```

A finite cubic/heavy-residue insertion `r Q(q)` is therefore in the same
on-shell span. Changing it by `delta r` can be compensated everywhere in
this four-scalar amplitude by

```text
delta B = -delta r/(n+2),       delta C = -3 delta r/(n+2).
```

For fixed complementary data the displayed covariant representative
does determine its `c_RH`. If those data are also varied, the rates fix
only the combined heavy-exchange direction. This is an exact degeneracy,
in addition to the finite-precision conditioning issue in the preceding
comparison. It is not an assertion of full covariant equivalence:
curved/source amplitudes can distinguish the lifts.

A heavy mass insertion proportional to `sum 1/(n-channel)^2` and higher
local EFT terms need separate accounting. The four conditions do not
force their finite values or the off-point subtraction remainder to
vanish. For a concrete algebraic check, the crossing-symmetric local
polynomial

```text
delta A(q) = zeta * product_{v in {72,118,176,328}} (S(q)-v)
```

vanishes at all four centers, but its forward `v_cross^2` coefficient at
`s=2+v_cross,t=0,u=2-v_cross` is `-25579520*zeta`. This is an allowed
polynomial ambiguity at the amplitude level, not a UV-realizability claim
or an exclusion of RATE4. A justified truncation or bound on such
complementary terms is essential to a positivity verdict.

There is also a limited useful decoupling statement. Hold `n`, `g`,
`kappa0` and the original independent functions fixed. If
the **complete** non-matter hard term and physical rate conversion tend
to zero on a specified compact nonforward domain containing all four
centers as `kappa -> infinity`,
then the same is true after RATE4 subtraction. `M^(-1)` is fixed, so it
preserves that limit, including uniformly on a compact domain where the
four shapes are bounded. This is a linear-algebra implication, not proof
of those quantum/infrared hypotheses, a limit uniform in `n -> infinity`,
or an omitted-loop estimate. Large raw finite coefficients must not be
confused with a bound on the matched amplitude or curved response.

## 7. Closure-driven next calculation: RATE4.REMAINDER

This construction completes the **four-direction first-order normalization
subtask**, not physical MATCHED status. No further user-selected numerical
coefficient is needed for that subtask. The first-order radiation-conversion
contribution is now bounded well below the matching-error target.

The remaining hard calculation has a finite **original one-loop graph
inventory**, distinct from arbitrary additional EFT matching data. The
[S297 graph count](../problems/P8/s6/continuation/s6_297/notes/source.md)
gives `sum(d_v-2)=4`; its five possibilities are four cubics, two cubics
and one quartic, two quartics, a cubic and a five-field vertex, or a
six-field vertex. The actual nonminimal higher jets either have too many
legs or cannot close the required species assignment when a metric is
inserted. This leaves the following assembly for the known part of `F`:

| Original hard sector | Retained source | No-double-counting rule |
|---|---|---|
| `C_tree/kappa` minimal dressing | S293 | Include its whole endpoint/OS terms once |
| `g^2/kappa` minimal heavy dressing | S294 | Its S290 bubble/endpoint is already included |
| `kappa^-2` light/graviton/ghost sector | S302 | Retain the full-D pole and old Gram finite terms |
| `kappa^-2` H/Proca/M1 spectator insertions | [S305](../problems/P8/s6/continuation/s6_305/FORMULATION.md) | Add only its spectator increment, not its augmented S302 reference a second time; no duplicate light-Phi insertion |

The no-internal-graviton contribution is already in `L_m`, not in this
table. The S305 fixed local gapped curvature piece is
`-(log n+2)*(S+12)/(640pi^2 kappa^2)` together with its finite Newton
insertion. Both are in the removed `S,1,T` span, as is its M1 scale
variation. Their nonlocal spectator remainder is not removed. This
inventory is not a bound on the assembled expression or a declaration
that independent higher-EFT matching vanishes.

The next admitted calculation is:

1. Assemble the complete `F` in the same reference, using the `D^rate`
   formula and bound above, and separating
   calculated terms, directions removed by `P`, and genuinely complementary
   finite data. Do not substitute the S302 minimal metric term or S296
   selected soft conversion alone for the complete inputs. Account for the
   distinct heavy-residue and heavy-mass issues above.
2. Establish or bound the particular subtracted projection required by
   the selected V test. Keep the earlier absolute `lambda` matching-error
   target, relative to the reference `b20=4lambda` and the positive S239
   first-loop shift below `4lambda*10^-203`. Errors in `L F`, the four
   samples, rate conversion and omitted orders need separate entries.
   The existing `15` norm controls only the four sample errors.
3. Identify the corresponding curved/source projection of the same
   finite lift and remaining data. Exact clock jets do not replace a
   finite-domain coefficient/state/error bound for B.

Success means a complete, physically justified bound in that scope, not
another small selected loop. An explicit inconsistency would exclude only
the stipulated candidate/ansatz. An unbounded complementary projection or
an interval straddling the decision margin leaves the verdict conditional
or inconclusive; it is not assigned zero. No new all-orders UV-construction
requirement is added. M/V/G/B/R and original P8 remain OPEN.

## 8. Reproduction and evidence level

The [read-only candidate diagnostic](../scripts/p8_rate4_candidate.py)
pins the two preceding diagnostics and 16 parent reports/source manifests,
checking 307 protected files before and after execution.

```sh
.venv/bin/python -B scripts/p8_rate4_candidate.py
.venv/bin/ruff check scripts/p8_rate4_candidate.py
.venv/bin/ruff format --check scripts/p8_rate4_candidate.py
```

The symbolic determinant and positive-half-line argument are checked,
along with the literal covariant contact vertex, Newton sign, heavy
residue identity, reference-switch factors and complementary-polynomial
coefficient. Independent `Fraction` and SymPy solves agree on 16 matching
coordinates at four mass ratios, including the exact original ratio.
Sixteen synthetic rate equalities and twelve off-center reference-shift
checks pass. Hard/soft redistribution, projector identities, untouched
imaginary parts, the conditional decoupling algebra and nine malformed or
incomplete exact inputs are checked. Sixty-four interval-box corners
check both coordinate enclosures and the direct correlated projection;
three malformed or incomplete interval packets are rejected.
Synthetic amplitudes and rates are not physical matching evidence.
The target-size estimate is inherited from
the pinned S239 bound, with the positive full-Born comparison checked here.
The additional same-D Born/pole cancellation, resolution logarithm,
rational radiation/projection bounds, one-loop valence inventory and
local spectator-axis cancellation are checked explicitly. The physical
integrability and graph-species arguments use the cited retained proofs;
the diagnostic does not rerun all of their frozen symbolic calculations.

The claim level is VERIFIED_N for this scoped algebra and written
construction, not CERTIFIED physical matching or FORMALIZED quantum field
theory. Frozen scientific sources and certificates are not edited.
