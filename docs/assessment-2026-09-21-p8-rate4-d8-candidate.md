# RATE4-D8 v1: conditional candidate and parent-bound witness contract

Date: 2026-09-21 (PDT). The user approved developing the separately named
RATE4-D8 candidate recommended by the
[RATE4.PARENT comparison](assessment-2026-09-21-p8-rate4-parent.md).
That approval authorizes this conditional specification. It does **not**
establish any of its new physical bounds or the existence of a common
quantum parent satisfying them. Original RATE4 v1 and the frozen P8
calculations are unchanged. Physical M/V/G/B/R/P8 remain **OPEN**.

**Result:** the D8 candidate now has a constructive first-order matching
map, complete-input interval interfaces, an order-by-order formal
normalization lemma, and a sufficient analytic witness for its entire
omitted local tail. The known rate-blind contacts are quantitatively
controlled **if** that tail bound is established. No physical matching
values, actual tail/pole/source bounds or higher-order error bounds are
filled in by this construction.

## 1. Conditional identity and unchanged normalization targets

Name: **QG2-H8A420-RATE4-D8, v1 (parent-bounded conditional candidate)**.
Retain the classical action, spectrum, original physical matter frame,
light mass/residue and heavy onepoint reference of
[RATE4 v1](assessment-2026-09-20-p8-rate4-candidate.md), with

```text
mu=1, kappa0=kappa=10^800, lambda=10^-600,
n=10^200/512+2, g=1/8192,
C_tree=-g^2[3/(n-2)-2/(n-2)^2].
```

There is **no fifth target**. At

```text
q1=(8,-2,-2), q2=(10,-3,-3), q3=(12,-4,-4), q4=(16,-6,-6)
```

retain exactly `R_i^target=1+hbar*2 Re L_m(q_i)/a_i`, with higher target
coefficients zero as in RATE4 v1. `L_m` is the complete S239 reference,
not a newly chosen scalar coefficient. `a_i` is the full unexpanded Born
amplitude, positive and below `650lambda`. The detector resolution is
`1/256`, with the same nonforward, fixed-resolution regulator-removal
order. These are stipulated physical normalization conditions, not
measured rates or a prediction of the unchanged original parent.

The candidate is conditional on a **single** quantum prescription supplying
the decomposition and bounds below in one reference and controlled domain.
It is not a demonstrated nonempty family of physical theories. Neither
minimal-subtraction zeros nor naturalness estimates constitute that proof.
The local/nonlocal EFT distinction is reviewed by
[Donoghue and Holstein, section II](https://arxiv.org/html/1506.00946);
it does not determine the numerical matching data needed here.

## 2. The bounded-EFT specification and its open obligations

In mass-one units let `S=s^2+t^2+u^2`, `U=stu` and separate the first-order
non-matter hard correction as

```text
F=K+p4+r Q+eta J+D_N T+E_tail,
p4=c0+cS S+cU U+cSS S^2,
Q=sum_a 1/(n-a), J=sum_a 1/(n-a)^2,
T=sum_cyclic(2-2a-bc)/a.
```

`K` contains the complete specified known non-matter reference:
S293, S294 (already including S290), S302 and the spectator-only S305
increment, each once. The S239 matter term stays separate. Known light
and graviton nonanalytic terms are not put into `E_tail` by a Taylor
expansion through their cuts. At higher orders, additional nonanalytic
terms belong in the full higher-order error or explicit calculation.
The finite unknown anchors parameterized by `p4,r,eta,D_N,E_tail` are
excluded from `K`; a component formula that displays such an anchor must
not count it again as a known loop contribution.

The new low-degree matching coordinates are `(c0,cS,cU,cSS)`, **not** the
old `(A,B,C,D_N)` basis `(S,H,1,T)`. The former heavy-residue combination
is now part of the bounded complementary input. Since `H=(n+2)Q-3`,
`r` combines equivalent flat residue directions; it does not identify
the individual curved coupling `c_RH`.

The declared sufficient first-order parent-bound hypotheses are

```text
|r|<=g^2/4, |eta|<=g^2 n/4, |D_N|<=1/(4kappa),
E_tail=sum_(2i+3j>=5) c_ij S^i U^j,
N64,D8=sum_(2i+3j>=5) |c_ij|64^(2i+3j) <= 2000lambda.
```

These are **obligations to establish**, not verified properties. They
restrict neither the signs nor the individual values to zero. In the
linear insertion convention `eta=-g^2 delta_n`; that identity is not a
resummed unstable-particle pole condition. The Newton bound concerns a
residual coordinate in the declared reference, not an already measured
Newton coupling. Reference changes must transform the data and their
bounds together. Scale 64 is an analytic norm radius, not a claimed cutoff.

Separate acceptance obligations cover the complete rate conversion,
known-hard evaluation, omitted orders, admissible V/G functional and
same-parent curved/source/state response. In particular, the following
are different statements:

| Statement | Status |
|---|---|
| This conditional candidate and its four targets are specified | DONE as an authorized definition |
| Four local normalization coordinates exist for complete first-order inputs | PROVED conditionally below |
| Parent-supplied pole/source and tail bounds satisfy the displayed hypotheses | NOT ESTABLISHED |
| Complete hard values and higher-order errors are available | NOT ESTABLISHED |
| A common physical quantum parent realizes the conditions | NOT ESTABLISHED |
| V/G gives a necessary-test verdict, or B is controlled | OPEN |

## 3. Constructive matching and complete-input intervals

Write `b4=(1,S,U,S^2)`, `M_ij=b4_j(q_i)` and let `D_i^rate` be the complete
first-order inclusive-rate conversion in the same hard reference. With
`F_rest=K+rQ+etaJ+D_N T+E_tail`, define

```text
d_i=Re F_rest(q_i)+a_i D_i^rate/2,
c=-M^(-1)d.
```

The four first-order rate coefficients are then exactly their retained
matter targets. The determinant and rate Jacobian are

```text
det M=10011648,
det[diag(2/a_i) M]=16*10011648/prod_i a_i > 0.
```

This is an existence and uniqueness statement for the **four real local
coordinates at fixed complete complementary inputs**, not their physical
evaluation. Real finite counterterms do not cancel absorptive parts.
Every unknown input remains required; there is no zero default.

Under `F_rest -> F_rest+b4 xi`, the fitted vector changes by `c -> c-xi`.
The full hard amplitude is unchanged at the four centers and off them.
Likewise the joint hard/soft redistribution
`Re F_rest(q_i) -> Re F_rest(q_i)+a_i z_i/2`,
`D_i^rate -> D_i^rate-z_i` leaves `d` unchanged. Changing only one side
does not preserve the physical prescription.

For `l4=(0,2,0,32)`, the direct coefficient weights remain

```text
w4=(-3579/13036,5084/9777,-3509/13036,232/9777),
sum|w4|=3544/3259.
```

The new [diagnostic](../scripts/p8_rate4_d8.py) implements the exact map
and interval versions. The latter return each coefficient's marginal
interval **and** the direct interval for `l4 c=-w4 d`. Reconstructing that
projection from the marginal coordinate boxes needlessly discards their
correlations. Good matrix conditioning bounds error amplification, not
the size or physical health of the inferred central coefficients.

One formal covariant representative uses the labelled-slot construction
of the remainder audit: apply the polynomial `p4(D12,D13,D14)` to four
scalar slots, divide by `4!`, covariantize the contracted derivatives in
the original metric, and multiply by `hbar V(X)`. It has at most eight
derivatives and no new quadratic or onepoint term. The switch retains
the reference-clock jets through order 1023. This is a representative
for the local four-point data, not a determination of the independent
curved kernel or a bound on a nearby clock tube. Those inputs stay separate.

### Formal continuation of the defining rates

There is a useful algebraic statement beyond first order. Suppose the
complete finite rate coefficient at order k is calculable, including
all lower-order counterterm insertions and fixed complementary data.
Omit only the four still-adjustable real contacts `hbar^k b4 c^(k)` and
call the remaining rate vector `rho_known^(k)`. At that order these new
contacts enter only through their Born interference, so

```text
c^(k)=M^(-1) diag(a_i/2) [rho_target^(k)-rho_known^(k)].
```

The same invertible Jacobian works at each order. Insertions of these
new contacts into loops or radiation, or products with an earlier
non-Born amplitude, enter later orders. Thus, **if every required finite
rate coefficient exists and is supplied**, normalization can be solved
recursively as a formal series. No truncation of the rest of the EFT to
four operators is assumed; all required other counterterms remain in
the complementary input.

This does not establish convergence, unitarity, a finite all-order rate,
the complete quantum parent, or the `lambda/4` higher-functional error.
The diagnostic checks two four-order synthetic examples with nonzero
imaginary amplitudes and their squares retained. They check the recursion,
not a computation of the physical higher loops.

## 4. A concrete sufficient witness for the infinite tail

The norm hypothesis can be supplied directly. Here is a separate,
sufficient route that makes the required analytic evidence explicit.
It is **not an additional mandatory assumption** of the candidate.

Suppose the residual analytic matching function `A(S,U)` extends
holomorphically to a neighbourhood of the closed complex bidisk

```text
|S|<=R^2, |U|<=R^3, R>64,
sup |A(S,U)| <= M_R.
```

It may contain the four low-degree words: those are absorbed into the
fitted polynomial and do not change the high-degree coefficients.
The double Cauchy formula gives
`|c_ij|<=M_R R^(-2i-3j)`. With `t=64/R<1`, absolute summation therefore
proves

```text
N64,D8 <= M_R W(t),
W(t)=1/[(1-t^2)(1-t^3)]-1-t^2-t^3-t^4
    =t^5(1+2t-t^3-t^4)/[(1-t^2)(1-t^3)].
```

The four subtracted terms are precisely the weights below five. The
geometric series proves the infinite sum; finite monomial tests are
not substituted for this argument. On the bidisk the supremum must
bound the whole function, not just finitely many sampled values.

For the explicit optional witness `R=128`,

```text
W(1/2)=29/336,
M_128 <= 23000lambda  =>  N64,D8 < 2000lambda.
```

The exact maximal sufficient supremum from this estimate is
`(672000/29)lambda`; 23000 is a convenient conservative integer.
The interface `tail_norm_from_bidisk` evaluates the implication, but
does not authenticate holomorphy or the supremum. It rejects `R<=64`:
a coefficientwise boundary estimate at equal radius would not produce
this absolutely summable majorant.

There is currently **no supplied physical `A`, holomorphic-domain proof
or `M_128` bound** for RATE4-D8. In particular, known light/gravity cuts,
a large heavy mass, a regulator or a prior on Wilson coefficients do not
supply this witness. A different parent proof of `N64,D8` or a tighter
direct projection bound can be considered on its actual merits; changing
the defining bound would require recording a new candidate version.

## 5. What happens to the two known invisible contacts

The new hypothesis addresses the earlier obstruction rather than
denying it. For the previous

```text
R5=3259SU-18645S^2+2463528U-152808S+21316224,
```

the only tail word is `3259SU`; its other terms belong in `p4`. In the
one-direction family `zeta R5`, the tail norm is therefore
`3259|zeta|64^5`. Under the D8 cap,

```text
|delta b20|=902256|zeta|
 <=(7048875/13669236736)lambda < lambda/1000.
```

Similarly the same-angle blind polynomial is

```text
Delta=(s-t)^2(t-u)^2(u-s)^2
     =S^3/2-20S^2-36SU+256S-27U^2+320U-1024.
```

Its tail is `S^3/2-36SU-27U^2`, with norm
`1796|zeta|64^5` for `zeta Delta`. Its forward coefficient is `64zeta`,
so the same cap gives

```text
|delta b20| <=(125/1883242496)lambda < 10^-7 lambda.
```

These are one-direction checks, not a decomposition-dependent bound on
every separately labelled coefficient in a sum. The general parent norm
bounds the **net** tail and the complete projection instead. Without that
new physical input the original unbounded contacts remain valid EFT
comparisons. The old `+/-8lambda` R5 controls violate the proposed D8 norm;
they have not been excluded by the four rates alone.

## 6. Total-coefficient interface and unchanged error target

Let `L` be an independently admissible functional agreeing with ordinary
forward coefficient extraction on the analytic terms, and let
`B_ref=L(tree+L_m)`. For exact first-order input,

```text
b_matched = B_ref + L K - sum_i w4_i[Re K(q_i)+a_i D_i^rate/2]
            + L P4(rQ+etaJ+D_N T+E_tail),
P4 f=f-b4 M^(-1)[f(q_i)]_i.
```

The implementation `matched_coefficient_interval` encloses this expression
from supplied intervals, then adds all declared remaining sample-rate
errors and the independent higher-L error. It requires the tree-plus-matter
reference interval, complete `L K` and four `Re K(q_i)` intervals, complete
conversion intervals, and each residual bound explicitly. An interval
containing zero or a negative endpoint is not itself a physical V/G verdict.
No gravitational contour allowance or functional admissibility is inferred.

The previous D8 bounds apply:

```text
|L P4 E_tail| <= [332221287/6998649208832] N64,D8,
|L P4(rQ+etaJ+D_N T)| < 10^-200 lambda
```

for the declared original-parameter caps and the required pole convention.
The second is a first-order bound, not an all-order heavy-sector theorem.
With the complete input allocations in RATE4.PARENT, the non-reference
matching uncertainty stays below
`(1620429/3259000)lambda < lambda/2`. Any uncertainty in `B_ref` is added
separately; the interval function does not discard it.

The test suite contains a deliberately **synthetic** interval centered
on `4lambda` that lies in `(3.5lambda,4.5lambda)` under those allocations.
This tests interval plumbing only. It is not the physically evaluated
RATE4-D8 coefficient: its known-hard central values were not calculated.
The genuine next calculation must supply them and all required premises.

## 7. Closure-directed next milestone: RATE4-D8.WITNESS

The candidate-definition and conditional normalization work is complete.
The next milestone is **parent-bound evidence**, not another approval of
the same conditional specification or another negligible known-loop term.

1. Identify the actual residual analytic matching function in the common
   reference and prove the D8 tail bound, either directly or using a
   sufficient witness such as section 4. If the source cannot determine it,
   state precisely which parent input is absent. Four exact rates do not
   supply it. Do not enter a convenient supremum as if it had been measured.
2. Tie the residue, linear mass insertion and Newton/source coordinates
   to the same parent's physical prescription and prove their enclosures.
   The unstable heavy field is not an external stable scattering state.
3. Keep the common-parent curved/source/state obligations separate:
   individual `c_RH`, the S336 curvature direction, a controlled covariant
   lift, original state and canonical boundary transport are still needed
   where the B test uses them. Flat matching does not choose them zero.

After those evidence gates, assemble the complete known-hard projection
and error budget with the actual admissible V/G prescription. S293/S294
selected forward bounds and S302/S305 finite sample bounds do not yet
constitute that full functional. Formal normalization also does not prove
quantum decoupling with canonical inputs fixed or a global quantum bounce.
The 32-row coverage/integration gate remains distinct.

No new user-selected Wilson coefficients, fifth target, model restriction
beyond the approved D8 proposal or all-orders UV-completion requirement
has been introduced here. A failure to establish a hypothesis leaves this
candidate conditional; it is not a refutation of every allowed P8 parent.

## 8. Reproduction and scope

```sh
.venv/bin/python -B scripts/p8_rate4_d8.py
.venv/bin/ruff check scripts/p8_rate4_d8.py
.venv/bin/ruff format --check scripts/p8_rate4_d8.py
```

The read-only diagnostic pins 310 inherited inputs before and after its
run. It checks sixteen independent inverse entries, 256 complete-input
interval corners, 32 formal rate equalities through four orders with
eight nonzero absorptive-square controls, three genuinely infinite
analytic fixtures, both blind contacts, 256 complete-projection fixtures
and 28 invalid/incomplete-input rejections. Exact algebra checks accompany
the written recursive-order and double-Cauchy proofs; no finite sample
test is claimed to prove an infinite physical tail bound.

Evidence level is **VERIFIED_N** for the scoped conditional construction,
not CERTIFIED physical matching. The historical RATE4.PARENT comparison
and its receipt remain unchanged, including their then-unadopted status.
See the new [validation receipt](validation/p8-rate4-d8-2026-09-21.json).
The frozen 115,487-test publication suite is not rerun for these root
diagnostics. Original P8 remains open.
