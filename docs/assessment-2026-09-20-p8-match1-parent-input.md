# MATCH-1 parent-input comparison: physical rates and curvature information

Date: 2026-09-20. This continues the
[finite-normalization audit](assessment-2026-09-20-p8-match1-renormalization-obstruction.md)
and the [closure-driven plan](p8-closure-plan.md). Original MATCH-1 and
M/V/G/B/P8 remain OPEN. No finite coefficient is assigned to the old parent.

## Decision and research admission

The most direct additional input is a set of **physical, finite-resolution
scattering normalization conditions**, supplied by a specified parent or
explicitly adopted as conditions of a separately named quantum-EFT family.
Use the name **QG2-H8A420-RATE4** for that parameterized family. At the time
of this comparison it named the proposed input route, not a verified
member: no four physical rate values or complete error bounds had been
obtained or stipulated.

**Subsequent decision:** the user approved developing the separately
normalized candidate. The [RATE4 v1 construction](assessment-2026-09-20-p8-rate4-candidate.md)
now stipulates its matter-calibrated rate targets and solves the four
first-order normalization directions. The comparison below records the
prior information boundary; its intervals remain conditional, and no
complete physical error bound or original-parent MATCHED status is added.

The calculation below addresses gate M, not the smallness of another known
loop. Its threshold is matching uncertainty below `lambda`, one quarter
of the original tree coefficient `4*lambda`. It yields two new results:

- A conditional, exact interval map from four finite-resolution rates to
  the MATCH-1 coefficient combination, with an explicit sufficient error
  allocation below that threshold.
- An information limit: the same low-energy data, at the previously
  proposed precision, do not tightly determine the individual `R H`
  coupling. The stable vacuum projection must not be mistaken for a
  controlled reconstruction of the curved parent.

The former is a usable matching interface, not a prediction of its inputs.
The latter rules out *using this four-value precision alone* as a
complete curved-matching argument; it does not exclude the model or a
better-constrained parent. A successful physical input packet advances M.
A packet with insufficient bounds is inconclusive. Without a packet, the
old retained-data route stays CONDITIONAL ONLY.

## 1. Comparison of concrete input routes

| Route examined | What it could supply | Missing bridge and decision |
|---|---|---|
| Quantize the retained S233/S238 heavy-scalar model more accurately | Additional calculable matter contributions | Its Einstein-gravity EFT still admits the finite direction established in the preceding audit. It cannot determine that direction from the old normalization conditions. Do not pursue as the matching-resolution step. |
| Specify four physical scattering normalizations | Decision-relevant finite information in a common infrared prescription | This defines a new, explicitly conditioned quantum-EFT family, or accepts data from an independently specified parent. It needs values, full errors and a separate curved matching argument. Recommended conditional route. |
| Import published string-theory gravity coefficients | Genuine microscopic matching in the theory actually studied | No identification with this light/heavy/Proca spectrum, full DHOST functions, physical matter frame or state has been supplied. Not a plug-in input for QG2-H8A420. |
| Set minimal-subtraction finite parts to zero or impose an order-one prior | A numerical representative | Adds an unsupported restriction on the physical family. Not accepted as matching evidence. |

The first limitation follows from the explicit surviving counterterm, not
from a claim that matter calculations are useless. Renormalization-group
running also needs boundary values: it transports a finite coefficient
but does not select it. The distinction between measurable or
microscopically matched local coefficients and calculable nonlocal effects
is explained in [Donoghue and Holstein, section II](https://arxiv.org/html/1506.00946).

A concrete microscopic comparison was checked:
[Calmet, Kiritsis, Kuipers and Lust](https://arxiv.org/html/2408.15146v1)
match curvature-squared terms in particular string compactifications.
Their coefficients depend on compactification data, and their displayed
matter sector does not supply this project's scalar/heavy/Proca matching.
**Our inference:** using those numbers here would require a new common
parent construction, including all induced matter operators; their
gravitational coefficients alone cannot determine the four required
residuals. No string embedding or naturalness assumption is adopted.

Physical normalization is therefore the narrower next route than importing
an unrelated microscopic example. It is important to distinguish two uses:
measuring or deriving the normalizations gives information about a parent;
declaring them defines a new candidate. Neither use follows from choosing
a subtraction convention. This distinction remains explicit in RATE4.

## 2. Specify physical inputs, not four bare amplitude numbers

Keep the original mass-one units, physical matter frame and the four
kinematic points

```text
(s,t,u) = (8,-2,-2), (10,-3,-3), (12,-4,-4), (16,-6,-6).
```

For this input route fix detector resolution `x_det=1/256` in the retained
mass-one convention. At each hard center use the inclusive differential
rate, summing unresolved radiation with the declared total-energy rule.
Remove the infrared regulator at fixed nonzero resolution and fixed
nonforward angle first. A finite experimental bin would require its
averaging error in the input packet; none is silently identified with a
point value. These are theoretical normalization observables, not existing
measurements of the hypothetical model.

The choice stays inside the retained compact domain. Its transfer margin
is `delta=1`, and `x_det<delta/192`, so the sharper existing real-emission
bound is available. This detector resolution is **not** the S302 analytic
reference `resolution/nu=1`. A finite conversion to that reference must
be explicitly subtracted or included in the error.

Inclusive rates are appropriate normalization inputs, but are not
crossing-analytic amplitudes on which to impose a dispersion relation.
This distinction is also important in the infrared literature; see
[Prabhu, Satishchandran and Wald, introduction](https://arxiv.org/html/2203.14334).
We use their distinction between inclusive observables and a fundamental
scattering construction, not an assertion that they establish this
parent's scattering theory or its G hypotheses.

### Positive Born normalization

Let `a_i` be the full, unexpanded original Born amplitude at point `i`:

```text
a_i = C + g^2*sum_channels 1/(n-channel) - T(q_i)/kappa,
C = -g^2*[3/(n-2)-2/(n-2)^2].
```

Here `T` is the Newton shape defined in the preceding audit, and the
original `n`, `g`, `kappa`, `lambda` are unchanged. Exact rational
arithmetic, independently evaluated in SymPy, gives

| Point | Strict interval for `a_i/lambda` |
|---|---|
| `(8,-2,-2)` | `(136,137)` |
| `(10,-3,-3)` | `(228,229)` |
| `(12,-4,-4)` | `(344,345)` |
| `(16,-6,-6)` | `(648,649)` |

In particular `0<a_i<650*lambda`. Normalize each inclusive rate by its
Born rate with the same flux, identical-particle and phase-space factors.
Do not normalize by a truncated gravity expansion or by the matter Born
amplitude alone.

## 3. Exact conditional rate-to-matching enclosure

Write `R_i^c` for this normalized rate after subtracting the selected
known finite soft/radiation conversion. Let `K_i` be the real known hard
amplitude in the same reference. Here `r_i` denotes the four-shape part
alone, and `d_i` the remaining real amplitude structures. The required
complete relation is

```text
R_i^c = 1 + 2*(K_i+r_i+d_i)/a_i + Q_i.
```

This equation defines which remainder must be controlled. It is not a
claim that the retained first-loop calculations have bounded `Q_i`.
`Q_i` includes every remaining rate contribution at the desired accuracy:
higher-order interference and squares, unaccounted radiation, reference
conversion and any relevant regulator/bin errors. Contributions moved
between `K_i` and the conversion must not be counted twice.

Suppose a supplied packet encloses `R_i^c` in `[R_i^-,R_i^+]`, encloses
`K_i` in `[K_i^-,K_i^+]`, bounds `|Q_i|<=q_i`, and bounds the omitted
amplitude structures by `|d_i|<=e_i`. Then the four-shape residual
has the rigorous interval

```text
r_i^- = a_i/2*(R_i^- - 1 - q_i) - K_i^+ - e_i,
r_i^+ = a_i/2*(R_i^+ - 1 + q_i) - K_i^- + e_i.
```

The previous exact weights `w_i` project these intervals to `C_pos`:
multiply each interval by `w_i`, choosing endpoints according to its
sign, and sum. At the original `mu=1` no extra mass factor occurs;
at general `mu`, divide by `mu^2`. Constant and Newton nuisance anchors
remain eliminated, not assigned values. Shared uncertainties may be
treated more sharply if their correlations are proved; independent
intervals give the conservative result used here.

### A sufficient, explicitly conditional budget

The following is one useful target allocation, not a result about the
actual parent's uncomputed errors:

1. The combined normalized-rate interval half-width and all remaining
   conversion/radiation/higher-interference errors, **excluding only the
   separately bounded hard square below**, are at most `1/10000`.
2. The entire complex hard correction whose square is omitted obeys
   `|delta A_hard|/a_i <= 1/200`. Its square then contributes at most
   `1/40000` to the normalized rate. This bound cannot be inferred from
   a small real interference or small known terms alone.
3. Known-hard subtraction uncertainty plus all omitted four-shape
   amplitude structures have total absolute error at most `lambda/100`.

These hypotheses imply

```text
half-width(r_i) < lambda * [650/2*(1/10000+1/40000)+1/100],

half-width(C_pos) < (243/320)*lambda < lambda.
```

Thus the apparent absolute `10^-600` scale does not require a
`10^-600` *relative* rate accuracy. The normalized-rate budget here is
`0.01%`, together with the separately stated hard and remainder bounds.
This is a theoretical precision allocation, not an experimental
feasibility assessment or evidence that any of its hypotheses hold.
No central rate values have been chosen, so no central value or sign
for `C_pos` is reported.

### What existing radiation calculations contribute

[S335's actual low-resolution bound](../problems/P8/s6/continuation/s6_335/notes/bounds.md)
gives a normalized remainder below `10^-787` for the **selected complete
one-real tree minus its leading soft term** at these kinematics and
resolution. That is ample for this particular piece of the proposed
budget. It is not the total `Q_i` bound. In particular, it does not supply
the finite virtual matching, uncomputed hard sectors, all-order
interference or all-multiplicity radiation control.

[S296](../problems/P8/s6/continuation/s6_296/notes/inclusive.md) and
[S297](../problems/P8/s6/continuation/s6_297/FORMULATION.md) also explicitly
distinguish a consistently truncated rate from an unqualified hard
amplitude square. The square in the above sufficient budget is an
**additional bounded error hypothesis**, not a loop square imported into
their certified first-loop result. Tightening the tiny S335 bound again
would not supply the missing hypotheses.

## 4. Why four accurate rates do not reconstruct the curved parent

Recall the four-shape matrix `M(x)` with `x=n/mu` and coefficient vector
`(A,B,C,D_N)`, where `B=2*g*c_RH/kappa`. Its stable projection determines
`2*A+2*(x+2)*B/(x-2)^3`, not each coordinate with comparable precision.

This identification of `B` with `c_RH` holds with the other heavy-residue
data fixed. The subsequent [RATE4 construction](assessment-2026-09-20-p8-rate4-candidate.md)
also exhibits an exact degeneracy when those data vary. The conditioning
bound below remains valid in the restricted identification; four elastic
rates alone do not resolve that additional curved-parent ambiguity.

The row extracting `B` alone is exactly

```text
v = F(x)*(-7681/24, 4595/6, -12457/24, 293/4),

F(x) = (x-16)*(x-12)*(x-10)*(x-8)*(x+3)*(x+4)*(x+6)/P(x),

P(x) = 3659*x^4 + 309616*x^3 + 815324*x^2
       + 30852016*x - 93987840.
```

For `x>=32`, `F(x)>0`, and

```text
sum_i |v_i| = (10069/6)*F(x),
lim_(x->infinity) sum_i |v_i|/x^3 = 10069/21954.
```

If each amplitude residual has an otherwise unrestricted independent
error interval of radius `epsilon`, the exact possible radius in `c_RH`
is

```text
radius(c_RH) = kappa/(2*g) * epsilon * sum_i |v_i|.
```

This is not just an upper-bound artifact of a numerical solver. The
error-box corner with signs matching `v_i` attains it; invertibility of
`M` gives the corresponding algebraic coefficient variation. Additional
parent constraints or proven correlations could shrink this set, but
they would be additional information. These algebraic variations are
not asserted to be realizable UV theories.

At the original parameters and the preceding target
`epsilon=lambda/15`, exact arithmetic gives

```text
9*10^793 < radius(c_RH) < 10^794.
```

For a scale comparison only, reducing that radius to `sqrt(kappa)`
would require `epsilon` below an exact threshold between `10^-996`
and `10^-995` (approximately `7.144589*10^-996`). The unit
`c_RH/sqrt(kappa)` is natural for the canonically normalized metric-heavy
mixing; bounding it by one is **not** asserted to establish a healthy
curved theory or adopted as a new physical condition.

The mechanism is transparent at large `x`:

```text
H(q;x) = 10/x + [S(q)+8]/x^2 + O(x^-3).
```

The constant column absorbs the first term, and the local scalar column
absorbs the leading nonconstant term. Only the small remaining shape
separates `B`. Their particular MATCH-1 combination is well conditioned
because its leading cancellation is the same; the individual curved
coupling is not. Changing numerical precision cannot remove this
information limit of data with fixed error intervals.

Consequently, do not demand a needlessly precise fit of individual
coefficients just to obtain `C_pos`. Conversely, do not export a precise
`C_pos` as a precise curved counterfunctional. B needs separately
matched curvature/source information or a direct bound on its actual
projection. Existing radiative directions such as S336's independent
finite term remain separate as well.

## 5. Input status and continuation boundary

| Required input for an actual RATE4 application | Current status |
|---|---|
| Parent identity, physical frame, mass/coupling definitions and finite-resolution observable | Original reference and a proposed observable specified; no completed new quantum parent |
| Four physical inclusive-rate intervals | Not supplied by the retained data or the examined microscopic example |
| Complete conversion to the common hard reference | Selected components available; full bound not established |
| Complete rate and four-shape remainder bounds | Sufficient thresholds derived; actual full bounds not established |
| Independent curved/source matching for B | Not supplied by the four-value projection; separate information requirement quantified |

This comparison therefore selects a **conditional physical-normalization
route**, not a new numerical representative and not MATCHED status for
the original parent. A proposed RATE4 member must explicitly supply its
normalizations and evidence for the bounds above. Declaring them would
be new model assumptions; deriving them from a microscopic candidate
would require that candidate's spectrum, frame, matching and errors.

The next admissible calculation is to establish one such complete input
packet or a separately justified direct bound on the needed projection.
Do not substitute more precise known-loop terms, four synthetic rates,
or imported curvature coefficients for that packet. If work proceeds
only by stipulating the rate conditions and remainder hypotheses, its
deliverable must remain a conditional classification, not closure of
the unchanged MATCH-1 route or original P8.

The user subsequently chose the separately normalized candidate. Its new
conditions and their first-order solution are recorded in the
[RATE4 v1 specification](assessment-2026-09-20-p8-rate4-candidate.md).
This explicit choice, rather than the recommendation alone, authorizes
the new model definition. The original parent remains unchanged and
conditional. The next obligation is the complete subtracted remainder,
not another request for arbitrary finite coefficients.

## Reproduction and evidence level

The [read-only input diagnostic](../scripts/p8_match1_rate_input.py)
provides exact interval-conversion and projection functions with no
default-zero errors. Missing complete-error arguments, floating-point
inputs, reversed intervals and invalid domains are rejected. It does not
authenticate physical evidence merely because numbers have been entered.

```sh
.venv/bin/python scripts/p8_match1_rate_input.py
.venv/bin/ruff check scripts/p8_match1_rate_input.py
.venv/bin/ruff format --check scripts/p8_match1_rate_input.py
```

It checks 268 protected inputs, including the preceding diagnostic and
14 parent reports with their source manifests, before and after the run.
Original SymPy and independent `Fraction` arithmetic agree on all four
Born values and 16 additional curvature-extraction weights. The
half-line sign, exact asymptotic conditioning, error thresholds and
original-parameter intervals are checked. Forty-eight synthetic interval
corners verify exact enclosure endpoints, and ten invalid inputs are
rejected. Synthetic fixtures have nonzero nuisance, imaginary, shape and
rate-remainder terms, with nonzero input intervals; they are not physical
matching data. The general light-mass normalization is checked as well.

The ledger level is VERIFIED_N for the scoped algebra and interval
argument, not CERTIFIED physical matching or a formalized QFT proof.
No frozen scientific source or certificate is modified.
