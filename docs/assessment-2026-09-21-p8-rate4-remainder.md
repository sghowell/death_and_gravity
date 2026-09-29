# RATE4.REMAINDER: a precise matching obstruction and a bounded gapped sector

Date: 2026-09-21. This advances the user-approved
[RATE4 v1 candidate](assessment-2026-09-20-p8-rate4-candidate.md) under the
[closure-driven plan](p8-closure-plan.md). It does not change its targets,
the original parent, or any frozen scientific packet.

**Decision: CONDITIONAL ONLY.** The four physical normalization conditions
do not bound the complementary first-order projection to the required
absolute `lambda` accuracy. A local contact with at most ten derivatives
preserves all four first-order conditions without changing their four matching
coordinates, but changes the analytic `b20` coefficient by an arbitrary
multiple of `lambda`. This is an identifiability obstruction for the
stated EFT family, not a UV-complete counterexample or a model exclusion.

There is also a positive quantitative result: the complete **gapped
H/Proca spectator increment**, after the specified four-sample subtraction,
contributes less than `10^-1003*lambda` to its analytic `b20` projection.
The massless spectator, full metric functional, complementary finite data
and higher orders are not included in that bound. Physical M/V/G/B/R and
original P8 remain OPEN.

## Admission and decision threshold

This calculation addresses M, the complete matching remainder required
before V. The reference tree coefficient is `4lambda`, with
`lambda=10^-600`; the allotted absolute matching uncertainty is `lambda`.
The previous construction bounded the first-order radiation conversion
but did not determine complementary EFT data. The present task tests that
data sufficiency before attempting an expensive full numerical assembly.

Success would be a complete justified projection bound. A contradiction
would exclude only its specified ansatz. An unbounded complementary
projection triggers the plan's CONDITIONAL ONLY stop rule: further
known-sector refinement cannot repair missing physical normalization.
This audit reaches the third outcome. It does not label the full hard
calculation complete or claim that every possible additional parent input
has been exhausted.

## 1. Separate the complete expression before bounding it

In the original mass-one units let `b=[S,H,1,T]`, `M_ij=b_j(q_i)` and
`w=(2,2(n+2)/(n-2)^3,0,0) M^-1`, exactly as in RATE4 v1. Its proven
sample norm is `sum |w_i|<15` for `n>=32`.

The first-order hard term must have the form

```text
F = F293 + F294 + F302 + DeltaF305 + F_complement.
```

Here S293 includes its complete quartic endpoint/OS terms, S294 already
includes S290, and `DeltaF305` is only the additional H/Proca/M1 spectator
increment. The S239 no-internal-graviton matter reference is separate;
the S302 reference is not added twice through S305. S297 supplies the
original one-loop graph inventory. No independent finite term in
`F_complement` is silently assigned zero.

For an independently admissible full functional `L`, the required result is

```text
R[F] = L F - sum_i w_i Re F(q_i) - sum_i w_i a_i D_i^rate/2.
```

The `b` span cancels. An error budget must still contain the known hard
projection, the complementary projection, evaluation/subtraction errors,
the conversion error, and all relevant omitted orders. The diagnostic's
`conditional_remainder_interval` requires each of these explicitly. Its
interval arithmetic cannot authenticate their physical completeness.

The results below on analytic local or gapped terms use their ordinary
forward Taylor coefficient. They remain useful for any full prescription
agreeing with that coefficient on those terms. They do **not** construct
such a prescription for the full massless-gravity amplitude.

## 2. Heavy residue and mass are different problems

Set `Q=sum_a 1/(n-a)` and `J=sum_a 1/(n-a)^2`. The exact identity
`H=(n+2)Q-3` shows that `P Q=0` at every point, not merely at the samples.
A finite heavy residue change therefore disappears from this flat
subtracted projection. Its separation from the curved `R H_heavy`
coefficient remains unidentified; their covariant lifts need not agree.

A heavy mass displacement is different. Define

```text
K_n = 6/(n-2)^4 - sum_i w_i J(q_i),
K_3 = -sum_i w_i s_i t_i u_i.
```

The diagnostic derives both rational functions, proves their signs and
the following bounds on the entire half-line by polynomial coefficient
positivity after `n=32+z`, and independently cross-checks them with
Fraction arithmetic:

```text
0 < K_n < 10^6/n^6,
0 < -K_3 < 10^5/n,
lim n^6 K_n = 8591792/3659,
lim n K_3 = -8591792/10977.
```

Thus neither a mass insertion nor an independent six-derivative `stu`
contact is identically removed. The latter also changes the four fitted
coordinates when the four rates are kept fixed.

For example, expanding `g^2 Q(n+delta_n)` gives `-g^2 delta_n J` at first
order. **If** a parent supplies `|delta_n|<=n/4`, its first-order
subtracted effect is below `10^-388*lambda` at the original parameters.
This is a sufficient conditional estimate, not an imposed mass condition,
an expansion to all orders, or a heavy-resonance error bound. The actual
finite heavy mass/pole condition is still absent.

## 3. A lower-degree contact invisible to all four rate conditions

For identical scalar scattering with `s+t+u=4mu`, parity-even symmetric
local on-shell polynomials are generated by

```text
S=s^2+t^2+u^2,       U=stu.
```

This follows from the symmetric-polynomial generators: `e1=4mu`,
`e2=(16mu^2-S)/2`, `e3=U`. At fixed `mu=1` the complete contact basis
through eight derivatives is `(1,S,U,S^2)`. Its evaluation matrix at the
four RATE4 centers has determinant `10011648`, so it has no nonzero
contact-only kernel. Through ten derivatives the additional word is `SU`.
The five-dimensional space has the following one-dimensional kernel:

```text
R5 = 3259 S U - 18645 mu S^2 + 2463528 mu^2 U
                  - 152808 mu^3 S + 21316224 mu^5.
```

It vanishes at **each** `q_i=mu*(8,-2,-2), mu*(10,-3,-3),
mu*(12,-4,-4), mu*(16,-6,-6)`. Yet at forward crossing coordinates

```text
R5(2mu+v,0,2mu-v)
    = 18900480 mu^5 - 902256 mu^3 v^2 - 74580 mu v^4.
```

Adding `hbar*zeta*R5` therefore leaves all four first-order normalized
rates unchanged, with **no** compensating change of `A,B,C,D_N`, but
changes `b20` by `-902256*zeta` at `mu=1`. The heavy poles/residues,
Newton pole and first-order absorptive parts are unchanged. The physical
Born-tree radiation conversion is unchanged at this order: a finite
order-hbar vertex in that radiation correction costs another order.
Equality of the full all-order rates is not asserted; higher-order
matching and existence remain separate obligations.

### Explicit local covariant representative

This is more than an arbitrary list of four sample values. A constructive
local representative is obtained from four labelled scalar slots. In flat
space use the pair operators
`D_ij=-(partial_i+partial_j)^2`, which have Fourier eigenvalues
`(p_i+p_j)^2`. Substitute `(D12,D13,D14)` for `(s,t,u)` in `R5`, expand
into contracted scalar derivatives, and take the four slots coincident.
Multiply by `zeta/4!`. The 24 assignments in the four-field vertex each
give `R5` on shell; the factorial cancels exactly.

Replace those contracted derivatives by symmetrized covariant derivatives
in the original physical metric, restore the volume density, and multiply
by the already specified switch `V(X)` and `hbar`. This chooses one finite
local off-vacuum lift; it does not claim that all covariant lifts agree.
It is real and parity even, has four scalar fields and at most ten
derivatives. In ordinary dimensions `zeta` has mass dimension `-10`.

There is no quadratic or onepoint term. At the vacuum `V=1+O(X^1024)`,
so the relevant four-field vertex is unchanged. At the reference clock
the order-1024 zero of V preserves the action jets through order 1023
for any finite coefficient. This does not prove nearby-domain bounds,
nonlinear stability, a quantum state, or a UV realization of either lift.
The statement concerns perturbative local EFT comparisons; the added
higher-derivative term is not treated as an exact fundamental theory.

### Comparison with the actual decision margin

If this one-dimensional ambiguity alone receives the entire `lambda`
budget, the required bound is

```text
|zeta| <= lambda/902256.
```

No existing RATE4 condition supplies it. If other errors consume part of
the budget, the available radius is smaller. Correlated parent information
could instead bound the combined projection directly; an independent
bound on every coefficient is not mandatory.

For illustration only, `zeta=+/-8lambda/902256` gives opposite shifts of
size `8lambda`. They put the selected tree-plus-S239 reference on opposite
sides of zero; this is **not** a sign verdict for the missing full amplitude
or its gravitational inequality. On the retained compact physical window,
`|S|<=768` and `|U|<=4096` imply
`|R5|<=31478479488`, hence these comparison amplitudes have absolute
magnitude below `3*10^5*lambda`. This is not a relative perturbativity or
high-energy/cutoff bound. No finite value is adopted for the candidate.

## 4. More data require a stated truncation, not just a fifth rate

At the off-center physical point `(10,-2,-4)`, `R5=-37140096`; a rate
condition there would see this direction. But it would not determine the
entire complementary EFT.

Even restricting contacts to at most ten derivatives leaves seven
independent on-shell shapes `(1,S,U,S^2,SU,H,T)` at fixed n: the two
rational shapes have distinct heavy and massless poles and cannot equal
polynomials. Four conditions therefore leave three dimensions, not just
the displayed contact-only direction. Adding the independent double-pole
shape `J` gives another. A fifth rate does not fix this whole space.

There is also an angular blind direction for **any number** of energies
at the existing `t=u` angle:

```text
Delta=(s-t)^2(t-u)^2(u-s)^2,
Delta(2mu+v,0,2mu-v)=64mu^4 v^2-32mu^2 v^4+4v^6.
```

Thus more same-angle samples cannot replace an operator truncation and a
tail bound. A finite set of data can become sufficient after a justified
finite-dimensional restriction with a controlled remainder; no such
restriction is inferred here from a loop factor or the large heavy mass.

## 5. The complete gapped spectator increment is small after matching

[S305's full insertion](../problems/P8/s6/continuation/s6_305/notes/insertion.md)
separates a fixed Newton term, the local H/Proca curvature term
`-(log n+2)(S+12)/(640pi^2*kappa^2)`, and the nonlocal gapped remainder.
The first two are in the removed `(S,1,T)` span. Do not bound the large
raw Newton coefficient and then mistake it for a surviving remainder.

For all complex channel variables `|a|<=16`, the denominator of each
radial integral is at least `4nu-16`. With original `n>=10^12` and
Proca squared mass `m=10^6`, S305's exact radial moments imply

```text
|A_non| <= 16[(1/210)/(4n-16)+(3/14)/(4m-16)] < 10^-6,
|H_non| <= 16[(68/35)/(4n-16)+(36/35)/(4m-16)] < 10^-5.
```

The polynomial source bounds `|N2|<=344`, `|N0|<=108` and `pi^2>9`
therefore give the full three-channel bound

```text
|F_HV,non| < 3[344*10^-6/144+108*10^-5/6912]/kappa^2
           < 10^-5/kappa^2.
```

This also holds on the complex forward circle `|v|=1`, since the varying
channels `2+v,2-v` have modulus at most 3 and the transfer is zero. The
gapped integrals are holomorphic on and inside the circle. Cauchy's
coefficient estimate bounds their `v^2` coefficient by the same number.
Combining that coefficient and the four samples with `sum|w|<15` gives

```text
|L F_HV,non - sum_i w_i F_HV,non(q_i)|
    < 16*10^-5/kappa^2 < 10^-3/kappa^2 = 10^-1003*lambda.
```

This is the entire specified Gaussian H/Proca **increment at first loop**,
not an additional copy of S302, and not the massless M1 increment.

## 6. What is and is not in the current error ledger

| Contribution | Current bound or missing input |
|---|---|
| Retained first-order radiation conversion | `<10^-783*lambda`, already assembled in RATE4 v1 |
| Complete gapped H/Proca subtracted increment | `<10^-1003*lambda`, above |
| S302 metric contribution at the four samples | Weighted bound `<10^-989*lambda`; not a bound on its full functional |
| S305 massless M1 at the four samples | Weighted bound `<10^-998*lambda`; not a bound on its full functional |
| S293/S294 known crossing-center coefficients | Existing selected bounds; their full physical-sample subtraction still needs assembly |
| Heavy mass displacement | `<10^-388*lambda` at first order **if** a parent establishes `abs(delta_n)<=n/4` |
| Complementary local `R5` direction | `-902256*zeta`; no physical bound |
| Other complementary directions and omitted orders | No complete bound supplied |
| Full massless-gravity V/G functional | Its pole/cut/infrared prescription and contour error remain open |
| Same-parent curved/source/state response | Not determined by these flat conditions |

The metric sample estimate uses S302's `10^9/kappa^2` bound at `delta=1`
and the norm 15. For M1, the four centers have `2<=tau<=6` and
`|log tau|<2`, so S305 gives `|M1|<6/kappa^2`. Neither compact sample
bound permits a forward limit through massless cuts. The earlier
[S297 contour boundary](../problems/P8/s6/continuation/s6_297/notes/forward.md)
continues to apply.

Primary-source cross-checks support this boundary, not a numerical bound
for our candidate: [Alberte et al.](https://arxiv.org/abs/2007.12667)
explain the obstruction to standard spin-2 pole-subtracted positivity;
[Chang and Parra-Martinez](https://doi.org/10.1007/JHEP08(2025)175)
analyze loop-level forward singularities and use additional assumptions
for an infrared-finite four-dimensional construction. Their numerical
results or hypotheses are not imported into RATE4.

## 7. Curved matching and the next admitted task

Even fixing all flat shapes in a chosen truncation would not fix their
covariant lifts. The exact residue/`R H_heavy` degeneracy above and
[S336's independent curvature contact](../problems/P8/s6/continuation/s6_336/notes/matching.md)
are separate examples. S336 first changes a radiative amplitude at this
order, not the flat four-point amplitude. Its coefficient must not be
identified with the known S347 scalar-loop contribution. Nor does this
example prove that every such direction changes the homogeneous FLRW
background; B needs its own specified response projection and state.

The next task is **RATE4.PARENT: compare a complete additional matching
input**, before adopting one or returning to known-sector refinements.
It must supply either a common quantum-parent matching prescription with
controlled errors, or physically anchored conditions for a justified EFT
truncation together with a bound on its omitted projection. The latter
must address off-angle information, any retained heavy-mass freedom, and
the distinct curved/source data actually needed for B. It must also state
how the full V/G functional will be justified.

The minimum acceptance test is a bound on the **combined** complementary
projection and omitted terms within the remaining `lambda` budget. An
unqualified fifth rate, an arbitrary zero for zeta/chi, or a naturalness
guess is not that input. This is a specific research/design decision,
not a request for the user to invent Wilson coefficients and not a new
requirement to construct an all-orders UV theory. No additional parent
conditions have been adopted in this audit.

## 8. Reproduction and evidence level

```sh
.venv/bin/python -B scripts/p8_rate4_remainder.py
.venv/bin/ruff check scripts/p8_rate4_remainder.py
.venv/bin/ruff format --check scripts/p8_rate4_remainder.py
```

The [read-only diagnostic](../scripts/p8_rate4_remainder.py) pins RATE4 v1
and its inherited manifests, checking 308 protected files before and
after execution. It checks 24 independent Fraction/SymPy comparisons,
four half-line sensitivity bounds, the complete contact-basis rank,
mass-restored crossing identities, all 24 contact assignments, eight
first-order rate-null equalities, exact original-parameter majorants,
32 synthetic budget corners and ten invalid/incomplete-input rejections.
The gapped Cauchy step and the covariant/loop-order construction above
are written arguments using the pinned component proofs, not assertions
that finite fixture checks formalize quantum field theory.

Claim level: **VERIFIED_N** for this scoped derivation, identifiability
obstruction and bound. This is not a new numbered frozen packet, a full
physical matching certificate, a P8 closure verdict, or a repeat of the
115,487-test publication regression. The existing frozen results and
unrelated P4/P9 work are preserved.

The [validation receipt](validation/p8-rate4-remainder-2026-09-21.json)
records the actual diagnostic outputs, source hashes and successful
replays of the three preceding matching diagnostics.
