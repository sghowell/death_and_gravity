# P8: quantitative preparation cost and the low-frequency limit

This continues the [fixed-source and local-vacuum audit](assessment-2026-09-07-p8-fixed-source-and-local-vacuum-obstructions.md).
The new [S6.29 cost theorem](../problems/P8/s6/matching/variable/response/prepared/cost/FORMULATION.md)
turns existence of a selected regular tensor sector into explicit source
and frequency budgets. It does not complete original P8.

## What is now quantitatively controlled

Keep the actual full S6.20 two-metric action, its g-metric matter frame
and its conserved external tensor source. All four physical initial
data vanish at u=-1/50. The source is applied on
J=(-1/50,-1/100), and the spatial momentum parameter K lies in [1,4].
Only the preparation interval's punctured coefficients are evaluated
at delta=0 to choose the source. The same source is then used for every
positive delta; there is no later reset of relative-mode data.

The output is measured in the full Frobenius-normalized wave frame.
A regular light target is Lv with the specified two-component origin
normalization. It is not raw canonical outer data or unit physical
g-metric data. These distinctions enter the cost proof.

The full four-dimensional loading Gramian has the continuous enclosure

    9/(4*10^12) I <= C_K C_K* <=16/25 I.

The lower bound is obtained by a validated matrix-ODE calculation at
K=5/2 followed by an explicit continuous-momentum perturbation theorem,
not by a parameter grid. At the midpoint, the code uses 64 adjacent
closed Taylor steps, order18, and outward 256-bit Arb arithmetic.
A whole-step state tube and normalized derivative recurrence bound the
actual Taylor remainder. Four strictly positive pivots certify the
shifted Gramian; exact Fraction intervals independently repeat the
positivity calculation.

The continuous extension uses a scaled canonical chart whose heavy
rotation is skew. Its duration-resolved energy bound is
exp(11(t-s)/20), its momentum-generator derivative is below1/5000,
and its source norm is below1/2000. These give a source-map variation
at most3*10^-6 across the full momentum interval. The actual wave map
then transfers the Gramian bound to the true target coordinates.

For any full target y, the minimum L2 source exists and satisfies

    (5/4)||y|| <= s_min(y) <= (2*10^6/3)||y||.

Actual symplectic source moments strengthen the lower bound to
19||v|| for every nonzero regular light target, and to900 for the unit
even target. These are sufficient lower bounds, not near-optimal costs.

## A concrete source, and the distinction between minimum and infimum

A degree11 polynomial physical f field matches six endpoint jets on
each side. Recovering g through the unsourced f equation and then
recovering sigma through the actual sourced g equation prepares all
four target data exactly. The source and its first derivative vanish
at both endpoints. Its zero extension is C1 and H0^2, but is not
claimed to be C-infinity with support strictly inside J.

For a nonzero target y, the explicit construction has

    ||sigma||2 <2*10^8||y||,
    ||sigma_uu||2 <4*10^14||y||.

These bounds retain all four source-operator derivatives, both clock
rescalings and the true endpoint map. An independent rational engine
reconstructs the operator by composition, solves the Hermite conditions
in a power basis and converts to Bernstein bounds. It agrees with the
separate primary calculation over the entire K interval.

The spectral optimization is posed with an L2 source budget but an
H0^2 source class. That feasible set is not closed in L2: its optimal
value need not be attained. The proof gives the exact unique L2
optimizer and a one-parameter dual. For every budget strictly above
the minimum L2 control norm, smooth approximation with an exact
endpoint correction proves equality of the H0^2 infimum and L2
minimum. At the exact minimum budget, H0^2 feasibility must be checked
separately. A generic constant minimum-norm control illustrates why
the endpoint traces cannot be presumed.

## What the source spectrum does and does not allow

The Fourier transform is applied to the actual physical-clock source
sigma, extended by zero, with energy measure domega/(2pi). It is not
applied to a weighted canonical forcing component. For any temporal
band |omega|<=Omega<=100,

    E_in <= (1/3)||sigma||2^2,
    E_out >= (2/3)||sigma||2^2.

This is a statement about a squared source norm, not positive matter
energy. In dimensional variables the frequency is omega/tau and the
squared norm is multiplied by M^4/tau^3. The ratio is unchanged.

For unit regular v and source budget S=10^6, the unconstrained L2
spectral optimizer has squared norm at most2*10^12/3, strictly within
that budget. Its optimal tail cost, equal to the H0^2 infimum, is

    722/3 < inf E_out <=4*10^12/9,

with the stronger lower bound540000 for the even target. The broad
interval is not presented as a sharp numerical optimum. The result
does not compute the optimizer or the actual low-pass loading moments.

An exact sinc expansion gives a further, replayable interface: nine
actual loading moments suffice for a finite-rank band Gramian whose
operator error is at most1/(3*11!). Those moments remain uncomputed.
A Fourier-low-pass projection has a noncompact past and future; using
it with newly imposed zero data inside J would change this problem.

## Fixed-source propagation and remaining P8 gates

For a fixed source with norm at most S and limiting normalized loading
error e, the unchanged positive-delta system has output error at most

    42e+delta(8600||v||+12600000 S)

relative to Lv. The comparison with the exact analytic prepared target
has8400 in place of8600. The estimates extend from smooth sources by
L1/Duhamel continuity, and include actual delta-dependent endpoint
maps. They are not raw g-field derivative error estimates.

For unit v, S=10^6, e=0 and delta=10^-17, the bound is below1/1000.
It is not small for all delta<=10^-9 at that source budget.
The explicit polynomial source has a larger sufficient budget and a
separate, correspondingly smaller delta example. None of these
dimensionless restrictions is called an EFT cutoff.

The matching advance is therefore specific: this selected sector has
controlled, quantitatively bounded preparation, with its spectral
limitations exposed. It is not a statement that generic slow sources,
a vacuum or arbitrary initial data select that sector. It supplies
neither a rolling gap nor the actual off-shell C/D operator and matter
dictionary. Local-vacuum absence, finite-gravity dispersion remainders,
omitted operators/loops and complete-background health remain separate
obligations. The completed scoped photon objective is unchanged.
Original P8 remains open; no new user authority is required to continue.

## Verification record

The frozen report is
`7c0795c9333641039c5737a42f0dee60859af4c132ff017dc3972f1c0136f596`.
Its21-source manifest includes119 exact residuals,29 strict continuous
margins,256 literal rational-to-Arb coefficient bridges,291 independent
Fraction source/Hermite comparisons,18 independent energy comparisons
and34 rejected-domain controls. The two separately authored scientific
audits cover the physical source and validated Gramian arguments.

Fresh ordinary verification passed126 tests in156.69 seconds; a separate
ordinary certificate replay passed in144.05 seconds. The full P8
regression passed2960 tests in585.02 seconds with the reviewed exact-GCD
adapter and its128 self-checks,6509 exact descents and4542 ordinary-domain
fallbacks. Both hash and SymPy factor RNG seeds were zero; the host's
problematic faulthandler watchdog was disabled, with no mathematical
assertions omitted. Ruff and whitespace checks are clean. Frozen
ancestors and unrelated P4/P9 work remain untouched.
