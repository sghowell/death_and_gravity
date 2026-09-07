# P8: the actual off-shell reduction and its missing CD structure

This continues the [quantitative preparation audit](assessment-2026-09-07-p8-quantitative-preparation-cost.md).
The new [S6.30 reduction theorem](../problems/P8/s6/matching/variable/reduction/FORMULATION.md)
tests the original operator/matter dictionary directly. It does not
complete original P8.

## Main result

The unchanged variable-coupling two-metric parent does not reproduce
the original CD/M1 action through its leading stationary-metric
derivative expansion. Three distinct calculations establish what is
missing, without substituting a selected tensor solution for an
off-shell action match.

First solve the hidden metric's own algebraic equation, keeping the
physical g and both scalar fields fixed. Its positive proportional root
is f0=r(phi)^2 g with r^3=-beta1/beta4. The full ten-component potential
Hessian is invertible on the stated local positive branch, including
the trace and lapse components. Eliminating the first metric correction
gives the explicit four-derivative action

    kappa [B_mn B^mn-(tr_g B)^2/3],
    B=G[r^2 g],  kappa=M^4 r^3/(4 beta1).

The calculation retains every clock derivative and uses the physical-g
contractions. The constant-r limit reproduces the previously checked
physical-source Schur kernel. A direct canonical-matter variation fixes
the source sign as j=+2Pi_TT for the stated tensor convention.

After covariant integration by parts, the scalar-tensor part is on the
quartic-Horndeski locus A1=2F2_X, A3=0. The separate operator
kappa(Ricci^2-R_B^2/3) remains. Its coefficient is variable, so a
curvature identity involving the Euler density cannot discard it.

## Two independent off-shell mismatches

A literal tensor variation of the physical metric derives all
Christoffels, Ricci components and clock contractions independently.
The curvature-square term has a fourth-order tensor coefficient kappa;
every original quadratic-DHOST action is second order on the same test.
Thus the retained action is not exactly CD modulo boundaries and
scalar-only clock changes. This is not a ghost verdict on the parent:
a finite truncation and the complete differential system are different
objects.

The original clock is theta=tau*u, not the parent's canonical phi.
Applying the actual clock map and granting the favorable leading
normalization M_*^2=5M^2 gives, at the bounce,

    retained F2/M_*^2=-1/2,  retained F2_X/M_*^2=0,
    CD F2/M_*^2=-X_theta/2,  CD F2_X/M_*^2=-1/2.

Both have A1=0 there. The values agree at X_theta=1 but not their
coefficient jets on the open tube. The quantity
I=X(A1-2F2_X)/G_T is invariant under scalar-only clock changes and
constant normalization: it is 0 in the retained scalar-tensor normal
form and 1 for CD throughout that center tube. G_T here is not the
full response after including the retained higher operator.

Exact matching therefore requires an omitted F2_X coefficient of
magnitude M_*^2/2, and a weighted remainder
abs(Delta A1)+2abs(Delta F2_X)>=M_*^2. Even allowing ten-percent errors
in both normalized G_T and I leaves a necessary weighted remainder
of at least 81M_*^2/100. This is an explicit coefficient-norm test,
not a universal norm on every effective description.

## Why the small background lapse error is insufficient

With delta=c-2, the first correction gives the retained time-time metric
component 4c-4 at the center, versus actual c^2. Their difference is only
delta^2. Nevertheless the first correction's second u derivative tends
to 64, and the retained metric defect's fourth u derivative tends to 10752.
The latter is strictly larger than 10752 for every admitted c>2.

The proof checks the full nonlinear canonical-clock conversion as
well. It concerns the same physical-T metric component differentiated
as a function of the clock label, not a tensor component in a newly
transformed time coordinate. The clock-labeled fourth derivative has
the positive lower bound 107520000/5755201. Physical-time derivatives
carry the corresponding powers of 1/tau.

Thus delta alone does not produce a uniformly small C2 first metric
correction or C4 retained metric error on a fixed neighborhood. A
small value on a shrinking interval is compatible with this result;
the explicit inner expansion is recorded. Neither statement rules
out every differently controlled, nonadiabatic or finite-band inverse.

## Matter and the remaining matching obligation

An explicit metric redefinition that cancels the curvature-square
operator generates chi-curvature interactions. A further rewrite using
specified leading Einstein equations generates a positive Y^2 contact
as well as mixed clock contacts. Keeping the old free-M1 interpretation
requires retaining the full physical metric, matter and source map;
renaming the transformed metric does not satisfy that requirement.

The next omitted action coefficient is identified structurally, with
its constant-r six-derivative tensor symbol checked. Its full variable-
clock operator content and contribution to the missing CD structure
remain uncomputed. The target's inverse-X coefficients also prevent
naive total-derivative counting from being a hierarchy on X near 1.
Any higher-order repair must supply and bound the missing coefficients,
alongside the source/state and higher-operator errors.

This closes a concrete shortcut: background coincidence plus a large
center algebraic mass is not a controlled CD/M1 action match. It does
not close the original UV classification, prove a whole-row exclusion,
or forbid other parents. The scoped photon result remains unchanged.
Original P8 is open; no new user authority is needed for the next
discriminating calculations.

## Verification record

The complete 18-source child has been reviewed and frozen. Its freshly
rebuilt report has SHA-256
`1095385a6c9da4dded4e6a5234bae202457f9198fd1cd52967c91801906d5977`.
It records 138 exact identities, 12 positive algebraic margins and 26
rejected-domain controls. The written continuous proofs are complemented
by 68 independent Fraction profile/clock comparisons and separately
authored literal tensor, full-matrix potential and clock-jet audits.

All 97 new tests passed with ordinary SymPy in 134.97 seconds; the
separate ordinary certificate replay passed in 120.97 seconds. The full
P8 regression passed all 3,057 tests in 621.40 seconds. Only that broad
regression used the checked exact-GCD adapter: 6,509 exact descents and
4,552 ordinary-domain fallbacks, after its 128 ordinary-reference checks.
Both random seeds were fixed at zero; the host's problematic
faulthandler timer plugin was disabled. No cached ancestor verdict was
substituted for the ordinary certificate rebuild.

The result is CERTIFIED at the repository's written-proof plus exact-
algebra level, not FORMALIZED or an established controlled EFT.
Frozen ancestors and unrelated P4/P9 work remain untouched.
