# P8 S6 — Conditional UV-consistency and matching contract

Adopted 2026-09-04 after the user delegated the choice of UV assumptions,
including additional research. This is a new conditional research branch,
not a revision of S0--S4 or a claim that S5 establishes a Wilsonian cutoff.
No UV completion is asserted. The original 32 linear row verdicts stand.

## 1. Recommended admissibility class

Investigate a **local, causal, unitary UV theory with a Lorentz-invariant
Minkowski vacuum**, coupled to gravity in the prescribed physical matter
frame. Lorentz invariance may be spontaneously broken on the bounce. This
vacuum requirement is an additional restriction of this branch, not a theorem
that every consistent cosmology must have such a vacuum. No string embedding,
swampland conjecture, naturalness prior or observational prior is imposed.

Treat the DHOST action as a low-energy EFT. Do not assume its finite derivative
truncation remains exact at arbitrary energy or is closed under quantum loops.
Higher operators/heavy fields are allowed in an explicitly recorded matching
model. Their effects on the bounce, degrees of freedom, causal cone and
interaction estimates require quantitative error bounds before that model
earns a verdict. They do not silently change the old function-group ladder.

The first computations retain the exact old action on the open clock tube.
An off-tube extension is permitted as a **separately named candidate**, with
the same disabled function groups and matter frame. Smooth coefficient
functions are allowed by the original contract. Global real analyticity of
those functions across X=0 and X=1 is a separate, stronger branch. Analyticity
of a vacuum scattering amplitude is not the same assumption. A smooth splice
is not presumed to arise from a UV theory: it must pass the matching gate.
Future corrections on the tube require a new candidate and a stated error
budget; none are authorized to overwrite an existing certificate.

## 2. Observables and order of limits

**V: vacuum scalar test.** Establish the vacuum equations, positive kinetic
residue and non-tachyonic spectrum first. Initially use one scalar with
physical mass m>0 to avoid conflating scalar massless limits with the graviton
problem. Verify an actual M->infinity limit that holds its canonical mass and
interactions fixed; deleting tensor diagrams on a rolling background is not
this limit. For its elastic amplitude require crossing, unitarity, analyticity
away from specified poles/cuts and sufficient complex-energy boundedness for
the subtracted dispersion relation. Subtract known light poles and, when
using improved bounds, the explicitly computed low-energy cuts. Strict
positivity requires nonzero absorptive weight; a tree-level zero is not alone
an exclusion. Every truncation/loop error needs a bound for a strict verdict.
Passing any finite set of necessary inequalities is not UV completion.

**G: finite-gravity test.** Keep the massless graviton. Our chosen conditional
route is a fixed-negative-transfer dispersion relation with an explicit
Regge description/remainder bound on the relevant high-energy contour.
Specify the infrared regulator or dressed observable, the light spectrum,
subtractions, contour and order of regulator removal for each calculation.
Combine the graviton pole with the dispersive/contour contribution *before*
t_transfer->0-. Do not impose a positive sign on a pole-subtracted coefficient
in isolation. Report the residual gravitational contribution R_grav and its
bound, including light-mass dependence; a presumed order-one coefficient in
1/(M^2 M_Regge^2) is not a certificate. If a controlled relation has the form
b2=I_positive+R_grav with R_grav>=-Delta_grav, the necessary bound is
b2+Delta_grav>=0. Without a justified Delta_grav it stays UNTESTED, not passed
or failed. The Regge route is conditional, not universal quantum gravity.

**B: connection to the bounce.** A vacuum test counts as evidence about the
bounce only after a common parent EFT/matching construction relates them.
Specify its field/derivative domain, state or observable, canonical maps,
heavy spectrum and cutoff; control omitted operators, loop corrections and
any threshold/branch changes on the relevant domains. A connected set of
coefficient functions, a finite-time hard tree bound, and a healthy vacuum
are each insufficient. A global scattering history connecting the two
background solutions is not assumed. A local-to-stationary amplitude argument,
if used, must bound its errors and establish the analyticity it needs.

Do **not** apply a flat-space inequality independently at every bounce time.
Boost-breaking amplitude analyticity is an extra condition, not a consequence
of the low-energy sound speeds. Cosmological correlator/time-delay routes
remain alternatives only after their hypotheses are matched to this bounce.

## 3. Verdict discipline and first deliverables

Keep separate outcomes: algebraic vacuum existence; leading scalar necessary
test; finite-gravity necessary tests; controlled bounce matching; full UV
completion (not furnished by positivity tests). An inapplicable theorem
excludes no witness. Failure of one specified extension excludes neither
other extensions nor its whole ladder row.

1. S6.1: determine how much a tube-fixed witness actually fixes at X=0.
   Construct a D-only smooth off-tube family; independently derive its
   canonical vacuum amplitude; compare with the globally analytic branch.
2. S6.2: construct/test a common matching model, starting with the positive
   vacuum member. Locate any failure of gap, kinetic or derivative control;
   record failures as failures of that matching ansatz only. Quantitative
   Regge/IR information is evidence to derive or parameterize, not another
   arbitrary sign to ask the user to choose.
3. S6.3: apply further scalar/mixed and finite-gravity dispersion tests only
   where applicable. M1 interaction control remains an independent P8 task;
   the D-only S6.1 calculation does not address it.

This contract resolves the previous *choice-of-assumptions* blocker. Missing
matching estimates are now research work, not a reason to ask the user to
choose a numerical correction. The assumptions are deliberately conditional
and may eventually exclude the selected ansatz without closing all of P8.

## 4. Evidence boundary

The source review is in [notes/literature.md](notes/literature.md). S6.1's
smoothness and analytic-identity arguments are written proofs; exact replay
checks their algebraic inputs and finite regression identities, not an
all-orders formal proof or the external dispersion theorems. New files live
under `s6/`, outside the earlier source globs. No old proof hash is rewritten.
