# P8: first closure-driven matching decision audit

Recorded 2026-09-20. This is an assessment and elementary algebraic
cross-check of existing results, **not a new numbered theorem/certificate**.
It implements the first screen in the [closure plan](p8-closure-plan.md).
Original P8 and V/G/B remain open.

## Decision

The present known-sector calculations do not determine the finite matching
data needed for a model-specific positivity verdict. At least two already
identified finite directions affect its local coefficient with nonzero
sensitivity. Completing another known scalar-loop contribution cannot
bound those independent directions.

This establishes which input to prioritize. It does not prove that the
selected parent is inconsistent, that admissible UV theories realize
arbitrary coefficients, or that a particular physical coefficient is free
once its microscopic parent has been specified.

## Exact local sensitivity screen

[S302's formulation](../problems/P8/s6/continuation/s6_302/FORMULATION.md)
and its [actual assembly function](../problems/P8/s6/continuation/s6_302/src/p8_vacuum_affine_minimal_gravity_finite/assembly.py)
retain the finite contribution

    [alpha*(s^2+t^2+u^2) + beta*mu^2
       + 16*pi^2*delta_kappa*T0] / (16*pi^2*kappa^2).

Here mu is the light mass squared. At fixed t set
s=2mu-t/2+v and u=2mu-t/2-v. The coefficient of v^2 in the regular
local polynomial is exactly

    delta b20_alpha = alpha/(8*pi^2*kappa^2).

Beta contributes zero to this coefficient. The Newton/pole term has been
held separate, **not set to zero or counted as another local contact**.
Taking the coefficient of this polynomial requires no massless forward
limit of the full amplitude and establishes no physical dispersion theorem.

[S290's forward functional](../problems/P8/s6/continuation/s6_290/notes/forward.md)
gives, for the heavy mixing term f2(a)=h/(n-a),

    b20_h = -d_a^2[(a+2mu)*h/(n-a)] at a=2mu, divided by kappa
          = -2h*(n+2mu)/[kappa*(n-2mu)^3].

For the independent finite R H direction, h=-2g*c_RH. Thus

    delta b20_RH = 4g*(n+2mu)*c_RH/[kappa*(n-2mu)^3].

This agrees with the existing
[S294 source identity](../problems/P8/s6/continuation/s6_294/src/p8_vacuum_affine_whole_mixed_heavy_gravity_sector/cuts.py)
and [matching qualification](../problems/P8/s6/continuation/s6_294/notes/matching.md).
Both sensitivities are strictly positive for kappa,g,mu>0 and n>2mu,
including the original parameters.

For any fixed real known coefficient K, the formal expression K+w*alpha
with w=1/(8*pi^2*kappa^2)>0 attains either sign if alpha is left unrestricted.
For example alpha=(1-K)/w and alpha=(-1-K)/w give +1 and -1. These are
**algebraic insufficiency controls**, not proposed physical values. They
need not preserve perturbative control or be realizable by a UV parent.
Their role is to prevent a bound on the known term from being reported as
a bound on the sum without the missing coefficient assumption.

This is only a two-direction projection, not an exhaustive counterterm
basis or the complete b20. The physical soft/detector conversion, other
finite terms, higher orders and gravity contour contribution remain
separate. Individual known pieces are prescription-dependent; only a
consistently matched observable can support a physical verdict.

## What a useful matching bound would enable

Suppose a future physical matching argument supplies valid enclosures
alpha in [a_-,a_+] and c_RH in [r_-,r_+], a known contribution in [K_-,K_+],
and a total residual error at most E for the *same* target coefficient.
Then the positive weights above give the conservative enclosure

    b_- = K_- + w_alpha*a_- + w_RH*r_- - E,
    b_+ = K_+ + w_alpha*a_+ + w_RH*r_+ + E.

The rectangular enclosure is conservative if the two matching coordinates
are correlated; no independence of physically realizable values is assumed.
If instead only their combination C_pos is needed, bound it directly.

Only after the G hypotheses are independently justified and an allowance
Delta_grav>=0 is established does the necessary condition become

    b20 + Delta_grav >= 0.

At that stage:

- b_+ + Delta_grav < 0 excludes the tested candidate/domain under those
  hypotheses;
- b_- + Delta_grav >= 0 passes this necessary inequality, not UV completion;
- otherwise the enclosure is inconclusive.

Currently the matching intervals, complete residual E and gravitational
allowance have not been supplied. This is a **decision specification**, not
an evaluated P8 inequality. It cannot be populated with guessed naturalness
scales or the smallness of the known loop pieces.

## Why S347 and the fold results do not resolve this screen

[S347](../problems/P8/s6/continuation/s6_347/notes/scope.md) determines the
selected scalar-loop degree-six real-TT curvature coefficient in one
off-shell jet convention. Its normalized magnitude below 10^-207 does not
bound alpha, c_RH or the independent added parent chi. Nor does a real-null-TT
projection supply the full Ricci/trace/off-shell metric information needed
for a curved background.

[S336](../problems/P8/s6/continuation/s6_336/notes/matching.md) independently
shows why retained flat and soft data leave a hard radiative matching
direction undetermined. It explicitly treats parent matching as research,
not as a user choice of chi=0.

[S276](../problems/P8/s6/continuation/s6_276/README.md) excludes a global C1
extension of a particular lapse sheet. [S277](../problems/P8/s6/continuation/s6_277/FORMULATION.md)
excludes smooth passage through a specified classical fold family and
constructs different classical endpoint solutions. Neither supplies a
no-go for the original prepared quantum state or refutes the corrected
S275 finite hybrid. They cannot be promoted to a shortcut closing B.

## Keep the gravity target finite where justified

A primary-source recheck supports retaining the graviton pole and its
contour contribution together: their finite remainder need not be positive.
See [Tokuda, Aoki and Hirano, sections 3.1-3.2](https://arxiv.org/html/2007.15009v2).

An important scope economy is available: the finite-energy sum-rule
construction in [Caron-Huot and Tokuda, section 2.1](https://arxiv.org/html/2406.07606v2#S2.SS1)
requires control on its finite complex-energy contour, rather than a proof
of an asymptotic Froissart-like bound. Their assumptions still need to be
verified for P8, with an absolute remainder; the result is not automatically
transferable to this parent. This motivates a bounded contour target, not
an assumption that the missing high-energy contribution is small.

The existing [S298 conditional budget](../problems/P8/s6/continuation/s6_298/notes/regge.md)
has not supplied its U,L,epsilon,B inputs and bounds only a selected part of
the sum rule. It is not yet the required total Delta_grav.

## Validation and next action

A fresh exact SymPy check on 2026-09-20 confirmed the crossing-local
coefficient, zero beta sensitivity, heavy-residue differentiation, positive
weights on n=2mu+d with d>0, and both formal sign controls. Original SymPy
was used without the full-regression adapter. This elementary cross-check
does not replay the parent certificates or verify the physical hypotheses.
All 16 local links in the two new documents were checked, and all 251
source hashes recorded by the S336-S347 certificates still match. The
unrelated seven P4/P9 files also match the pre-existing preservation record.

Next: **MATCH-1**, a common-parent enclosure of the decision-relevant finite
combination, with an explicit physical input and residual budget. If the
same missing input remains unbounded, label the route conditional and
compare a concrete alternative; do not continue an unbounded series of
known-loop refinements. No frozen scientific file, physical parameter,
claims-ledger result or P8 closure status is changed by this assessment.
