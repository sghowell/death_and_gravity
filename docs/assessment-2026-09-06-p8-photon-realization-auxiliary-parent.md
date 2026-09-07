# P8: actual photon realization and auxiliary-parent applicability

Recorded 2026-09-06, following the
[cosmological-strength and symmetric-vacuum checkpoint](assessment-2026-09-06-p8-photon-cosmology-vacuum-branch.md).
P8(a)'s precisely stated photon/global-flat-FLRW objective remains complete.
The new photon example strengthens it without changing its prerequisites.
Original P8 remains open because the surviving DHOST rows do not yet have
an applicable controlled UV-matching/positivity verdict.

## A.19: a state and a geometry that solve the actual equations

The [new thermal construction](../problems/P8/a/fields/maxwell/thermal/FORMULATION.md)
uses the two physical transverse Maxwell polarizations. A positive
quasifree occupation covariance with `n(k)=1/(exp(b_T*k)-1)` is Hadamard:
its difference from vacuum is smooth by the ultraviolet exponential and
infrared integrability estimates. Conformal transport gives a state on
each interior point of the chosen smooth open flat-FLRW spacetime.

The thermal charge and physical temperature are

    Q=hbar*pi^2/(15*b_T^4),
    k_B*T_physical=hbar/(a*b_T), c=1.

Here b_T is a conformal-time length, not the photon finite counterterm
beta_M or the older scalar prescription parameter. The example explicitly
sets beta_M=0 and Lambda=0, with no added classical radiation or separate
nonzero gravitational curvature-squared term. The actual density retains
the anomaly:

    rho=Q/a^4+31*hbar*H^4/(480*pi^2).

Both independent Einstein components and stress conservation are checked.
This is not an order-reduced equation or a formal approximate solution.
With `delta=kappa*hbar/(8*pi^2*tau^2)` and `lambda=31*delta/180`, the
normal-clock variables `x=s/tau`, `y=-tau*H_normal` satisfy

    y'=2*y^2*(1-lambda*y^2)/(1-2*lambda*y^2),
    y(0)=2, a(0)=1,
    a^4=4*(1-4*lambda)/(y^2*(1-lambda*y^2)).

For every `0<delta<=10^-8`, exact continuous bounds place the entire
past interval `-1/100<=x<=0` strictly within the anchored A.18 radiation
C3 tube. The four response coefficients per lambda are
`(16/25,1728/25,155136/25,14770176/25)`; no sample scan or unidentified
asymptotic remainder supplies this embedding. Q is chosen from the
constraint and b_T from Q, so the matter state is part of the solution.

The regular low branch is `0<y<1/sqrt(2*lambda)`. Its exact clock
primitive gives a finite endpoint `1/16<x_end<1/4`. The endpoint scale
factor is positive for each fixed delta, but curvature diverges and the
same metric admits no C2 continuation through that endpoint. The
comoving geodesic therefore has finite remaining proper length. Neither
the endpoint nor another positive algebraic branch is admitted as smooth
initial data. Fundamental EFT validity at this high-curvature endpoint
is not proved.

This particular solution has positive effective energy density in the
focusing sense. Its explicit endpoint is consequently **not** a new
QEI-only argument. A.16–A.18 remains the independent all-Hadamard theorem;
this optional example provides an actual quantum/SEE realization of its
short-history radiation neighborhood. It proves neither generic-beta
existence nor observed-universe applicability.

## S6.13: vary the full auxiliary parent before claiming a vacuum

The [specified auxiliary action](../problems/P8/s6/matching/trimetric/FORMULATION.md)
has two dynamical vierbeine e,v and an auxiliary physical matter vierbein
u. With `Q=p_g*e+p_f*v`, its interaction density is

    -2*det(u)*(B+tr(u^-1*Q)).

The source-paper choice omits separate g/f cosmological terms. For a
constant regular triple the full dynamical equations include

    E_e=-2*p_g*det(u)*u^-T,
    E_v=-2*p_f*det(u)*u^-T.

If either link is nonzero, its equation cannot vanish. A source acting
only on u cannot fix that contradiction. In particular a constant
potential merely replaces B by `B+epsilon*V/2`; checking the u equation
alone would falsely accept a flat vacuum. The theorem does not divide
by B, use a heavy-mass formula around a nonsolution, or assume a
physical-time path from vacuum to bounce.

Where exact source-free elimination is available, `u0=-3Q/B` and the
effective potential has restricted geometric beta coefficients. These
do not equal the earlier beta2-only parent. A formal constant-potential
cancellation also fails on resummation: its nominally vanishing
coefficient is actually `9B/32`. Generic canonical sources change both
the actual matter metric and the Lorentz constraint at first order;
the frozen leading composite metric is not an exact substitution.

The action can genuinely be changed: add
`-2*beta4g*det(e)-2*beta4f*det(v)`. For positive q, the choice
`p_g=p_f=q`, `B=-6q`, `beta4g=beta4f=-q` has `e=v=u=I` as a full flat
vacuum with zero total density. Its relative flat quadratic mass is
`q*(1/G+1/F)>0`. The physical metric is h=eta, with no factor-four clock
rescaling. This countercontrol prevents generalizing the no-vacuum
result to all auxiliary parents. Rolling matching, stability and cutoff
control for the extension require separate work.

## Current literature check: no replacement for the matching obligation

The existing [S6 source audit](../problems/P8/s6/notes/literature.md)
already distinguishes vacuum from cosmological dispersion hypotheses.
Three additional primary-source checks sharpen the next-step decision:

- Ye–Piao's cosmological amplitude decomposition, published section 3.2,
  equation (20), remains an assumption with unknown background corrections
  and omitted particle production. Its section 3.1 removes the graviton
  pole through compactification. These do not provide the adopted
  finite-gravity or rolling matching certificate.
  [Published paper](https://link.springer.com/article/10.1140/epjc/s10052-020-7973-z)
- Bellazzini et al. provide a genuinely useful new IR-finite,
  detector-resolution-dependent amplitude framework. Sections 4.2 and 5,
  equations (4.10), (5.13), and (5.17), retain a scaling-limit functional
  unitarity condition, all-spin smearing positivity and finite-coupling
  `O(G*M^2)` corrections. Those corrections are not a supplied numerical
  error bound for this bounce. The framework is a candidate observable
  route after matching, not a current row verdict.
  [arXiv:2512.13780v2](https://arxiv.org/html/2512.13780v2)
- The newer emergent-geometry bounce uses a wider acoustic cone and an
  Einstein lapse passing through zero; see equations (2), (5), (13) and
  the discussion following (15). It changes the regular, prescribed-frame
  subluminal contract rather than supplying a witness within it.
  [arXiv:2602.21642v2](https://arxiv.org/html/2602.21642v2)

These are applicability findings, not refutations of the papers.
The useful next calculation is whether the surviving auxiliary extension
can match the physical tensor/matter cones once the actual source is
retained. Adding more vacuum inequalities alone cannot decide the old
rows: S6.1 already preserves the entire D tube while changing the vacuum
scalar diagnostic. No unproved UV assumption is selected to force closure.

## Verification

| Report | SHA256 | Hashed files | Exact residuals | Targeted tests |
|---|---|---:|---:|---:|
| [A.19](../problems/P8/a/fields/maxwell/thermal/certificates/thermal-see.json) | `d77f1a4ef293c53c70bb7fd540b473c557c36018dfbea558ea1af9afc147375c` | 17 | 27 | 66 |
| [S6.13.AUXILIARY](../problems/P8/s6/matching/trimetric/certificates/auxiliary-parent.json) | `062ea71b4fcb139af3ea727f013eeeb3a340b6e37207cb872e372eacce09631e` | 15 | 49 | 33 |

Both new targeted suites and standalone read-only certificate CLIs pass
with the ordinary interpreter. The separate root audits have twelve
tests each. S6.13 also replays 64 independent Fraction full-matrix
derivative directions and three coefficientwise polynomial identities.
All 32 new hashes and pinned ancestor replays are checked without writes.

Together with the previous checkpoint's 1,753-test adapted exact run,
the 99 ordinary new tests cover all 1,852 tests in these checkpoints.
This is coverage across runs, **not** a claim that an ordinary combined
1,852-test command passed. See the
[exact-runner audit](p8-exact-regression-runner.md). Subsequent cone work
is outside this checkpoint's inventory and verdicts.
