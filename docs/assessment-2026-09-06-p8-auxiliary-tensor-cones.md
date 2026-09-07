# P8: the actual matter cone rules out the symmetric auxiliary route

Recorded 2026-09-06, following the
[actual-photon and auxiliary-vacuum checkpoint](assessment-2026-09-06-p8-photon-realization-auxiliary-parent.md).
This continuation studies the positive-link auxiliary extension that
survived the previous flat-vacuum screen. It does not change any old
certificate, the prescribed physical matter metric, or P8's closure scope.

## What the new theorem adds

The [S6.14 formulation](../problems/P8/s6/matching/trimetric/cones/FORMULATION.md)
retains the actual matter metric `h=u^T eta u`, the full auxiliary source,
positive Einstein coefficients and positive link coefficients. Separate
endpoint cosmological terms are allowed. For homogeneous coframes define

    P_g=p_g*a_e/a_u, P_f=p_f*a_v/a_u,
    c_e=a_u*n_e/(n_u*a_e), c_f=a_u*n_v/(n_u*a_v).

These c values are gravitational tensor propagation speeds in physical
h units. The full u equations, varied before the FLRW restriction, imply

    P_g*(c_e-1)+P_f*(c_f-1)=n_h/2,
    n_h=epsilon*(rho+p).

Positive link weights and strictly positive actual null stress require
`max(c_e,c_f)>=1+n_h/[2*(P_g+P_f)]>1`. Neither B nor the endpoint terms can
change this identity. It uses no Hubble quotient and remains valid at a
bounce. Only at least one full-parent cone is constrained, not both.

That result alone would not exclude an asymmetric light-only EFT. A
positive quadratic countercontrol has high-frequency squared speeds
`1/4,4`, but its heavy-relative limit has light squared speed `29/101`.
This is a logical countercontrol, not a full auxiliary-parent solution.

The exact exchange-symmetric branch gives the stronger useful conclusion.
With `e=v=r`, equal Einstein coefficients G, equal positive links q and
equal endpoint terms, write `u=diag(N,A,A,A)*r`, N,A>0. Then

    c_T=A/N=1+A*n_h/(4q)>1,
    G_T,h=2G*N/A^3, F_T,h=2G/(N*A).

Both kinetic and gradient coefficients are positive, but their ratio gives
a faster-than-matter cone. The common tensor is physical and decouples
from the exchange-odd tensor even on the time-dependent background. It
therefore survives elimination of the relative tensor. A large relative
mass cannot by itself fix this common channel.

## The source and the clock are part of the result

For the literal external probe action `integral sqrt|h| j*gamma_u`, exact
auxiliary elimination gives

    gamma_u=gamma+A*j/q,
    L_probe,eff/sqrt|h|=j*gamma+A*j^2/(2q).

The local contact is retained. It does not remove the propagating common
response. The derivation uses a unit-normalized TT polarization and the
explicit physical-stress dictionary; the probe is not the homogeneous
canonical scalar source. Homogeneous canonical matter has no independent
linear TT anisotropic response. Merely isotropic background stress would
not establish that condition for a general medium.

The actual time conversion is `dT=N*n*dt`, with physical scale `a_h=A*a`.
Dropping it would give the wrong light cone. Canonical tensor normalization
also retains the time-dependent pump term and boundary term. A positive
relative algebraic mass is not promoted to a proved adiabatic gap or a
stationary scattering pole.

The two project tracks use different curvature conventions. Here P8(b)
has `R_B=-6*(Hdot+2H^2)`, `G_B00=+3H^2` and Einstein action `-G*R_B/2`.
P8(a)'s FK signs cannot be copied into this action without conversion.

## A genuine rolling solution and controls outside the theorem

The same positive-link beta4 extension has a nonempty local rolling
family. With B=-6q, both endpoint coefficients -q and a free canonical
scalar, choose common-r proper time and any A0>1. Then

    rho=p=12q*(A-1)/(epsilon*A), N=A/(6A-5),
    H_r^2=2q*(A^3-1)/(3G), H_r>0,
    A'=-6*H_r*A*(A-1)/(6A-5), a'=H_r*a,
    psi'=N*sqrt(2rho).

The analytic ODE supplies a local solution, and all three full coframe
Euler maps and the actual scalar current equation vanish. At A0=2 the
physical tensor speed is 7. Arbitrarily near-vacuum rolling points also
exist, with speeds above one tending to one. The example is expanding,
not a CD bounce, a global existence theorem, or a scalar/vector/cutoff
certificate. Its role is to show that the theorem concerns real sourced
solutions, not an inconsistent algebraic fixture.

Two mixed-sign controls keep the assumptions explicit. The action
`p_g=-2,p_f=1,B=3,beta4g=2,beta4f=-1` has a genuine flat vacuum with
`e=v=u=I` and positive relative quadratic spring 4. Thus a healthy flat
tensor mass does not force positive links. A separate mixed-sign
auxiliary-only point obeys both u equations with positive null stress
and speeds `1/2,1`; it fails the dynamical e/v equations and is labelled
accordingly. Neither control is a full healthy rolling UV construction.

## Quantitative meaning and remaining closure boundary

For `delta=A*n_h/(4q)>0`, the exact coefficient gap is
`F_T,h-G_T,h=G_T,h*delta*(2+delta)`. If a separately named correction has
physical coefficient errors bounded by E_G,E_F and corrected kinetic
coefficient remains positive, repairing subluminality requires

    E_G+E_F >= G_T,h*delta*(2+delta).

This is a necessary correction threshold, not an estimate of actual
omitted operators or loops. Keeping the old canonical free-chi bounce
velocity and additionally identifying the vacuum Planck/mass scales
gives `c_T-1>=A/[100*(m_0*tau)^2]`. Changing a bookkeeping epsilon cannot
erase a fixed canonically normalized source. No old background or field
dictionary is silently identified with the one-scalar example.

The new result rejects the exact symmetric canonical auxiliary route
against P8's **prescribed subluminal matching requirement**. It is not
a universal UV-inconsistency theorem inferred from superluminality.
Mixed links, asymmetric controlled light reductions, new operators and
quantum mechanisms remain outside the corresponding hypotheses.

| Original obligation | Current evidence |
|---|---|
| P8(a), photon/global-flat-FLRW specialization | Scoped objective complete; A.19 is an additional actual-state realization |
| P8(b), frozen 32-row linear-principal classification with prescribed canonical matter | Complete within the stated operator/tube/tail contract |
| Applicable vacuum and rolling matching for a surviving UV candidate | Not established; additional named routes now excluded, not whole rows |
| Controlled finite-gravity positivity for that common candidate | Not established; no numerical gravitational/IR/loop remainder is invented |

P8 therefore remains open. The original problem is not completed by
changing the closure criterion to count an inapplicable positivity test
or a failed selected parent as a classification of every UV possibility.

## Verification

The [S6.14 report](../problems/P8/s6/matching/trimetric/cones/certificates/actual-tensor-cones.json)
has SHA256
`ab280acf8d20fb2ae18d75c265deecd59832e993ba0338e9e5c8f7691c0a7e6e`.
It pins 15 source, proof and test files, with 28 main exact residuals,
four independent coefficientwise Fraction identities, two independently
reconstructed physical-clock pullbacks and six bridges to the primary
calculation. Two literal source-jet directions and four rolling fixtures
are also replayed. The separately authored audit contributes 16 tests.

All 46 new targeted tests pass with the ordinary interpreter, independently
of the opt-in exact GCD adapter. The standalone read-only certificate
replay, source hashes and Ruff checks also pass. The combined P8 regression
passes **1,898 tests in 424.23 seconds** using the documented
[exact regression runner](p8-exact-regression-runner.md): 6,509 exact
Gaussian-to-real GCD descents, with every returned cofactor identity
checked. This is an adapted exact run, not an ordinary unmodified full-suite
claim. Its diagnostic long-test stack dump was followed by successful
completion; no skipped or timed-out test is counted as a pass.
