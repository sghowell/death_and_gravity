# Polynomial massive vacuum: exact tree realization and quantum preparation

## Model and fixed parameters

Use canonical +--- fields Phi,H in dimensionless units tau=1:

    L=(dPhi)^2/2+(dH)^2/2-V,
    V=Phi^2/2+mu^2 H^2/2+G H Phi^2/2+lambda4 Phi^4/24,
    lambda=10^-600, gamma=1024/10^800,
    D=2lambda/gamma, mu^2=D+2,
    G^2=2lambda D^3, G=+sqrt(G^2),
    lambda4=G^2(3D-2)/D^2.

These vacuum parameters are fixed; the label lambda4 is not the
derivative Wilson coefficient lambda. Restoring tau gives light
mass 1/tau, heavy mass mu/tau and cubic G/tau.
No cosmological action, state or counterterm is overwritten.

## Established result

For D>1 the polynomial potential has a unique global zero
minimum and positive canonical mass matrix diag(1,mu^2).
The selected D>2 puts H above the two-light-particle threshold.
The exact on-shell tree Phi-Phi amplitude equals the complete
S6.108 heavy-exchange amplitude; its actual light contact and
tree remainder therefore transfer at precisely that order.

At every fixed finite positive Euclidean regulator, H integrates
out exactly into a field-independent determinant and a nonlocal
quartic action for Phi. This is not a finite derivative expansion.
The action is coercive by a positive inverse-kernel bound.
A generic two-site action, its actual full Hessian and literal
four-by-four determinant check the exact Schur reduction.

Topological power counting determines the required local
counterterm structures at arbitrary loop order. None is higher
than quadratic in H, so the Gaussian structure is preserved.
This is formal perturbative renormalizability, not an all-energy
nonperturbative continuum existence theorem.

For constant Phi, evaluate the full one-loop Hessian on the
classical stationary H. The reduced kernel is positive at every
Euclidean momentum. After subtracting its constant, Phi^2 and
Phi^4 Taylor coefficients once at Phi=0, the resulting convergent
one-loop potential remainder is nonnegative and satisfies

    delta^3 Phi^6/[1536 pi^2 (1+A Phi^2)]
       <= V1_ren(Phi) <= A^3 Phi^6/(192 pi^2)
       < 10^-618 Phi^6   for Phi != 0,
    delta=lambda4-3G^2/mu^2 > 0,
    A=(lambda4-G^2/mu^2)/2.

At Phi=0 all bounds equal zero. The subtractions fix a
constant-background curvature and quartic, not the quantum
on-shell mass, residue or forward coefficient b2.

Retaining the momentum dependence gives the sharper bound

    V1_ren <= [(delta/2)^3+(G^2/mu^2)^3/mu^2] Phi^6/(48 pi^2)
           < 10^-815 Phi^6 for Phi != 0.

This does not require Phi to be small; a uniform relative
higher-loop error at arbitrary Phi is not inferred.

## What remains outside the result

The exact tree forward coefficient is 4lambda, but its full
renormalized loop correction and absorptive/cut errors are not
bounded here. A light loop cut lies inside the huge tree
polydisc; tree analyticity cannot be transferred to loops.
The bare unstable H pole must not be squared and integrated
through a dispersion contour without a consistent treatment.

This vacuum model is not a propagating common affine bounce
parent. Minimal Einstein coupling with its canonical kinetic
terms obeys the classical NEC and cannot itself give a flat
bounce. No full V, G or B verdict, original P8(b) classification,
absolute cosmological source or original P8 closure follows.
