# Entire connection-source pullback and instantaneous contact

The formal source used to define connection observables is distinct
from the retained physical Proca/heavy source. Require projective
invariance G^T J=0 for the original 64-component source.

At a fixed retained history and finite cell write

    Gamma=Gamma_bar+E K z,       b=(E K)^T J.

The entire source term is J^T Gamma_bar+b^T z. It cannot be replaced
by a source linear in a preferred subset of the new fields.
With the source-independent normalized complement density and
the continued Fresnel phase, completion of the full square gives

    Z_complement[J | retained history]
      = exp[i J^T Gamma_bar-i b^T A^-1 b/(2lambda)].

The source-dependent auxiliary equation is lambda A z+b=0, whose
root is -A^-1 b/lambda. The positive constraint density and delta
Jacobian still cancel on this regular root. In the reduced action,
however, the full quadratic source term remains. Normalizing it
away with a source-dependent factor would change the observable.

The source-functional conditional mean and connected contact are

    <z>_J=-A^-1 b/lambda,
    <zz^T>_J,connected=i A^-1/lambda.

These derivatives include an instantaneous contact. They are not
a new propagating mode, nor an assertion that a Dirac-reduced
zero operator has ordinary time-ordered fluctuations. They refer
to the explicitly specified source-dependent functional, with its
changed auxiliary constraints and contact convention.

For the original full connection the contact kernel is

    (i/lambda) E K A^-1 K^T E^T.

Its entire 64-by-64 matrix is retained. The 60-quotient inverse
decomposes as

    Mnew^-1=K A^-1 K^T+L eta L^T.

The first term has zero retained trace. The second belongs to the
retained trace mass and is not used to integrate out the dynamical
Proca field. Both full Euler identities are checked.

## After averaging over retained fields

If a retained finite functional, with its state and ordering,
has separately been specified, differentiate the conditional
identity before performing its retained average. Then

    <Gamma>=<Gamma_bar>,

and the two-point source-functional derivative contains both
the full composite Gamma_bar correlations and the finite-cell
delta contact i<lambda^-1 E K A^-1 K^T E^T>.
There is no general identity Gamma_bar(<Q>)=<Gamma_bar(Q)>.
The independent centered Gaussian example <q^2>=variance while
<q>^2=0 is a simple countercheck, not a new P8 state.

## Derivative coordinates and the remaining momentum equation

In the finite-jet diagnostic, an original source J^T y becomes
J^T(z+f(q,v)). After eliminating z the whole retained Lagrangian is

    Lred+J^T f(q,v)-J^T A(q,v)^-1 J/2.

The remaining kinematic momentum equation is obtained by
differentiating this full expression with respect to v.
The inverse contact and f_v are not optional. In the explicit
f2=v^2-q example, the second velocity derivative of the term
linear in J2 is 2. That source-dependent contact would be lost
by imposing the unforced constraint before variation.

The same rule applies to the original connection's full
metric/clock jet-dependent center and to its unchanged
T=W+B(u,X)du shift. All earlier fixed profiles, physical source
squares, temporal boundaries and CTP preparations remain.
No mean-field factorization or solved interacting mean follows
from the conditional Gaussian identity.
