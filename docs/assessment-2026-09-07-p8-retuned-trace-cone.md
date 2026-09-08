# P8: mass retuning fixes kinetic signs but fails the matter-cone test

Original P8 remains open. The scoped photon objective, linear CD/M1
classification, original matter frame and adopted V/G/B contract are
unchanged. This follows the
[curved Ricci-difference checkpoint](assessment-2026-09-07-p8-ricci-difference.md).

## Exact new result

[S6.41](../problems/P8/s6/matching/affine/kinetic/retuned/FORMULATION.md)
adds a separately named source-centered mass term for T=V+U, with
mu=11/9, before testing its positive curl. The center Tstar is the
complete original stationary connection trace, not an external profile.
The linear source and Tstar² counterterm are both retained.

At zero curl, all 64 stationary equations and the complete original
CD/M1 reduced action match exactly. The new projective quotient has
rank 60 throughout the original tube and an explicit infinity-norm
inverse bound below 6900. A dense independent inverse and determinant
check supplements the generic block and Woodbury identities.

The derivative-bearing trace source aligns with the original scalar
kinetic null. In the exact nonlinear spatial chart it is

    Tstar_normal=(4p²-1)*K_hat+6p*J3*s³.

There is no lapse velocity. A uniform relative bound below 19/1080
also keeps the full ten-velocity metric/vector/matter block at rank.
This is a primary-degeneracy result, not a complete nonlinear Dirac
constraint theorem or an inference from an unreduced kinetic sign.

The actual rolling source is Tstar_0=-d*n, d=15H/(4h), with no spatial
component. After retaining and solving lapse, shift, vector temporal
and free matter, the three scalar kinetic directions are positive.
An exact polynomial bound gives positive principal gradient energy
for 0<zeta<=1/2000 at every finite u!=0 and q>0. The regular first-order
center uses J, not Theta, as its lapse denominator; the q>6 center
velocity chart and both transverse vector polarizations also pass
their sign checks.

## Why this is not an accepted witness

The original acceptance condition requires K-G>=0 as well as positive
K and G. It is not a new requirement. The original two-scalar witness
has K0=G0 exactly, so both modes saturate the physical matter cone.
For every nonzero curl, the complete new high-frequency characteristic
polynomial is

    (1-c²)*[(1-c²)²-kappa*c²], kappa=d²/(2J)>0 (u!=0).

Thus one scalar remains luminal and the other two have

    c_plus²=1+kappa/2+sqrt(kappa²+4kappa)/2 > 1,
    c_minus²=1+kappa/2-sqrt(kappa²+4kappa)/2 in (0,1).

Their product is one. A separate two-dimensional principal minor of
G-K has determinant -f²/4<0, with f=d/Theta. These are physical modes
after all constraints, not the phase speed of a massive field at a
finite momentum. The curl cancels from the principal limit; setting it
to zero first instead removes the extra mode and is a different rank
limit. A general positive isotropic retained mass changes the mixing
strength but not this mechanism.

The literal candidate therefore fails the original matter-cone gate
despite favorable energy signs. This is not an exclusion below a
justified EFT cutoff or a general metric-affine/UV no-go. It supplies
no coupled heavy gap, nonlinear health theorem or original P8 closure.

## Next calculation

The next candidate should address the source mixing explicitly rather
than infer causality from a positive isolated vector mass. A concrete
local test is to curl W=T-B(phi,x)*dphi instead of T, with
B=-d(phi)*(1+x)/2. This is a new operator, not a reinterpretation of
the rejected one. The preliminary exact normal-source identity is

    Wstar_normal=(4p²-1)*(K_hat-3H*s)+(3/2)*s*Q_lower.

The unchanged lower-coefficient ODE has Q_lower=Q_lower,x=0 on the
clock. This motivates checking whether all first-order mixing vanishes,
then deriving the actual canonical vector equations and nonlinear
source bounds. The new [working scope](../problems/P8/s6/matching/affine/kinetic/aligned/FORMULATION.md)
does not yet carry a certificate, a heavy-sector verdict or a UV claim.
No user intervention or new external authorization is needed for it.

## Verification

The frozen report covers 18 sources, 48 named exact identities comprising
4898 scalar entries, 24 continuous/interface proof checks and 34 rejected
inputs. The final scientific suite passes 70 tests in 25.19 seconds.
The full ordinary suite passes **89 tests in 143.64 seconds**, without
the broad GCD adapter, and the separate seeded CLI passes.

The full **4241-test P8 regression passes in 724.81 seconds**, with no
certified checkpoint excluded. The exact GCD adapter passes 128
original normalized-tuple comparisons and records 5597 domain fallbacks
and 6509 exact descents. This uses the unchanged per-test SymPy RNG
seed and exact adapter from the
[documented recipe](assessment-2026-09-07-p8-constant-two-trace-family.md#verification).

One negative-control test was corrected before freezing: omitting the
curvature time boundary changes the old null rather than its rank by
itself; it breaks the shared null after the new mass term is included.
A source-expression factorization stall was diagnosed and resolved by
exact numerator expansion instead of unnecessary multivariable
factorization. Neither change modifies any frozen ancestor or weakens
the intended-action identities.

This is exact symbolic verification with written proofs and independent
internal calculations, not proof-assistant formalization or peer review.
