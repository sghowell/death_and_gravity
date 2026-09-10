# P8 S6.132 — complete one-loop fermion four-scalar increment

Original P8 remains OPEN. Keep the prospective GY14-unbroken
reference boundary and the explicit fermion local/field
ledger of S6.131. This is not an import of the old scalar
two-loop result into the new candidate.

## The selected calculation

Compute every one-loop fermion 1PI four-Phi box, the entire
local reference, the complete two-derivative term and a
bound on all higher momentum degrees. Combine this increment
with the reference tree amplitude only, using the S6.131
pole field normalization. Old scalar-loop graphs and the
conversion of their different reference scheme remain
outside this checkpoint.

There are six cyclic labelled boxes with both orientations,
N=6 active color/flavor copies, Y=y^2 and Q=16pi^2.
The continued Euclidean vertex Gamma_F4 contributes
A_F=-Gamma_F4 to the scattering amplitude.

## Complete low degrees and high-degree bound

At reference scale nu=mF,

    Gamma_F4(0)=24NY^2 Ibar/Q -64NY^2/Q.

The entire dimensional local reference is removed before
using four-dimensional integral norms. A second derivation
from the Dirac trace agrees with the frozen determinant.

The complete Lorentz- and S4-invariant quadratic-momentum
vertex is

    Gamma_F4,degree2 = [2NY^2/(QmF^2)] sum_i p_i,E^2.

It is fixed independently by a constant-background two-point
derivative. On the equal-mass-one shell its value is the
constant -8NY^2/(QmF^2). Thus degrees zero and two do not
contribute to the forward center second coefficient;
they are not discarded as off-shell operators.

On the complex forward disc |s-2|<=5, the full finite
remainder after these terms obeys, for mF>=36,

    |R_F(s)| < E_F = 3 times 10^8 Y^2/(QmF^4).

All higher momentum degrees are included. Exact permutation
identities make the convergent Euclidean series holomorphic
in the invariant s across the square-root coordinate
thresholds. Cauchy on the unit disc gives |b2_F|<E_F.

## Canonical dictionary and numerical scope

If kappa=1-f'(1), the selected normalized amplitude is

    (A0+A_F)/kappa^2.

If the potential quartic L+v4 is used in the displayed
tree, the residual box must instead be A_F+v4. The
same quartic threshold must not be counted twice.
The separate formal one-loop expansion is

    A0+h[A_F+2f'(1)A0],
    b2=4lambda+h[b2_F+8lambda f'(1)].

At the same numerical boundary E_F/(4lambda)<10^-604.
The pole-field contribution and its interval are explicit.
The selected normalized second coefficient is strictly
above 4lambda by less than a relative 10^-204; the
formal one-loop increment is separately positive.

These statements concern the selected tree plus one-loop
fermion functional. They are not an all-loop resummation,
the complete quantum model, a full scalar-loop scheme
conversion, a strict V/G test, a common bounce parent,
a whole-row exclusion or original P8 closure.
