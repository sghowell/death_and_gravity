# Independent first-sheet denominator certificate

For each actual graph let alpha_e>0, with the first three edges light.
The spanning co-tree sum gives U. The sum over spanning two-forests
separating the two external labels gives P; it is zero if both labels
are at the same vertex. Isolated vertices are retained in the forest
enumeration. The number of edges and components enforces acyclicity.

Independently construct the weighted Laplacian with conductances
1/alpha_e, remove one ground vertex, and form the external unit
source/sink vector b. Literal determinant and adjugate identities give

    U = (product alpha) det Lap_red,
    P = (product alpha) b_red^T adj(Lap_red) b_red.

Both are checked against their independent graph sums for all 32 cases.
U has positive coefficients and P has nonnegative coefficients.
Let potential_num=(product alpha) adj(Lap_red)b_red, with ground zero,
and J_e=(potential_num_u-potential_num_v)/alpha_e.
These are polynomial current numerators. Native identities check
div J=bU, sum alpha J^2=UP.

When the external vertices differ, exhaustive path search supplies two
edge-disjoint paths. Send half a unit along each, obtaining f_e in
{0,+1/2,-1/2}, div f=b. If the external vertices coincide set f=J=0.
Native checks include both divergence identities and sum alpha fJ=P.
Consequently the exact polynomial identity is

    U^2 sum(alpha)/4 - U P
      = sum_e alpha_e (U f_e-J_e)^2
        + U^2 sum_{f_e=0} alpha_e/4 >= 0.

Thus P/U <= sum(alpha)/4; no unproved graph-resistance theorem
is needed to establish this finite certificate.

With actual light mass one and heavy mass-squared M>=1,
F=U(sum_light alpha+M sum_heavy alpha)-sP.
For |s-1|<=2, Re s<=3, so

    Re F >= U sum(alpha)/4 + (M-1)U sum_heavy(alpha).

Multiplying the gap at s=3 by U gives exactly three times the
displayed nonnegative square polynomial. This identity is checked
independently for every graph. A small neighborhood of the closed
disc keeps a strictly positive coefficient of U sum(alpha).

This is a common first-sheet interior denominator bound, not a
UV-boundary integrability theorem. The local sunset and the proper
subtraction-dependent graphs need the separate operations proved
in the following notes.
