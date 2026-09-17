# Fixed-pair factorization and the conditional angular estimate

Fix i,j. Removing their cubic branch and its internal edge from every
full tree containing it leaves every allowed tree with the four Phi
leaves, the other N-2 physical h leaves and one amputated h root.
Conversely attaching the specified cubic branch reconstructs exactly
one original tree. The covariant mixed vertex is symmetric, and the
labeled source has no extra Bose factor at this amplitude stage.

The remaining graph sum is therefore the COMPLETE current R_ij.
Its only off-shell external leg is its root, because the other leaves
are the original physical on-shell waves. The proof in ward.md gives
Q_ij^T eta R_ij=0 at every finite N, away from other internal poles.

Apply the inherited exact cubic vertex and propagator:
A_ij=V_hhh(eps_i,eps_j,-trace_reverse(R_ij)/Q_ij^2)/sqrt(kappa).

S311 parametrizes any conserved symmetric R by its six spatial entries.
For each plus/cross pair, the coefficient after multiplication by
ab/(a+b)^2 is a polynomial P(z,r)/(1+r^2)^m with z=a/(a+b).
Every monomial is bounded in absolute value by1 for0<=z<=1,r>=0.
The exact absolute coefficient budgets are530,128,88,268. Unit-Frobenius
TT tensors have plus/cross coefficient l1 norm at most1, and spatial
rotations preserve the Frobenius norm. It follows pointwise that

|A_ij|<=530*(a+b)^2*||R_ij,spatial||_F/(sqrt(kappa)*a*b).

All24 coefficients and their denominator/degree checks are retained
from the frozen S311 input. The estimate is uniform in the EXPLICIT
pair angle for a specified current norm. R_ij itself depends on the
pair kinematics. No uniform norm on R_ij is proved here; its other
propagators can become singular in intersecting or simultaneous
collinear limits. The displayed inequality is not an integrated
all-N remainder estimate.

## Complete three-real calibration

There are no disjoint two-pair matchings on three labels. Thus
A3=A_no_pairs+A_01+A_02+A_12, and
5116=3814+3*434. The no-pair evaluator suppresses only pure two-leaf
h currents; it retains every other higher vertex and graph.

The exact full amplitude identity is checked at the rational physical
three-real state of S315, first at diagnostic couplings and separately
at the ORIGINAL parameters. The three pair terms are all nonzero.
Every complete hard remainder is also compared entry by entry with
the independent frozen S310 lower-multiplicity evaluator.

At N4, pair classes overlap in three disjoint-pair intersections.
Their remainders now have two off-shell external roots. The single-root
Ward argument cannot be applied independently to both; obstructions.md
is not a license to omit these intersections. The next quantitative
step must control the complete interacting remainder and its overlapping
singular regions, not recursively declare every descendant transverse.
