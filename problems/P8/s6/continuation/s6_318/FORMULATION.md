# Uniform weighted pure-soft currents at all finite multiplicities

## Unchanged source and physical domain

Use the same selected four-dimensional covariant action, canonical
g=eta+2h/sqrt(kappa), original parameters and gauge-complete pure-soft
temporal recursion as S317. In the original physical domain kappa=1e800,
mu=nu=1 and the total radiated positive energy is at most1/8.
No interaction, matching coefficient or admissibility assumption changes.

Let S contain n physical future null TT leaves with energies w_i>0.
The leaf polarizations have unit Euclidean-Frobenius norm. Work on the
generic finite-tree domain, away from exact internal propagator poles.
Write Q_S=(W,Q_sp), W=sum_i w_i, v=Q_sp/W and delta^2=1-|v|^2.

For the temporal current H_S with spatial block A, define

N(H_S,Q_S)=max(||A||F,||A*v||/delta,|v^T*A*v|/delta^2).

At a TT singleton delta=0 and both numerators vanish, so its seed norm
is the leaf Frobenius norm. Generic larger subsets have delta>0.

## The bound

At every finite n,

N(H_S,Q_S)<=C_n*W^n/[product_i(w_i)*kappa^((n-1)/2)],

C_1=1,
C_n=64*sum_(labeled root partitions pi with k>=2 blocks)
       (k+1)!*32^(k+1)*13^k*product_(A in pi)C_|A|.

The coefficients have the explicit all-n envelope

C_n<=2*(2*10^10)^(n-1)*n!.

The weighted norm includes more than the field norm: it retains the
small longitudinal components needed to close the next recursive step.
Its angular bound covers arbitrary three-dimensional directions and
nested collinear hierarchies, not just a planar or fixed-energy family.
The W^n/product(w_i) factor explicitly retains arbitrary soft-energy
hierarchies. Unnormalized leaf norms multiply the right-hand side.

## Proof boundary

Future-ray geometry transfers child weighted norms to the parent frame.
An arbitrary transverse Rosen metric makes every constant transverse
Einstein vertex vanish. A proper transverse rotation fixes the remaining
even/odd grading. The frozen all-order vertex majorant extends to exact
polynomial coefficient-l1 norms. Complete-source conservation and the
temporal propagator then cancel the dangerous root powers.

Finite labeled induction proves the current estimate. A nonnegative
auxiliary counting-series barrier proves the factorial envelope; it
does not claim convergence of the physical all-N perturbation series.

No value is assigned to an excluded exactly-collinear propagator point.
The full hard four-scalar/N-graviton amplitude, overlapping all-N soft
subtraction, inclusive probability and real-virtual completion remain
separate. Hard/evanescent matching, the interacting state, unitarity,
absolute complex Regge and the common-parent bounce are not established.
Original V/G/B/P8 remain OPEN; scoped P8(a) is unchanged.
