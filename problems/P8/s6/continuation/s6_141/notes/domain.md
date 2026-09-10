# Joint complex domain with an unshifted chord

Rotate a word to start at one boson endpoint.
Assign real Euclidean momenta q and l to the
two arcs at that endpoint. The boson momentum
is q-l. After the first arc's two external
momenta sum to P_arc, the second endpoint
carries l+P_arc and q+P_arc. Its difference
is still q-l. Each fermion propagator therefore
has q or l plus a prefix sum of external legs;
the boson is not shifted by complex external
momenta. No change of integration contour is used.

For real k and Hermitian Euclidean Clifford
matrices, D(k)=m+i slash k satisfies
D(k)^dagger D(k)=(m^2+k^2)1. Therefore
||D(k)^-1||=1/sqrt(m^2+k^2). For a complex
shift delta, ||slash delta||<=||delta||_1.
The Neumann identity gives a shifted inverse
norm at most 20/[19 sqrt(m^2+k^2)] when
||delta||_1<=sqrt(m^2+k^2)/20. Use the larger
bound 2/sqrt(m^2+k^2).

On |s-2|<=1, let E^2=s/4 and P^2=s/4-1.
Choose the four incoming Euclidean momenta
(iE,P), (iE,-P), (-iE,-P), (-iE,P).
Their squares are -1, they sum to zero, and
the invariants are s,0,4-s. Analytic square-root
branches exist on a neighborhood of this disc.
Both |E| and |P| are below one. Every leg has
L1 norm below two. The two-leg arc prefixes
have norm below four, and hence below 18.

Scale all external momenta by complex zeta.
For S_q=m^2+q^2 and S_l=m^2+l^2 choose

    R(q,l)=min(sqrt(S_q),sqrt(S_l))/360.

All fermion shifts on |zeta|<=R are strictly
below the appropriate sqrt(S)/20. This proves
holomorphy with a strict margin for every pair
of real loop momenta. For m>=720, R>=2.
This is a joint two-loop argument; no subgraph
radius is extrapolated past a transfer threshold.

The six inverse norms and |Tr A|<=4||A||
give 256/(S_q^(3/2) S_l^(3/2)).
For the gauge chord each of the four Euclidean
gamma indices has endpoint norm one; taking
absolute values of the sum costs at most four.
The exact color contraction is Cf times the
identity. Scalar and gauge bounds therefore
add with coupling N Y^2 (Y+4aCf).

The chord is the positive real denominator
(q-l)^2+b, with b=1 for scalar exchange and
b=0 for gauge exchange. The massless denominator
is an upper bound for both. Its diagonal
singularity is integrable in four relative
dimensions; the exact angular/radial majorant
below also controls the common zero-momentum
endpoint. A pointwise infinity at q=l has
measure zero and is handled by that integral.
