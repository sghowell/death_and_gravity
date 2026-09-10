# Exact anisotropic-box integral bound

Let h be the number of heavy edges, N=4+h, and B>=1.
The positive box integral is

    J(B)=integral_[light alpha<=1, heavy alpha<=1/B] U(alpha)^-2 d^N alpha.

Partition the box by all N! descending edge orders. In each
order choose the co-tree whose sorted ranks are lexicographically
first. Native code independently checks componentwise dominance
against every co-tree. Thus this monomial dominates all other
co-tree monomials in that ordered sector, and in particular

    U >= alpha_i alpha_j.

Write p_r=-2 at the two selected ranks and zero elsewhere.
For a suffix beginning at rank j define

    a_j=sum_(r>=j)(p_r+1).

All a_j are strictly positive for every selected graph and
sector. The proof therefore covers all shrinking parameter
flags, including disconnected edge sets.

The box has a further disjoint partition: its first k variables
exceed 1/B, and all trailing variables are below 1/B.
The first k must be light, so k ranges from zero through
the number of light edges before the first heavy edge.

The ordered small-variable integral is exactly

    B^(-a_(k+1)) / product_(j>k) a_j.

Independent nested integrations verify this formula for
every supported small power pattern. Multiplying by B^h
and changing each large variable to t_r=B alpha_r cancels
all powers of B, since sum_(r=1)^N(p_r+1)=N-4=h.

The remaining ordered integral is

    integral_(1<=t_k<=...<=t_1<=B) product t_r^p_r dt_r.

Set u_r=ln t_r and z_r=u_r-u_(r+1), with u_(k+1)=0.
The z_r are nonnegative, their sum is at most ln B, and
the exponent becomes sum A_j z_j with

    A_j=sum_(r<=j)(p_r+1).

Every A_j is nonpositive and at most two are zero. Both
facts are exhaustively checked. Integrate negative-exponent
coordinates to infinity; each supplies 1/(-A_j). If z
prefix exponents vanish, their simplex volume is
(ln B)^z/z!. This gives a positive upper bound, not an
equality after enlarging the integration domain.

Summing all pieces yields

    B^h J(B) <= C0+C1 ln B+C2 (ln B)^2.

Each coefficient is an explicit sum of positive rational
inverse products of tail and nonzero prefix exponents,
including the z! factor. All 258480 edge orders and 526320
cutoff pieces are checked. No floating optimizer is used.

The eight possible large-variable power patterns are also
integrated independently with their lower endpoint fixed
at one. The report records those exact formulas and the
logarithmic-coordinate identities behind their upper bounds.

The box is an auxiliary integration estimate, not a
physical momentum cutoff or a heavy-propagator expansion.
