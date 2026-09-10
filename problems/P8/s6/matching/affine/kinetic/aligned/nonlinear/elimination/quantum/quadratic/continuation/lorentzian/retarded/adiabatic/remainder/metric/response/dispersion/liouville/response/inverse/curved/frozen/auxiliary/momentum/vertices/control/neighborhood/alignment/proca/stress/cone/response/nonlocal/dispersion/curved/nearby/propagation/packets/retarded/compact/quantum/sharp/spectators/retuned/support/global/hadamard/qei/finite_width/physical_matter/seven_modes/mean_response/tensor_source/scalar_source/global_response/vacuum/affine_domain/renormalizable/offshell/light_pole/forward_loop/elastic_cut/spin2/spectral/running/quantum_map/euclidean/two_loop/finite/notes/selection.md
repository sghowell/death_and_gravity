# Exact finite-subsector selection and amplitude weights

Use every S6.120 heavy refinement whose full local graph
contains no superficially divergent loop subgraph. This
selects 88 of 192, with counts by heavy-edge number h:

    double bubble: h=1:2, h=2:16, h=3:22;
    wineglass:     h=1:4, h=2:20, h=3:24.

No tadpole-insertion refinement passes. All 104 other
refinements are explicitly rejected by the finite-subsector
API; their subtraction work is not claimed complete.

For each selected graph, deleting a pair of edges gives a
spanning tree precisely when its remaining N-2 edges connect
all vertices. Because N=V+1, this connected complement has
V-1 edges and is automatically acyclic. Its edge-pair product
is a monomial in U. This independent construction exactly
reproduces the frozen first graph polynomial in all 88 cases.

The strict small-scale exponents are also checked for every
ordered sector, not inferred solely from the 1PI inventory.
Thus possible disconnected shrinking edge sets are covered
by the direct integration argument.

For a refinement with h heavy edges, the absolute coupling
factor is lambda4^(3-h) g^h, where g=G^2. The full quartic
kernel has exactly its local term and three heavy exchanges;
there are no additional factors from expanding a derivative
series. The expansion uses the frozen S6.120 Wick weights.

The graph bound treats the three pair channels symmetrically,
so it applies to all external-label assignments and vertex
permutations represented by a family. Sum each representative's
choice bounds with family weights 3/4 for double bubbles and
3 for wineglasses. These already include the vertex factorial,
parallel-edge symmetry and all external-label assignments.
No further factor of 4!, 3! or an identical-particle phase-space
factor is inserted into this virtual graph calculation.

Weights in an absolute bound are positive. This does not
assert positivity of the signed finite-subsector correction.
