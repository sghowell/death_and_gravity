# Why all named through-two-loop coefficients are unchanged

The k-th new Yukawa vertex has 8k+1 scalar legs and two
fermion legs, k>=1. Its total valence is 8k+3 and its
weighted excess (valence-2)/2 is 4k+1/2.

For any connected perturbative graph with E source
endpoints, I internal edges and V vertices, the halfedge
and Euler identities give

    L_eff=1-E/2+sum_v[(n_v-2)/2+j_v],

where j_v is the counterterm loop order. In the physical
cubic-map pullback add its nonnegative map degree p.
All ordinary interacting nonvacuum vertices have
nonnegative weighted excess. In particular a linear H
reference counterterm first occurs at j>=1 and has
weight at least -1/2+1=1/2. A vacuum counterterm cannot
be attached as an extra disconnected vertex to a
connected correlator.

One new vertex therefore implies these integer lower
bounds on the first POSSIBLE changed graph:

| Total source endpoints | Loop lower bound |
|---|---|
| 0 | 6 |
| 1 | 5 |
| 2 | 5 |
| 4 | 4 |
| 8 | 2 |
| 10 | 1 |

The single-source bound need not be sharp. Nor is a
possible four-point graph necessarily a nonzero b2:
its external momentum structure still matters.
The bounds suffice to protect the named vacuum,
single-H, two-Phi and four-Phi coefficients through
loop two. They do not identify arbitrary higher-point
functions at that order.

Counterterm contraction cannot evade the counting.
For a subgraph with e external legs, Is edges, Vs
vertices and original counterterm orders summing to J,

    sum[(n-2)/2+j]=(2Is+e-2Vs)/2+J
                 =(e-2)/2+(Is-Vs+1+J).

The contracted local vertex has exactly this same
weight. The conclusion thus applies to complete
forests and finite counterterm insertions, not only
unsubtracted diagrams.

The cubic map has F(Psi)=Psi+O(Psi^3). Hence
f_R(F)-F begins at degree nine as well, and transformed
physical sources only increase the loop grade.
The proof concerns the regulated formal perturbative
coefficients; it does not bound a physical omitted tail.

Independent mixed-halfedge tests and Gaussian/Grassmann
diagnostics reproduce the different source-count orders.
They do not replace the connected graph proof.
