# Exhaustive quartic skeletons and full heavy refinements

The unchanged Gaussian-heavy reduced action is S2+S4 with a full
symmetric quartic kernel. Its three heavy channels are rational
propagators, not a derivative expansion. Connected four-point
graphs satisfy 4V=2I+4 and L=I-V+1. At L=2 this forces V=3, I=4.

Enumerate multiplicities of the six pairs 00,11,22,01,02,12.
Their sum is four, each vertex has degree at most four, the
graph is connected and deleting any non-self internal light
edge leaves it connected. The remaining half-edges are external.
This finite exhaustive search gives nine labelled skeletons:
three each of double bubble, wineglass and tadpole insertion.

For external counts e0,e1,e2, the number of assignments of four
distinct external labels is 4!/product(ej!). There are 72 over
the nine skeletons. With the 1/3! vertex expansion and normalized
quartic Feynman rule, the weight per fixed external assignment is

    1/[3! product_parallel(n!) product_self(2^l l!)].

The summed family weights are 3/4, 3 and 3/2, totaling 21/4.
These are positive combinatorial weights, not the signs or
values of renormalized amplitudes.

An independent zero-dimensional effective-action check uses
K(phi)=K+L phi^2/2. The two-loop background action is

    Gamma2(phi)=L/[8 K(phi)^2]-L^2 phi^2/[12 K(phi)^3].

The factors come from three pairings at the quartic vertex
and six fully connected three-line pairings between two cubic
vertices. The one-particle-reducible double-tadpole pairing is
excluded. Its fourth background derivative at zero equals
(21/4)L^3/K^4. This is an independent counting diagnostic,
not a substitute for continuum momentum integrals.

At each quartic vertex choose its local term (0), or one of
the three pairings (1,2,3) of its four distinctly labelled
half-edges. A heavy choice splits the vertex into two cubic
H Phi^2 vertices joined by one full heavy edge. Self-loop
half-edges retain separate endpoint labels during splitting.
There are 4^3=64 refinements per family, 192 in total. Every
refinement still has two loops and four external light legs.

The representative of each family suffices for the denominator
theorem because a vertex permutation relabels graph variables,
and an external permutation only permutes the three pair cuts.
The coarse bound treats all pair cuts with the same coefficient.
It therefore covers all 72 external assignments, not only the
particular forward labels used to display the refined graphs.

Refinement can produce a graph reducible through a heavy edge.
This does not contradict the original light-1PI criterion in
the exact reduced action. Heavy tadpoles and self-energy
subgraphs must not be discarded as an accidental consequence
of changing representations.
