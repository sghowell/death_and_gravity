# Local ultraviolet subgraphs and vertex-aware forests

The full-heavy refinement restores ordinary local quartic
Phi^4 and cubic H Phi^2 vertices. Enumerate every nonempty
internal-edge subset, keep its incident vertices, and test
connectedness and one-particle irreducibility. For E internal
propagators and L independent loops in four dimensions,
the superficial momentum degree is omega=4L-2E.

Retain loop subgraphs with omega>=0, including the whole graph
when appropriate. External light and heavy legs are counted
from the actual local valences, subtracting twice each internal
edge of the corresponding species. This detects subgraphs
with heavy external legs as well as purely light ones.

Over all 192 refinements the distinct
(light external, heavy external, omega, loops) tuples are

    (0,1,2,1), (0,2,0,1), (0,2,0,2),
    (2,0,0,1), (2,0,2,1), (2,1,0,1),
    (2,1,0,2), (4,0,0,1), (4,0,0,2).

These describe allowed local subtraction structures, not their
calculated coefficients. In particular a heavy tadpole is not
silently removed by a supposed heavy on-shell prescription.

A restricted forest is a set whose elements are pairwise nested
by edge inclusion or have disjoint vertex sets. Nested subtraction
operations must be applied inside out. Edge-disjoint subgraphs
sharing a vertex are not treated as independent here.

The bare double bubble has two one-loop bubble subgraphs and
the whole graph. The two bubbles share their middle vertex,
despite having disjoint edge sets. The six allowed forests are
empty, either single bubble, the whole graph, and the whole
graph with either bubble. The pair of bubbles is excluded.
The bare wineglass and tadpole insertion each have four forests.
The full refinements have one, two, four or six compatible forests.

This is an exhaustive combinatorial inventory. A renormalized
amplitude still requires a specified subtraction operator,
counterterm coefficients, the fixed finite scheme conversion,
and convergent bounds on the resulting boundary integrals.
The inventory alone is not the R operation evaluated on a graph.

Finally, a light-1PI four-point inventory is not the full two-loop
S matrix: the corresponding pole mass, residue and LSZ corrections,
including primitive two-loop self-energies, remain necessary.
