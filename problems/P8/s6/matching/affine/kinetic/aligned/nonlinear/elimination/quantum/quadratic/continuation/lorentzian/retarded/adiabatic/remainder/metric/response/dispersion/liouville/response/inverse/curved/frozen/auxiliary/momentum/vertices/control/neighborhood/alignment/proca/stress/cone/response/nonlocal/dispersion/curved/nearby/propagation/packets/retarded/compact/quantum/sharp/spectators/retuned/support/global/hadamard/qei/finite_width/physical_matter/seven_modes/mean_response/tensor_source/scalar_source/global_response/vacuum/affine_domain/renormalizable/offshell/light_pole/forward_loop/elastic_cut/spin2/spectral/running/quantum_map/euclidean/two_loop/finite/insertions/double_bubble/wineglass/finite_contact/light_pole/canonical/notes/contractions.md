# Independent full-label marked-core and insertion census

Start with all nine labelled two-loop light quartic skeletons from
S6.120, not only one family representative. Assign the four distinct
external labels to the actual external half-edges in every allowed way.
There are 72 assignments. At each of the three vertices select the local
term or one of the three exact heavy pairings. Thus 4608 labelled
full Phi/H graph records are generated, with 1152 double bubbles,
2304 wineglasses and 1152 insertion records.

The raw Wick weight per external assignment is
1/[3! product parallel multiplicity factorials product self-loop
2^n n!]. The heavy pairing expansion changes neither that weight
nor the light external labels. No heavy inverse is expanded.

A logarithmic one-loop proper core of type Phi4, H Phi2 or H2
consists exactly of two parallel light edges between two actual
local/cubic vertices. Selecting these cores independently is checked
against the full frozen UV-subgraph enumeration for all 192
representative refinements. The remaining light quadratic and heavy
one-point cores belong to the separately controlled inner OS grouping.

Contract each selected core to a marked local vertex of its boundary
type. Preserve all remaining external labels, vertices and typed
light/heavy edges, including multiplicities. Also preserve the original
coupling monomial L^nL G^nG. Canonical graph keys permit only internal
vertex relabellings within identical vertex-color/external-label classes.
Therefore a key equality preserves the full momentum-space cograph,
not merely a value obtained by setting external momenta equal.

There are 1152 marked contractions: 288 Phi4, 576 H Phi2 and
288 H2 cores. They give 144 distinct labelled cograph classes.

Independently generate the complete one-loop four-point graphs:
two quartic vertices joined by two light lines, six external-label
assignments, and 16 heavy refinements, for 96 records of weight 1/4.
At every local quartic vertex insert 3L^2/2 in place of L;
at every cubic vertex insert LG/2 in place of G;
at every heavy edge insert G^2/2 as an H2 mass vertex, splitting
the heavy propagator into two. One common factor
Ibar=I0/(16pi^2) is suppressed in these coefficients.

These are exactly the inherited one-loop local reference coefficients.
All 144 independent insertion keys and rational Wick weights match
the marked-contraction keys, with no missing or extra graph.

Signs are checked independently as well. A raw two-loop reduced
graph has sign (-1)^nL. A single subtraction adds a minus sign.
In the independent one-loop insertion construction, replacing an
L or G vertex keeps the original one-loop sign, while a heavy
mass insertion adds a minus sign from differentiating its inverse.
All signed coefficients match for each full labelled cograph.

Why this checks integrands: a proper light bubble at zero external
momentum is exactly the same regulated Ibar for every selected core.
Its vertex factors are momentum independent. Contraction leaves
precisely the marked cograph's momentum-dependent propagators and
external assignments. Equal graph keys and signed Wick sums therefore
give equal regulated counterterm integrands, up to harmless loop
momentum relabelling. The equality is before regulator removal.

The Wick-factor correspondence can also be understood by marking
a subset of contractions in the labelled expansion: collapsing those
contractions leaves the core and cograph contractions, and inserting
them reverses the operation. Dividing by the original vertex/edge
factorials gives the weights summed here. No unverified symmetry-factor
convention or zero-dimensional substitution is used to replace the
full-label comparison.
