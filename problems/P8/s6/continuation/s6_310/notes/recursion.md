# Complete labeled rooted-tree recursion

Fix external Phi leg0 as the amputated root. Each other external leg has
a distinct bit label and species Phi,H,h. A current with label set S
has momentum P_S equal to their all-outgoing sum. Singletons are the
matching external scalar value1 or graviton polarization; a species
mismatch contributes no graph.

Cut the vertex next to the root of a nontrivial current. Its subtrees
partition S into at least two nonempty unordered blocks. The recursion
enumerates each such set partition exactly once, assigns a species to
each block and uses the allowed full multilinear vertex. Its vertex
momenta are -P_S and the subtree P_B. The vertex is symmetric under
interchange of equal-species arguments, so no extra vertex or subtree
factorial is applied.

For a scalar current of mass-squared m2, divide the amputated source
by -(P_S^2-m2). For a graviton source U with upper indices, the lower
current is -(eta U eta-eta tr(eta U)/2)/P_S^2. Ten symmetric basis
evaluations reconstruct U; off-diagonal entries have a factor1/2 because
their basis matrix has two unit entries. The final external-Phi root
is not multiplied by a propagator.

Every recursive call uses proper label subsets. Cutting and reattaching
the root vertex are inverse operations on connected labeled trees.
Induction on the number of external labels therefore proves both
completeness and absence of duplicate topologies. A zero numerator is
still counted as a graph; the counting does not depend on a chosen
polarization accidentally annihilating an interaction. No singular
propagator value is evaluated.

An independent exponential generating series forgets all Lorentz and
momentum algebra. Let x count labeled external Phi, y labeled external
h, and P,H,h the three rooted-current series. A symmetric k-subtree
vertex contributes its usual1/k! in the EGF. With exp(h) representing
the symmetric metric attachments, use

P=x+P(exp(h)-1)+gPH exp(h)+C P^3 exp(h)/6,
H=H(exp(h)-1)+gP^2 exp(h)/2,
h=y+[(P^2+H^2)/2+gHP^2/2+C P^4/24]exp(h)+exp(h)-1-h.

Truncate only the formal label degrees to x<=3,y<=2; iterate until all
coefficients stabilize. The tree degree budget guarantees that formal
higher metric valences cannot change these retained coefficients.
Multiplying [x^3 y^N]P by3!N! yields

N=0: C+3g^2+3,
N=1: 5C+21g^2+21,
N=2: 38C+198g^2+198.

Thus the complete six-point inventory has434 graphs. These count
topologies, not contractions inside a vertex. Independent runtime
counts with the contact or heavy sector switched off are198,236,396.
The diagnostic switches never alter the original-parameter public
amplitude or any frozen scientific source.
