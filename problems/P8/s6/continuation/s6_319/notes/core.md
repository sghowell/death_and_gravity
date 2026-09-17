# Unique hard cores and disjoint energy charges

## Unique maximal-soft contraction and the three hard cores

In any full tree, remove each maximal connected rooted branch which
contains only real graviton leaves on its outer side. Maximal branches
are disjoint and cover every real graviton leaf. Replace each by its
complete temporal current, including all its internal root partitions.
The partition into maximal branches and remaining core is unique.
Summing the internal branch trees first is essential: the S318 bound
belongs to their complete sum, not to one branch graph.

To classify the remaining graph, retain only Phi edges. This is a
forest, with four labeled external degree-one Phi vertices, internal
degree-two vertices, and degree-four Phi4 vertices. The forest identity
4=2*c+2*v4 implies either v4=1,c=1 or v4=0,c=2.
Here c counts nonempty Phi components. No Phi cycle is allowed.

Every H component is a path: H has degree2 at H2 metric vertices and
degree1 at H Phi2 vertices, and there are no external H legs.
Its two endpoints therefore lie on the Phi forest. A path joining
two already connected Phi points would create a cycle. If the Phi
forest is connected, no H path is possible. If it has two components,
at most one H path is possible and it joins those two components.

For the connected Phi4 case the remaining matter core is connected,
so any h connection with two matter attachments would also create a
cycle. All h branches off it are among the contracted pure-soft blocks.
This is the contact core, with four distinguished external Phi chains.

When an H path joins the two Phi components, the same argument forbids
an additional h connection between points of the connected matter core.
This is a single H path, with its two H Phi2 endpoints, and four
external Phi chains. The four Phi labels give three pairings.

When no H path exists, the two Phi components are joined by one h path.
Its endpoints are Phi2/h vertices. A branching vertex on that path has
two hard path legs; every remaining branch is a contracted pure-soft
block. This is the hard-graviton core, again with three pairings and
four external Phi chains. No extra fourth class can occur.

Degree-two vertices on the external Phi chains have at least one soft
block. Extra H-path metric vertices and extra hard Einstein vertices
likewise have at least one soft block. Endpoint and contact vertices
may have zero blocks. All their arbitrarily high metric valences are
retained. The endpoint distinction orders a path; identical attached
blocks are unordered labeled sets, not additional ordered permutations.

As a count-only check let x mark one contracted block and z=exp(x).
External Phi chains have factor1/(2-z), a potential endpoint/contact
has factor z, and either internal hard path has factor1/(2-z).
The resulting block-core EGF is
 C*z/(2-z)^4+3*(1+g^2)*z^2/(2-z)^5.
Composing x with the pure-soft graph EGF gives precisely S315's full
inventory. Independent literal root cuts with all non-singleton
pure-soft currents disabled check the uncomposed core counts.

## One inverse soft energy per block is enough for a coarse bound

Each light-scalar propagator on an external chain is paired with its
outer emission vertex. The propagator carries that external scalar
plus the nonempty cumulative radiation beyond the vertex, up to
orientation/sign. The S312 bound gives8/(3*w_cumulative).
Choose one directly attached soft block at that vertex (for example,
the one with smallest leaf label). Its energy w_block is no greater
than w_cumulative, so8/(3*w_block) is a valid upper bound.

Different emission vertices own disjoint directly attached blocks.
Therefore no block is charged twice. Enlarge the inverse-energy
product by adding1/w_block for every uncharged maximal block; this
increases the bound because0<w_block<=1/8<1 in mu=1 units.
Keep8/3 once per scalar propagator/vertex, not once per block.

For a block containing m original leaves, S318 gives
 ||H_block||F <= C_m*w_block^m/
 [product_i(wi)*K^((m-1)/2)].
The attachment contributes the additional canonical K^(-1/2).
After the added energy inverse, this is bounded by
 C_m*w_block^(m-1)/[product_i(wi)*K^(m/2)]
 <= C_m/[product_i(wi)*K^(m/2)].

Thus the complete product of charged and uncharged blocks leaves
only K^(-N/2)/product_i(wi). This enlarges an inequality; it does not
claim a factorization of the physical poles. The sharper ordered
eikonal identity is not needed here and remains important for future
leading-soft subtraction. No new inverse child invariant is introduced.
