# Loop degree including physical sources and counterterms

For a connected graph, count E external physical source
endpoints as vertices. Let N be the sum of the action vertices'
parent valences, V their number, j their total counterterm loop
weight, and p the total number of cubic replacements in both
the action and the physical sources. Each replacement adds two
halfedges. Thus

    I=(N+E+2p)/2,
    L_eff=I-V-E+1+j
         =1-E/2+sum_action[(n_base-2)/2+j_vertex]+p.

Parent interactions have valence at least two. Consequently
p <= L_eff-1+E/2. For Phi four-point functions through two loops
p<=3, not two. For Phi two-point functions p<=2; for vacuum
graphs p<=1 at two loops. A local vacuum counterterm has p=0
and is kept separately. Linear H-source counterterms also have
p=0 and loop weight at least one, so their weighted excess
j_vertex-1/2 is nonnegative and does not weaken the bound.
An H one-point graph with map dependence needs an H-Phi-Phi
vertex or the corresponding generated vertex; through two
loops it has p<=1.

The pure even-boson two-loop four-point action patterns are
three quartics, one quartic plus one sextic, or one octic.
There is no justification for retaining sextics but dropping
octics. With counterterms, generated valence n and loop weight
j obey the necessary support condition n+2j<=8. The program
lists 42 supported operator/weight/map-degree tuples for the
parent Phi, H, fermion, gauge and ghost classes. This is a
necessary support list, not an assertion that every listed
tuple is a connected diagram or a new diagram multiplicity.

The identity counts source replacements as well as action
replacements. For example, a free scalar four-source connected
graph can reach map degree three at two loops. A count which
uses only action vertices misses this boundary. All generated
terms remain in the full action; the bound limits the formal
calculation, not the physical EFT or its remainder.
