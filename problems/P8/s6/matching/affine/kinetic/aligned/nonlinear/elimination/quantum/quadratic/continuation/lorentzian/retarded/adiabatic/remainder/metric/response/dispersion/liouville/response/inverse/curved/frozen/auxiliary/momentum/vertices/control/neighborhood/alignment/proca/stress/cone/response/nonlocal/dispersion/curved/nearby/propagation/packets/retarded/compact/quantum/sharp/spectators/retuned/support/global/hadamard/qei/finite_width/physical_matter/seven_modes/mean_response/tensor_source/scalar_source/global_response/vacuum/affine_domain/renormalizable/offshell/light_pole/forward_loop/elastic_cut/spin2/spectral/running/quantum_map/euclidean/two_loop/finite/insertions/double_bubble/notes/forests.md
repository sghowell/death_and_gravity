# A fixed reference and the actual recursive subtractions

Write I for the regulated light bubble including its
loop measure and I0bar=I0/(16pi^2) for the entire
S6.113 regulated mass-one zero-momentum reference.
All algebra is first performed at the same regulator;
one must not separately take the limit of a divergent
factor or discard evanescent terms before cancellation.

The logarithmic one-loop core operation evaluates that
core at zero external momentum and therefore replaces
I by I0bar. This is the inherited prescription, not
a silently substituted pole-only subtraction.

For the shared-vertex two-loop cores, explicitly extend
the reference by the same zero-momentum recursive
operation. This fixes this family's higher-order local
reference before any error estimate. The one-loop
scheme alone would not fix every higher-order finite
counterterm.

For a graph with one divergent bubble and one finite
triangle T, its two forests give

    I T - I0bar T = (I-I0bar) T.

For vertex-disjoint bubbles the four forests give

    I^2 - 2 I0bar I + I0bar^2 = (I-I0bar)^2.

For shared-vertex bubbles, the pair of both small cores
is not a forest. Instead, the six actual forest terms are

    I^2 - I0bar I - I0bar I
        - I0bar^2 + I0bar^2 + I0bar^2.

The last two terms are nested small-core/large-core
subtractions applied inside out. Equivalently,
R'=I^2-2I0bar I has zero-momentum value -I0bar^2;
its recursive overall counterterm is +I0bar^2.
Only after this operation does the same square result.

The code enumerates each selected graph's actual forests
and checks their individual terms and sum. It does not
insert the forbidden shared pair to obtain the desired
answer.

Locality is checked before attaching external heavy
propagators. Put C(z)=-L+g h(z), h(z)=1/(M-z).
The shared-core overall terms sum over the three channels to

    (I0bar^2/4)[-3L^3 + 2L^2g sum h - Lg^2 sum h^2].

These are literal tree variations from

    delta L = 3L^3 I0bar^2/4,
    delta g = L^2g I0bar^2/2,
    delta M = Lg I0bar^2/4.

They are this graph family's contributions, not all
two-loop parameter counterterms. Disjoint counterterm
products can also contain h^3, as two ordinary heavy
mass-counterterm insertions on a propagator must;
this is not a new nonlocal interaction.

The total canonical-field convention is retained.
Finite-potential counterterms at lower order inserted
inside other one-loop graphs are additional two-loop
contributions and are not claimed to vanish here.
A new two-loop momentum-independent potential contact
has zero second forward derivative, but that fact
does not remove those lower-order insertions.
