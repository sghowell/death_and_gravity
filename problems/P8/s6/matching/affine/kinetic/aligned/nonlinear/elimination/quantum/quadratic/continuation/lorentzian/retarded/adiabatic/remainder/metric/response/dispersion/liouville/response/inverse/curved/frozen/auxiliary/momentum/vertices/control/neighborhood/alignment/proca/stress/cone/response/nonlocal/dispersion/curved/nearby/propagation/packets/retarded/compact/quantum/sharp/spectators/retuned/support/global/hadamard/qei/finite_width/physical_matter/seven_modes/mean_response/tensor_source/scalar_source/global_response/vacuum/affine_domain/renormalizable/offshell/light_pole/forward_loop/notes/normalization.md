# Full vertex and one-loop normalization

Integrate H exactly at a translation-invariant Euclidean regulator.
The light action is S_eff=Phi K_m Phi/2 + lambda4 sum(Phi_i^4)/24
-g J K_H^(-1) J/8, with J_i=Phi_i^2. Its field-independent heavy
determinant does not exhaust the quantum effective action.

The arbitrary-index quartic tensor, differentiating this actual action, is

Q_ijkl = lambda4 delta_ij delta_ik delta_il
 -g[delta_ij delta_kl (K_H^-1)_ik
    +delta_ik delta_jl (K_H^-1)_ij
    +delta_il delta_jk (K_H^-1)_ij].

Each term is obtained by distributing two derivatives to each J. Symmetry
of K_H^-1 combines the two assignments. Fourier transformation of the
three pairings gives the complete on-shell vertex
A0=-lambda4+g[h(s)+h(t)+h(u)], h(z)=1/(M-z).
No momentum channel is replaced by its local expansion.

Writing S_eff''=K_m+W, W is homogeneous quadratic in Phi. The quartic
one-loop term is -Tr(G_m W G_m W)/4. Its fourth derivative distributes
two external legs to each W in six ways. Cyclic trace symmetry pairs
the two assignments for each of the three channels, giving symmetry
factor 1/2 per bubble. The independent generic two-site calculation
checks all 16 entries of Q and all five independent fourth derivatives
of the actual trace. This finite check supplements the general index
derivation; it is not a claim that a two-site lattice alone proves the
continuum theory.

Consequently the full one-loop light four-point function is the sum of
three bubbles with **complete** nonlocal quartic vertices at both ends.
Expanding their products would yield the original pure-light bubble,
mixed triangles, boxes and heavy self-energy insertion; all are retained.
There is no six-light-field vertex in the exact Gaussian effective action,
so no additional one-loop six-vertex tadpole is missing.

External self-energy corrections are removed by the independently fixed
on-shell light mass and unit residue. One-loop mass/kinetic counterterm
insertions into one-loop internal propagators would instead be two-loop.
Heavy one-point counterterms contribute only local light quadratic and
vacuum terms after the Gaussian integration.

For each channel, combining the two light propagators with parameter x
and continuing as specified in notes/analytic.md gives radial measure
y dy/(32pi^2), including the bubble symmetry factor. The s-channel
vertices coincide in the forward configuration. The t-channel vertices
generally differ; the estimates below bound their product, not an
incorrect complex-conjugate product.
