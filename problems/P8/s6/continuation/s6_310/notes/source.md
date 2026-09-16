# Original source and scalar-forest ownership

S309 is the frozen immediate parent. Copy its source checks without
mutating any cached parent packet. S297 supplies the literal full source
jet audit; its six lower R-jet residuals are replayed here, not replaced
by the distinct clock-background germ. S304 supplies the canonical
Einstein normalization and the complete one-graviton tree.

At the inherited formal tree order the relevant covariant terms are
(kappa/2)sqrt(-g)R, the minimally coupled Phi and Gaussian H kinetic/mass
densities, -g sqrt(-g) H Phi^2/2 and C sqrt(-g) Phi^4/24.
The literal fourth-degree scalar derivative terms cancel. The vacuum R
coefficient has no scalar jet through weighted degree4. Its first r_Y
has four scalar fields, so the Ia contributions begin at six Phi; the
dependent a4/a5 differences begin at twelve. The higher heavy source
begins with H Phi^2 Y, containing four Phi, while the regular linear
Proca source contains at least four Phi. There is no H^3 or H^2 Phi^2
coupling. These are source facts already checked in the original chain,
not assumptions about a different minimally coupled model.

Here is a tree-specific exclusion argument that is stronger than total
valence alone. Remove all non-Phi edges and vertices with no Phi half-
edges from a connected full tree. The remaining Phi forest has C_phi
components. If d_phi(v) is the number of Phi half-edges at an interaction,
its external-Phi count is

E_phi = 2 C_phi + sum_v(d_phi(v)-2).

This follows by summing degrees in each finite tree component, including
its external endpoints. All relevant interactions have even d_phi>=2.
Thus any vertex with at least six Phi already forces E_phi>=6.
A higher linear H source with four Phi must terminate its H path at
another source: H^2 metric vertices do not terminate that path. Its two
ends cannot lie in the same Phi component, since that would make a cycle
in the full graph. Therefore C_phi>=2 and E_phi>=4+2=6.
Two four-Phi linear Proca sources similarly require E_phi>=8.
Quadratic-only spectators without external spectator legs require a
closed species line, and hence a loop; the same applies to ghost lines.

For E_phi=4 the possibilities reduce to one Phi^4 vertex, one H exchange
joining two Phi components, or minimal Phi/Einstein interactions joining
them gravitationally. Multiple H exchanges require either a cycle or
at least three Phi components. Attaching metric trees does not reduce
d_phi or the number of external Phi. No discarded higher source vertex
can re-enter merely because there are two external gravitons.

For six total external legs, sum_v(degree(v)-2)=6-2=4. Hence the finite
vertex inventory needs scalar-pair metric jets through h^3, the H Phi^2
and Phi^4 density jets through h^2, and Einstein h^3/h^4. Although the
implementation permits the H^2 h^3 metric jet as well, graph topology
prevents it from adding an extra four-Phi/two-h graph.
This is a formal tree source-ownership proof, not an all-loop
decoupling or quantization-limit theorem.
