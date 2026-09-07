# P8-S6.17.TREE: regular-flat HR-tree physical sign preservation

This gate enlarges the frozen star/bimetric background screens to a finite
acyclic tree of constant pairwise HR interactions, with separately coupled
and separately conserved NEC matter at any vertices. It is not a general
multigravity, ghost-freedom or UV theorem.

## Specified action and actual physical vertex

In four dimensions, orient each edge \(e:i\to j\) once and use

\[
 S=-\sum_v{G_v\over2}\int\sqrt{|g_v|}\,R_B[g_v]
 -2\sum_{e:i\to j}\int\sqrt{|g_i|}
       \sum_{n=0}^4\beta_{en}e_n\!\left(\sqrt{g_i^{-1}g_j}\right)
 +\sum_v S_{m,v}[g_v,\psi_v].
\]

The Einstein coefficients \(G_v\ge0\) and all real interaction
coefficients are finite constants. Choose the physical vertex \(r\)
with **\(G_r>0\)**. Zero-Einstein auxiliary vertices elsewhere are
allowed; the full equations at them must still hold. This positive target
condition is not inferred from a tensor Hessian or a propagating-mode map.
Reversing an edge replaces \(\beta_{en}\) by \(\beta_{e,4-n}\).
Endpoint terms are included. The normalization is the S6.16 action's
\(-2\beta\), not an unchanged copy of S6.5's \(-m^4\beta\).

Use P8(b)'s \(+---\), \(R_B=-6(DH+2H^2)\), EH \(-G R_B/2\), with
\(3G H^2=\rho\) and \(-2G DH=\rho+p\). Each actual metric is common
spatially flat homogeneous,
\(g_v=N_v^2dt^2-a_v^2d\mathbf x^2\), with positive finite smooth
lapse and scale on a connected regular interval, in the positive
square-root branch. All statements concern the actual physical metric
\(g_r\) and proper time \(dT=N_rdt\), not a composite or redefined frame.

Each matter sector couples to its own metric only and is separately
conserved on shell. Its isotropic null density is
\(n_v=\rho_v+p_v\ge0\); absent matter has \(n_v=0\). Distinct canonical
positive-field-metric scalars are sufficient. A single field coupled to
several metrics, or an exchange of matter energy between vertices, is
not covered by this assumption. Equal background scalar trajectories
do not make the fields the same field.

## The theorem

For every such solution and any two ordered physical times in the same
connected regular interval,

\[
 H_r(T_0)\le0\ \Longrightarrow\ H_r(T)\le0\quad(T>T_0).
\]

Strict negativity remains strict. Therefore an actual physical
contraction-to-expansion transition is impossible, including arbitrary
degenerate crossings and branch changes.

The proof first establishes a zero-divergence edge flux with reciprocal
weight \(N_v^2a_v^3\), then uses the injective incidence matrix of a
finite tree. At each time it selects the physical vertex's component of
edges with currently nonzero interaction polynomial. Internal dynamic
locks can be differentiated locally; boundary edges have zero null stress
at that point. The resulting fixed-subset kinetic coefficient is positive
because it contains \(G_r\). A locally bounded comparison coefficient
and the positive-part Gronwall argument cover all points. No global
derivative of a changing component membership is taken.

The theorem needs no sign of the beta coefficients, vacuum, subluminal
cone, auxiliary inverse, mass-gap condition, nonzero Hubble rate, finite
number of switches, finite polynomial-root decomposition or complete
proper-time tail. Identically zero interaction polynomials are harmless:
their edges are inactive, and the positive target Einstein term remains.

With an explicit actual-vertex-metric/proper-time dictionary to CD, the
necessary maximum Hubble error at \(T=\pm L\tau\) is
\(4L/[\tau(1+L^2)]\), hence \(8/(5\tau)\) at \(L=1/2\). This is an
actual full-parent background constraint, not a supplied map from an
approximate reduced action, matter/operator matching, or a cutoff estimate.

## Exclusions and evidence

Cycles, infinite graphs, nonpairwise or derivative interactions, variable
coefficients, nonflat/anisotropic or singular geometries, negative Einstein
coefficients, shared/nonconserved/non-NEC sources and a zero-Einstein target
without a separate positive-component guarantee are outside this theorem.
An explicit actual bouncing zero-target control explains the last guard.
A cyclic incidence circulation is only a failure of this proof's graph
step, not a constructed cyclic bouncing solution.

Exact symbolic action variation, independent Fraction/jet and Gaussian
incidence calculations, actual scalar backgrounds, source hashes and the
separately authored action audit support the [written proof](notes/proof.md).
No scalar/vector/tensor health, quantum, universal parent or P8 closure
claim follows. The completed scoped photon objective and frozen linear
classification are unchanged; original P8 remains open.
