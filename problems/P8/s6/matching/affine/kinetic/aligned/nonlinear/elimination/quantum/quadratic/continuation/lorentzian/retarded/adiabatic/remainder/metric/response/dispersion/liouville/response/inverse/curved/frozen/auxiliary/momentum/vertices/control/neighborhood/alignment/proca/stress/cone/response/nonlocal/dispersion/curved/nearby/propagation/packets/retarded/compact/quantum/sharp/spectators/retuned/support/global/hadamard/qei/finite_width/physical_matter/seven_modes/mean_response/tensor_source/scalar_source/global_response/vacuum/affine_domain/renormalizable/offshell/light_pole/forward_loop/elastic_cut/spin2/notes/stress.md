# Literal stress insertion and complete spin-two projection

The minimally coupled matter action contains canonical light and heavy
kinetic terms and the same real polynomial potential. The two-scalar
canonical stress vertex for a line of mass squared mA is

T_mu_nu(k,k+q)=k_mu(k+q)_nu+(k+q)_mu k_nu
 -eta_mu_nu[k.(k+q)-mA].

Choose a null projection vector e with e^2=e.q=0. Its contraction is
2(e.k)^2, independent of the line mass. The code verifies the full
four-component tensor before projection, with an explicit null vector
and transverse transfer.

The on-shell matrix element decomposes as
2P_mu P_nu F(t)+(q_mu q_nu-eta_mu_nu t)D(t)/2. The same projection
isolates 2(e.P)^2 F(t). It removes metric potential/mass vertices and
curvature improvements of the form R times a scalar operator. Such
improvements supply only the transfer/metric tensor structure, even
when their scalar coefficient depends on t.

At one matter loop, the only surviving spin-two terms are the two
mixed triangles with the stress insertion on the light or on the
heavy internal line. Contact metric variations of the cubic/quartic
potential are pure metric tensors. A local quartic vertex attached
to a quadratic stress bubble supplies no individual external P to
the loop: after shifting and covariant integration its tensors are
only eta_mu_nu and q_mu q_nu, so its spin-two projection is zero.
The code checks the projected shifted numerators for both situations.

Tensor integrals are projected in a covariant regulator, such as
dimensional regularization, before its removal. The trace factor
e^2 is identically zero in that calculation; no undefined divergent
integral is multiplied by zero after the fact.

The light kinetic counterterm supplies a constant F contribution.
Its coefficient was already fixed by the on-shell residue and is
not chosen again here. Mass and potential counterterms do not supply
spin-two terms. Renormalizable curved-matter scalar-curvature
improvements do not change this projection. Independent finite
higher-dimension curvature operators would specify a different
matching action and are not silently assigned values or bounded.

External LSZ residue is one in the same scheme. Additional external
self-energy factors must not be counted a second time.
