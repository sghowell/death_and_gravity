# Entire canonical spatial generator and retained boundary

Use the canonical convention {q,p}=1. For symmetric metric variables the
unrestricted sum pi^ij delta gamma_ij counts off-diagonal entries twice;
equivalently their independent canonical dual momenta are 2pi^ij.

For a smooth spatial descriptor xi at fixed clock, begin with the whole
pairing

    D[xi] = integral pi^ij L_xi gamma_ij + piW^i L_xi W_i
            + pm L_xi M1 + ph L_xi H + pN L_xi N + pT L_xi T.

Ordinary Lie derivatives act on the covariant metric/vector and scalar
fields. Integrating by parts gives D[xi] = integral xi^i H_i plus the
boundary flux, with

    H_i = -2 gamma_ij D_k pi^jk
          + piW^j (partial_i W_j - partial_j W_i)
          - W_i partial_j piW^j
          + pm partial_i M1 + ph partial_i H
          + pN partial_i N + pT partial_i T,

    flux^i = 2 pi^ij gamma_jk xi^k + piW^i W_k xi^k.

Here pi^ij is a contravariant tensor density of weight one. In its
covariant divergence the connection trace from its tensor indices and
the density correction cancel. Expanding the displayed covariant
expression agrees with the direct integration by parts in all three
directions. The vector Gauss term is essential, not removable because a
longitudinal decomposition was used in an earlier quadratic calculation.

On a torus the integrated divergence vanishes. With other boundaries it
must be retained along with the prescribed boundary conditions. The
pointwise generator identity is false without that flux.

For configuration tensors/scalars q write D[xi]=integral p.L_xi q.
Functional derivatives give delta D/dp=L_xi q and the adjoint action on
p for delta D/dq. Hence

    {D[xi],D[eta]} = integral p.(L_xi L_eta - L_eta L_xi)q
                   = D[[xi,eta]],

up to the retained boundary. The representation identity is obtained
by expanding the tensor Lie derivatives; for a density its weight
contribution obeys div[xi,eta]=xi.grad div eta-eta.grad div xi.
Independent noncommuting three-dimensional polynomial descriptors
check all nine entries for both weights 2/3 and one. These diagnostics
support, rather than replace, the continuum functional argument.

The whole parent Hamiltonian at fixed u is a scalar density. R,F,j and
all profile coefficients are scalar functions of u,N, and its metric,
matter, vector and Gauss contractions have their displayed tensor
weights. Under a spatial diffeomorphism the auxiliary equations
therefore transform covariantly. Uniqueness of the parent's one local
regular auxiliary root makes N,T scalar functionals on that branch.
The parent's Dirac bracket on retained variables is canonical, so the
displayed spatial action and algebra survive that local reduction.

The primitive boundary F_boundary=integral V I is an integrated spatial
density. Its Lie variation is a boundary. Its exact canonical momentum
shift therefore preserves D[xi] under the same boundary prescription.
Do not discard the separate temporal endpoint terms or source pullback.

The field generator is at fixed time. General time-dependent spatial
descriptors also transform the shift multiplier and involve its primary
constraints. Neither those familiar multiplier terms nor the local
argument above establish a BRST quantum measure, an anomaly-free finite
regulator or a globally solved residual momentum constraint.
