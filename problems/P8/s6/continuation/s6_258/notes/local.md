# Local flat-torus slice and whole weak-coercivity estimate

Work first on a spatial torus with period 2 pi L in every direction,
L>0 finite. This is an auxiliary local classical gauge domain, not a
replacement of the original infinite-volume quantum preparation.
A periodic divergence has zero mean.

## Local slice

Use periodic near-identity coordinate maps phi(x)=x+xi(x), fix the
mean displacement of xi to zero, and require det Dphi>0. Let the metric
be near a homogeneous positive multiple of I. For an integer k>=2 and
0<alpha<1 use metric coefficients in C^(k+3,alpha), coordinate maps in
C^(k+2,alpha) and gauge residual in mean-zero C^(k,alpha).
The positive determinant, inverse, products and pullback then define
a C1 gauge map. The extra metric derivative controls differentiation
of composition and continuity of its highest coefficient derivatives;
no differentiability on the same-regularity diffeomorphism group is
being presumed.

At the reference the derivative in xi is M0=-Delta I-grad div/3.
Its Fourier inverse from the gauge note is an isomorphism from
mean-zero C^(k,alpha) to mean-zero C^(k+2,alpha), by the constant
coefficient elliptic estimate on the fixed torus. The Banach implicit
function theorem therefore gives one small periodic gauge map for
sufficiently nearby metrics, unique in this local mean-displacement
chart. It also remains orientation preserving when the displacement
gradient is small. The argument does not evaluate a global chart
radius, solve global momentum constraints or fix translations as a
group at all orders.

Translations of an exact gauge representative remain gauge
representatives. The mean-zero complement used by the theorem is not
a Lie subalgebra: on a coordinate circle the vectors sin x partial_x
and cos x partial_x have zero mean but bracket -partial_x. Thus it is
incorrect to discard translations from a purported closed BRST
algebra. Their residual-group treatment remains separate.

## Whole operator estimate on the exact slice

Assume now that chi=div Q=0 exactly and ||Q-I||Linf,op <= delta.
Integrating the entire on-slice operator gives

    B(v,xi) = integral partial_j v_i Q^jk partial_k xi_i
                + (1/3)partial_j v_i Q^ij div xi.

The two integrations by parts create no unaccounted derivatives of Q:
symmetry and div Q=0 remove them. At Q=I,

    B(xi,xi) = ||grad xi||2^2 + ||div xi||2^2/3.

For E=Q-I, the first perturbation is bounded in absolute value by
delta ||grad xi||2^2. For the second, use ||E||F <= sqrt(3)delta and
|div xi| <= sqrt(3)|grad xi|. It too is bounded by
delta ||grad xi||2^2. The diagonal Cauchy identity includes all
off-diagonal gradient entries, so

    B(xi,xi) >= (1-2delta)||grad xi||2^2.

At delta<=1/8 this is at least 3||grad xi||2^2/4.
The full form is bounded on mean-zero H1 fields. The real nonsymmetric
Lax-Milgram theorem gives a unique weak solution on that complement.
For sufficiently smooth Q, the full strongly elliptic system and
elliptic regularity give the corresponding classical regularity.

Poincare on this torus is ||xi||2 <= L||grad xi||2. For Mxi=f this gives

    ||Mxi||2 >= 3||xi||2/(4L^2),
    ||M^-1||L2_to_L2 <= 4L^2/3.

This is a singular-value/coercivity bound, not a claim that the entire
operator is self-adjoint. The gap is explicitly not uniform as L tends
to infinity. The same weak estimate holds on a mean-zero Galerkin
subspace, because it is a restriction of the exact form, but this fact
alone supplies neither a closed discrete Lie algebra nor a quantum
regulator. The exact antisymmetric three-cell derivative in the packet
has summation by parts and a nonzero Leibniz defect.

Every homogeneous comparison slice of the unchanged current parent has
Q=I regardless of its lapse or homogeneous scale. That connects the
local reference to the evaluated classical interval without identifying
it with the fixed interacting quantum mean or proving a physical
nonlinear Cauchy theorem.
