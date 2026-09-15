# Full scaled geometry and cotangent reconstruction

## The invariant sublattice is not another truncation

The original six real nonzero Fourier inputs have wavevectors in
P Z^3. This is a closed Fourier subalgebra: products, derivatives,
projectors, and limits preserve its support. Equivalently its fields
are invariant under the three translations by 2pi/P.

The entire shape fixed-point equation is translation equivariant.
Its unique small solution inherits those translations. The full
mean-zero adjoint inverse, also unique and equivariant on its stated
complement, does too. Thus all generated homogeneous corrections and
arbitrarily high harmonics remain in P Z^3. Nothing is projected onto
the original six-mode support. The three residual mean translation
constraints are still retained; constants are not inverted.

Use y=P x and the scaled A2 norm with weight (1+|k/P|)^2. The initial
support has |k/P|=1, giving the exact factor 4/(1+P)^2 relative to the
old physical A2 row bounds. On the image ball r=8R this yields

    v_scaled <= 4r x 10^-478,
    tau_scaled <= 4r x 10^-423.

Physical first and second derivatives restore P and P^2. The full
physical A2 norms remain small enough for the original chart bounds.

## Complete curvature, including the quadratic terms

Keep Q=I+tau+Bf, with B(k) the full transverse projector and
trace B=2I+Pi0. In particular, the homogeneous eigenvalue is 3.
The complete determinant fixed point has norm bound
3(tau+f)^2+3(tau+f)^3. On the new ball f<=16 tau^2 it is a strict
self-map and its full derivative bound is less than 1/10.

Set h=exp(2v)Q^-1-I, so gamma=a^2(I+h). Full convolution and Neumann
bounds give

    ||h|| <= 8(v_scaled+tau_scaled),
    ||h-(2vI-tau)|| <= 100(v_scaled+tau_scaled)^2.

The second line includes the whole nonlinear determinant correction,
the full inverse series and all exponential cross terms.

Write e=||h|| in scaled A2. The inverse-metric gap is at most e/(1-e).
For each complete connection component use ||Gamma||<=9e; the linear
part is bounded by (9/2)e and the inverse-contact remainder by
(9/2)[e/(1-e)]e. Differentiating these products in A1 is covered by
the same A2 algebra estimates. The Ricci linear and remainder bounds
are respectively 6 Gamma_linear and
6 Gamma_remainder+18 Gamma^2. The final scalar contraction adds

    3 Ricci_remainder
    +9 inverse_metric_gap (Ricci_linear+Ricci_remainder).

This is below 10^4 e^2. The linear Ricci operator acting on the full
metric remainder costs at most 18 times its A2 norm. All these terms
are included in the reported bound.

The independent full 3D linear identity is
R1(h)=div div h-Delta tr h. Since tau is transverse and tracefree,
R1(2vI-tau)=-4 Delta v. Consequently the complete curvature obeys

    |R3| <= 4P^2 v_scaled
             +10^12 P^2(v_scaled+tau_scaled)^2 < 10^-327.

The factor a^-2 is at most one at real time. No nonlinear tensor
curvature is deleted; full exponential tensor fixtures have a
nonzero negative quadratic curvature. The proof uses the holomorphic
continuation of the REAL formal adjoint on complex phase inputs,
not conjugate transpose or asserted complex-metric positivity.

## Full cotangent lift, with the background retained

For every nonzero k in the sublattice the complete flat symbol and
inverse are

    M0=|k|^2 I+kk^T/3,
    M0^-1=|k|^-2(I-kk^T/(4|k|^2)).

The physical A0-to-A1 inverse norm is <=2/P. The full derivative-Q
adjoint contact remains in the Neumann correction, whose error is
bounded by 100||Q-I||physical_A2<1/10. Thus the full inverse norm
is at most 4/P. This is a full matrix inverse, not a scalarized
longitudinal or transverse projection.

The original metric and inverse physical A2 norms are less than 2.
The constant background a^2-I is now larger than the old tiny
bounce-centered metric gap, but is still small; it has no spatial
derivatives. The full DQ*, tensor projector and source bounds still
give

    base_A1=(100T+(1+P)Pi_v_fluct)/3
              +32(1+P)Pi_tau,
    D0 <= 180 base_A1
          +12(1+P)^2 W Pi_W
          +2(grad M1 + Pi_H grad H).

The retained matter momentum includes its nonzero 1/10 background
and is bounded by one. Every vector/Gauss and matter term remains.

The term Pi_v gamma^-1/6 is exactly pure trace pointwise. Its
tracefree projection vanishes, but its entire contribution to D0
above does not. With lambda_A1<=4D0/P the full tracefree momentum
therefore obeys

    pi_TF <= 5(32 Pi_tau +16 lambda_A1),
    shear <=48 pi_TF^2 <10^-370.

This is not deletion of a momentum mode or its constraint source.

Finally, compute all twelve invariant deviations from the actual
moving center. The trace and matter deviations include the complete
exp(-3v) density factors. All other electric, magnetic, mass,
gradient, heavy, Gauss and curvature invariants use the full bounds
above. Their maximum is below 10^-267.
