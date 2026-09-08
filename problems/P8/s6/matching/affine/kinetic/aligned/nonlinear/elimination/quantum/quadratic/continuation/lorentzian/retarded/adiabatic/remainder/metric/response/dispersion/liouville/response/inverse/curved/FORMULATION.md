# S6.73: prepared homogeneous tree-plus-Gaussian causal inverse

Keep the original physical clock on I=[-1/2,1/2], the S6.56
strict classical margin, the S6.55 selected all-order vector
state, the S6.67 finite prescription and the S6.60 fixed scalar
tadpole. The free classical matter field has zero prepared
charge perturbation. Independent lapse/log-scale sources and
forcing vanish on an initial neighborhood. No source-dependent
state or tadpole-profile reselection is permitted.
Fix m>=1000, L>0 and the classical margin parameter in (0,1/100].

Let f=(n,zeta), and let Kphys be the physical linear Gaussian
vector plus fixed-tadpole force response in order
(-delta rho,3 delta p), in the dimensionless gravitational
normalization used by S6.68. Define Q=64*pi^2*L^2*Kphys and
gamma=1/(64*pi^2*L^2)>0. Write E for the two classical metric
Euler equations after the zero-charge matter reduction and
the regular physical point-chart pullback.

## Normal form

Let I4 be four zero-past time primitives. With
S(u)=diag(1/h(u),1), h=(1+u^2)^3, and B_m the exact h=1
S6.72 massive block, define P=S B_m S, with the inner
multiplication inside the causal convolution. Then, on the
smooth prepared domain,

I4 Q = P + V.

V is a causal integral operator whose kernel, and every fixed
number of derivatives along its diagonal, are bounded by
C_j*(1+abs(log(u-v))) on the compact causal triangle. Each C_j
is finite; these constants are not numerically evaluated here.
Consequently V has a same-norm C0 extension, with short-interval
norm at most C_0*T*(2-log(T)) for 0<T<=1.

The proof uses the actual physical Hamiltonian current kernel,
its selected-state WKB expansion, the common dimensional
subtraction and the fixed fourth-order finite local matrix.
It does not replace the curved state by a flat vacuum.
The matrix P is a reference operator whose difference is now
controlled, unlike the definition-only preconditioner in S6.72.

## Coupled linear response

For every smooth prepared two-force g, the equations

E f + gamma Q f = g

have a unique smooth zero-past solution on I. The free matter
perturbation is recovered by s'=-3*ell*v-w*n, with
v=zeta-delta*n, delta=1/(2h). The inverse is causal and has a
finite C0 estimate against I4 g (hence also against g). No
numerical inverse constant or perturbative smallness estimate
is claimed. A weighted Volterra norm proves existence on the
whole fixed compact interval even when gamma^-1 is large.

This is a statement about the full retained *homogeneous,
prepared, linear* tree-plus-Gaussian-vector-plus-fixed-tadpole
response. It is not a claim about other quantum loops, arbitrary
initial covariances or freely prescribed higher-derivative
initial data. The forward logarithmic operator is not asserted
to be a C0 or C1 endomorphism. General continuous forcing is
interpreted in the zero-past distributional extension; smooth
prepared forcing gives the original smooth equations.

## Boundary

The proof supplies no nonlinear neighborhood of exact
semiclassical solutions, no quantum stability or cone spectrum,
no interaction/cutoff control and no V/G/B matching. A causal
finite-time inverse may grow strongly; it need not be a healthy
physical all-frequency propagator of the finite EFT.
Original P8 remains open, without a user-intervention blocker.

Evidence consists of independent exact vertex, dimensional,
time-primitive and literal classical-Euler identities together
with the written UV-extension and Volterra proofs. It is not
proof-assistant formalization.
