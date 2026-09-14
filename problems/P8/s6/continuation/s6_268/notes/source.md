# Whole local source and its time-dependent reduction

This packet introduces a named finite IR/UV/phase-localized regulator.
It does not assert an identity with the original quantum theory outside
its local classical domain, nor a new action-level counterterm.

Fix a finite torus of side2pi L and a finite symmetric set Lambda of
nonzero Fourier wavevectors. Count both signs in |Lambda|. There are
|Lambda| real spatial basis functions, eight physical configuration
channels per real basis function, and d=8|Lambda| canonical pairs.
For one wavevector pair, d=16, not8. Use orthonormal real Fourier
functions, with the corresponding exact torus volume factors.

Keep the S267 shape and cotangent reconstruction in their full
infinite-convolution form. In particular, generated homogeneous shape
corrections and unrestricted nonlinear harmonics are not projected
away. The adjoint inverse in the momentum lift is the full elliptic
inverse, with the full metric, both matter channels, electric/Gauss,
magnetic, vector-mass and boundary terms. Constants are not inverted:
all three residual translation constraints remain. This is a finite
variational restriction of reduced canonical variables, not a
finite-mode invariant submanifold of the nonlinear PDE, nor a
collocation representation of its diffeomorphism algebra.

Homogeneous reduced variables are held at explicit external reference
values/functions. They are not assigned a new quantum state. Their
background values and time dependence must remain in every symbol.
This does not describe quantum homogeneous backreaction or replace
the original R3 state by an allegedly identical torus state.

At u=0 the finite reconstructed canonical origin maps to the interior
of the S266 invariant box: its actual auxiliary root is the whole
source-pinned root, not the bare root N=1. Choose a sufficiently small
closed finite phase ball so that all twelve reconstructed invariants
remain strictly inside delta/2 and the spatial shape remains strictly
inside its S267 ball. This existence follows from continuity at the
origin. The inverse ghost is elliptic on the mean-zero complement
with strict coercivity. Restricting to smooth finite-parameter input,
elliptic regularity and the full implicit-function construction give
smooth parameter dependence, including the unrestricted reconstructed
harmonics. No same-regularity diffeomorphism Banach-group theorem or
infinite-dimensional strong Darboux theorem is being assumed.

The full canonical Hamiltonian after the retained Maxwell integration
by parts is algebraic in N and the normal vector T. There are no
remaining spatial lapse derivatives in this auxiliary equation.
The exact time-dependent coefficients are the original R,F,B,Fhat,j,
fixed profiles, Hclock and constants in canonical.current_data.
The strict S266 temporal and lapse pivots, a compact smaller phase
ball, and the parameter implicit-function theorem therefore yield
one common positive time interval around u=0. All actual coefficients
and the boundary primitive are varied before solving for N,T.
The u=0 polynomial comparison is not extended off that slice.
Neither the canonical radius nor this time interval has been
numerically evaluated in this packet.

Restore the common kappa in both the original canonical one-form
and Hamiltonian before passing to unit CCR. The S251 scalar map
C=sqrt(kappa) diag(I,a^3 I) obeys C Omega C^T=kappa a^3 Omega;
tensor polarization factors, vector normalization, both scalar
boundary shears and the central canonical map are also retained.
One can work in the S267 canonical coordinates with constant
physical kappa normalization; a^3 in the linear state dictionary
then comes from the density momentum, not a deleted contact.

Here is the general moving-chart identity, including possible exact
boundaries. If the fixed-time pullback is

    j_u^*Theta|du=0 = P dQ+d_z beta_u,

then on extended phase space the pullback action is

    P dQ - h_red du + d beta_u,
    h_red = H_parent(j_u)-Theta_parent(partial_u j_u)+partial_u beta_u.

This is obtained simply by separating the spatial and time pieces
of d beta. In the chosen S267 cotangent chart beta=0; the metric
gamma=a_hat^2 exp(2v)Q^-1 gives the scale term -H_hat Pi_v.
A moving matter mean gives its full analogous momentum term.
If an additional momentum boundary shear is used, its beta_u
term must be included. The original primitive time/spatial boundary
is already in the parent and is not counted twice.

The independent nonlinear point-chart fixture uses
x=a(u)exp(q), p_x=(p+b(u))/x. Its one-form is

    p_x dx = p dq+d[b(u)q]+[(p+b)a'/a-b'q]du.

Thus its Hamiltonian connection is -(p+b)a'/a+b'q. Omitting the
last term changes the equations even though the fixed-time
symplectic form is correct.

Let h_ref be the complete quadratic generator of the original
free reference mode flow, S_0=I, and U_ref its metaplectic lift.
The SAME state psi0 is the original prepared pure nonzero-mode
state evaluated at u=0, not a new instantaneous minimizer.
Define the entire interaction-picture difference

    g(u,z)=(h_red-h_ref)(u,S_u z).

This includes all constant, linear, mixed-quadratic, H/light and
higher terms. It need not start at cubic order. The moving-chart
connection is in h_red once; the free-flow connection subtracts
h_ref once. It is invalid to identify the full Hamiltonian Hessian
with the selected free reference without checking that identity.

Choose an invariant smooth time-independent cutoff chi on a
smaller common interaction-picture phase ball. It equals1 on a
nonempty core, takes values in[0,1] and is compactly supported
inside the regular domain. Extend

    g_ext=g(u,0)+chi[g(u,z)-g(u,0)].

For F_red equal to the spatial average of the WHOLE normalized
physical volume exp(3v)R_full(u,Nstar)^(-3/4), first pull it back:
F(u,z)=F_red(u,S_u z). Then extend
F_ext=F(u,0)+chi[F(u,z)-F(u,0)].
Both backgrounds are retained scalar identities. F_ext is strictly
positive by convexity and compactness of the positive domain.
These definitions agree with the whole local symbols on the core
but explicitly change their continuation outside it.
