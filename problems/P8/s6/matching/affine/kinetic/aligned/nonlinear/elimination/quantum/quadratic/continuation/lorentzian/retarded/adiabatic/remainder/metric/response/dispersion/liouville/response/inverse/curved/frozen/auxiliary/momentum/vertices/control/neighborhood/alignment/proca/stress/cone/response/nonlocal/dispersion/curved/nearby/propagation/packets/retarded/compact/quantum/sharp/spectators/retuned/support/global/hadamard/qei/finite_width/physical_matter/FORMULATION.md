# S6.102: actual quadratic physical matter density and pressure

This checkpoint reconstructs the actual minimally coupled scalar matter
observable of the separately named epsilon=1/200 constant-Proca action.
It is not S6.101's positive test energy. Original P8 remains OPEN.

The source signature is -+++. All formulas use the inherited normalized
M=tau=1 units. The overall action factor is kappa=M^2 tau^2; physical
stress is M^2/tau^2 times the displayed normalized stress. This kappa
is not the affine connection deviation that older sources also denote
by that letter. The quantum covariance normalization remains hbar/kappa.

For the actual clock foliation let N be the hatted lapse, h=(1+u^2)^3,
e=[(h-1+N^-2)/h]^(-1/4), and U=e^3. Let
dc=P_psi/hat_volume-ell, ell=1/[10(1+u^2)^6], and
G=hat_g^ij partial_i chi partial_j chi. The physical unit-normal
density and one-third spatial trace are exactly

    rho=(ell+dc)^2/(2e^6)+G/(2e^2),
    p=(ell+dc)^2/(2e^6)-G/(6e^2).

They are scalar frame contractions, not coordinate T_00. The shift
has already canceled through the canonical normal velocity.

Write N=1+n1+n2+O(phase^3), dc=dc1+dc2+O(phase^3), G=G2+O(phase^3).
The common quadratic kinetic contribution is

    K2=ell dc2+dc1^2/2-3ell dc1 n1/h
       -3ell^2 n2/(2h)+3ell^2(1+3h)n1^2/(4h^2),

and rho2=K2+G2/2, p2=K2-G2/6. Thus quadratic physical stress needs n2,
although the stationary quadratic Hamiltonian needs only n1.

The complete lapse force and spatial momentum constraints are
reconstructed from the actual action, with no division by Theta.
At the bounce u=0 all 105 independent pairs of the fourteen physical
phase channels give exact density, pressure and omitted-n2 Hessians.
Three nonzero-output mixed pairs are also reconstructed at each of
u=-1/100,0,1/100. The general invariant formula is exact at quadratic
order; these fixtures are controls, not a sampled proof for all angles.

An independent closed center calculation reproduces all 588 entries
of the three 14 by 14 matrices. The scalar map to the existing quantum
packet retains BOTH canonical boundary shifts and its inverse-k
symplectic normalization. The actual principal density Hessian is
indefinite, including a scalar configuration value -134/135 and
negative tensor/Proca momentum entries. This obstructs an automatic
positive-square transfer of the previous test-energy inequality.

It does NOT imply negative full classical matter energy, a ghost,
an instability, an arbitrary-amplitude perturbative regime or an
ultraviolet no-go. Nor is this a conserved renormalized quantum
source: second-order mean fields, connection variations, full state
construction and renormalization still matter.

The center zero-output Hessians describe the coordinate spatial
average of the pointwise density jet. They are not the complete
local bilinear kernel, or the energy integrated with the perturbed
physical volume. The fixed spatial constraint gauge itself can be
pseudodifferential. No new quantum locality claim is made.

Exact algebra and written arguments are source-pinned, not
proof-assistant formalized or independently peer reviewed.
