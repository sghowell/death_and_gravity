# Physical-volume Green operator and field algebra

S6.98 gives a regular density Hamiltonian H_rho and
M=Omega Hess(H_rho), polynomial of degree two in q=k^2/a^2.
Coefficients are real smooth functions of u, with no Theta, Lambda
or spatial inverse in this original phase. Replace q by -Delta/a^2.
Set H=Hess(H_rho) and

    Q=a^-3 (Omega partial_u+H)=a^-3 Omega(partial_u-M).

The fibre pairing is the positive Euclidean Hermitian pairing on the
trivial complex rank-four bundle; integration uses a^3 du dx.
H is a symmetric even-spatial-derivative matrix. Its time coefficients
commute with spatial derivatives. Thus integration by parts gives

    <f,Qg>-<Qf,g> = integral partial_u(f^T Omega g)
                              + spatial boundary terms

for real fields; use conjugate transpose for complex fields.
The time coefficient is antisymmetric, which compensates the
adjoint sign of partial_u. No missing derivative of a^3 occurs:
the factor a^-3 in Q cancels it BEFORE integration by parts.
This Q is the negative normalized action Hessian. That sign is chosen
to match advanced-minus-retarded and the physical positive CCR.

Let U(t,s) be the complete original density transfer. Relative to the
PHYSICAL source volume the kernels are

    G_ret(t,s)=-theta(t-s) U(t,s) Omega,
    G_adv(t,s)= theta(s-t) U(t,s) Omega,
    E_Q=G_adv-G_ret=U(t,s) Omega.

For example the coordinate integral defining G_ret has integrand
-U(t,s) Omega a(s)^3 f(s). Acting with Q gives f: the jump product
is -Omega^2=I and the coincident volume factors cancel.
The homogeneous terms vanish by U'=MU. Applying G_ret to Qf and
integrating the source derivative instead gives f for compact f,
using d_s U=-UM and H=-Omega M. The same calculation with the other
endpoint proves both advanced identities.

The matrix estimate in S6.98 notes/growth.md is for ALL of U, even
though that checkpoint's main observable was chi. Each entry has
the same entire-type bound. Apply its pinned distributional support
argument entrywise. Hence both Green operators have physical causal
support. Smooth compact sources give smooth solutions: spatial
Schwartz estimates follow from the polynomial transfer/derivative
bounds; time derivatives follow from the complete differential
equation. Causal support then gives the usual smooth support class.
Uniqueness follows from the same Fourier Cauchy evolution.

The physical metric is globally hyperbolic with Cauchy surfaces
u=constant. The scale is smooth and at least one, u runs over all
real times, and causal spatial displacement is bounded on every
finite interval by integral du/a. Thus Q and its formal adjoint
(which is Q) have the required Green operators on this spacetime.
A degenerate highest-spatial-order symbol does not invalidate this
direct construction or justify using normally-hyperbolic theorems.

Canonical symplectic conservation gives the two-time field
commutator i*hbar U(t,s) Omega/kappa. It agrees with i*hbar E_Q/kappa
under the physical-volume convention. The state below is therefore
a state of this spacetime field algebra, not merely unrelated
oscillators. Projection onto chi gives the same normalized Kubo
response kappa*a(s)^3*i/hbar as the frozen causal-support calculation.
