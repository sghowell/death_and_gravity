# Canonical density evolution and normalization

Use the actual S6.91 phase X=(v,chi,p,P), with weighted generator
M_X=Omega Hess(H2)-diag(0,0,3H,3H). Its coefficients are evaluated
along the S6.88 analytic constrained solution, not frozen at u=0.
The physical comoving source density is V=R^3 N e^3.

Set Z=D_R X, D_R=diag(1,1,R^3,R^3). Its conjugate density momenta
are pi_v=R^3 p and pi_chi=R^3 P. The literal Hamiltonian density is
H_rho=R^3 H2(v,chi,pi_v/R^3,pi_chi/R^3). Differentiating D_R gives
M_rho=(D_R'+D_R M_X)D_R^-1=Omega Hess_Z(H_rho).
The native identity retains D_R' and the weighted original damping.
Consequently M_rho Omega+Omega M_rho^T=0.

The fundamental solution U(t,s,k) exists uniquely on I for every
real k-vector. Its generator is a smooth polynomial in |k|^2,
with degree at most four in the Cartesian momentum components.
No inverse Legendre Hessian and no 1/k chart are needed to define it,
including at k=0. The k=0 extension is a smooth Fourier multiplier;
no separate homogeneous constrained perturbation theorem is claimed.

Differentiate U Omega U^T. It solves Q'=M_rho Q+Q M_rho^T
with Q(s)=Omega. The constant Omega solves the same equation;
ODE uniqueness gives U Omega U^T=Omega. The generator and transfer
are real and even in k-vector. This also fixes the -k argument
in position-space Fourier commutators without an unspoken reality
or momentum reversal assumption.

For k=|k-vector|>0 the exact S6.90 packet phase is Y=E Z:
Y=R^(3/2)(-kp/(2q), k chi, 2qv+Bp/2, P), q=k^2/R^2.
E Omega E^T=k Omega; the packet phase is not a unit canonical
density phase. Both inverse maps, this scaled symplectic form,
and E V e_4=R^(3/2)N e^3 e_4 are checked exactly.
The observable is O=e_2^T Z=chi, the linear unitary-gauge form
of the relational scalar chi-chi_background(phi).

Restore an arbitrary positive constant action prefactor kappa:
the action is kappa times the normalized quadratic action.
Physical canonical density momentum is kappa*pi; using Z as above
therefore means [Z_a(x),Z_b(y)]=i*hbar*Omega_ab*delta(x-y)/kappa.
The normalized source term is +kappa*integral V J O, yielding
Z'=M_rho Z+V e_4 J. If instead the source term is +integral V j O,
then j=kappa J and its response is smaller by 1/kappa.
No earlier physical mass/time/Planck-scale dictionary is assumed.
