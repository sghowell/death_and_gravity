# CD/M1 nonlinear reduction and uniform invariant bounds

The following written lemmas accompany exact symbolic checks. They are not
formalized in a proof assistant and do not establish physical interaction control.

## 1. Full nonlinear chart, including the extra boundary

Write `R=1/d^3` (this R is a coefficient, not spacetime curvature),
`T=1+R(N^-2-1)` and `omega=-log(T)/4`. On the entire frozen tube,
`9/10<=T<=11/10`; lapse is positive and bounded away from zero. Hence the
spatial rescaling and all its finite derivatives are smooth at every finite u,
including the bounce. Its background value is the identity. No inverse map
requires a nonzero Theta or Lambda.

The pinned covariant-to-ADM derivation gives, per `N*sqrt(h)`,

`F+T*(K_ij K^ij-K^2)/2+L0*K+c*K*V+Dv*V^2+T*R3/2+E_acc*a_i*a^i`,

where `c=-R/N`, `L0=-T_u/N`, `Dv=-3c^2/(4T)` and
`V=-D_t N/N^3`. Derivatives with subscript u hold N fixed. The source of L0
is `2*sqrt(X)*F2_phi`; it vanishes identically for D/M0 but not for CD/M1.

For `h_ij=e^(2omega)*hat_h_ij`,

`K_old^i_j=hat_K^i_j+[cV/(2T)+W]*delta^i_j`, `W=-T_u/(4NT)`.

Completing the exact full kinetic square leaves

`T*(hat_K_ij hat_K^ij-hat_K^2)/2+(-2TW+L0)*hat_K`

`-3TW^2+3L0W+3L0cV/(2T)`.

The V-squared and V-times-hat-K terms cancel by Ia degeneracy. The last,
linear-V term does **not** cancel. For the spatial terms, integrating the
conformal curvature Laplacian by parts differentiates `N*e^omega*T`. The
remaining coefficient is

`E_acc+T*(N*omega_N)^2+2*(T+N*T_N)*(N*omega_N)=0`.

The `T_N` contribution is essential. These identities use the full Ia
coefficients, not only their background Taylor series.

Define the local coefficient function

`Psi(u,N)=integral_1^N [-3*R*T_u/(2*n^4*T^(7/4))] dn`.

The integration variable is the lapse argument; all time coefficients are at
the same spacetime point. This is not integration over physical time, spatial
points or past history. The primitive is smooth and is fixed by `Psi(u,1)=0`.
The residual density is `sqrt(hat_h)*Psi_N*D_t N`. The shift-aware identity

`partial_u(sqrt(hat_h)*Psi)-partial_i(sqrt(hat_h)*N^i*Psi)`

`=sqrt(hat_h)*(Psi_u+Psi_N*D_t N+N*Psi*hat_K)`

reduces it, modulo a time/spatial divergence, to
`-sqrt(hat_h)*Psi_u-N*sqrt(hat_h)*Psi*hat_K`.

Therefore the complete gravitational Lagrangian density per `sqrt(hat_h)` is

`N*[A*(hat_K_ij hat_K^ij-hat_K^2)/2+P*hat_K+U+T^(3/4)*rho/2]`,

`A=T^(1/4)`, `P=-T_u/(2NT^(3/4))-Psi`,

`U=T^(-3/4)*[F+9T_u^2/(16N^2T)]-Psi_u/N`, `rho=R3[hat_h]`.

## 2. Metric and matter Legendre transform

Keep all six independent symmetric metric velocities. With
`p^i_j=Pi^i_j/sqrt(hat_h)`, their conjugate momenta satisfy

`p^i_j=A*(hat_K^i_j-hat_K*delta^i_j)/2+P*delta^i_j/2`.

The Legendre map is invertible since A is nonzero. It is not positive definite
before imposing constraints: the usual gravitational trace direction remains.
The free physical matter action, in these variables, is

`L_chi/sqrt(hat_h)=T^(-3/4)*(D_t chi)^2/(2N)-N*T^(-1/4)*z/2`,

with `z=hat_h^ij*partial_i chi*partial_j chi`. Thus its momentum is
`r=pi_chi/sqrt(hat_h)=T^(-3/4)*D_t chi/N`.

The full Hamiltonian density per `sqrt(hat_h)` is

`h=2N[tr(p^2)-(tr p)^2/2]/A+NP*tr(p)/A-3NP^2/(4A)-NU`

`-N*T^(3/4)*rho/2+N*T^(3/4)*r^2/2+N*T^(-1/4)*z/2`.

`legendre_checks()` obtains this independently by differentiating all seven
velocities with the correct off-diagonal metric multiplicities. The retained
spatial constraint is

`-2*hat_h_ij*hat_D_k Pi^jk+pi_chi*partial_i chi=0`.

It must be solved before interpreting this as physical scalar/tensor vertices.

Set `tr(p)=-3H+sigma`, `r=l+eta`, `l=1/(10d^6)` and let shear2 be the
traceless momentum square. Then

`tr(p^2)-(tr p)^2/2=-3H^2/2+H*sigma-sigma^2/6+shear2`.

The invariant polynomial before lapse elimination is
`A0(N)+B(N)*sigma+C(N)*rho+L(N)*eta+D(N)*sigma^2`
`+E(N)*shear2+M(N)*eta^2+Z(N)*z`.
Here A0 is the invariant-independent piece, **not** the kinetic coefficient A.

The newly derived lapse jets obey the independent old-action bridges

`h_N|bg=0`, `h_NN|bg=-2J`, `h_Nsigma|bg=2Theta`,

`h_Nrho|bg=-Lambda/2`, `h_Neta|bg=-w`, `w=l*(3delta-1)`.

In particular `J=Sigma_total+3Theta^2-w^2/2`; forgetting the matter Legendre
subtraction would give the wrong constraint denominator.

A separate boundary regression has

`A0_NN[full]-A0_NN[omitted]=9HRR_u+3(R_u^2+RR_uu)=18(u^2-1)/d^8`.

At the bounce, Psi is zero but Psi_u is not. Direct integration there gives
`Psi_u=24-6N^(-3/2)-18sqrt(N)` and the omitted-boundary J residual is 18.
All first-jet bridges still pass if the boundary is dropped; checking only them
would miss the error. The independent literal bounce expressions verify every
one of the 40 recorded jets through fourth N derivative.

## 3. Stationary expansion and local units

Assign weights one to sigma, rho, eta and two to shear2, z. Before eliminating
lapse h has invariant degree at most two. Its N jets are ordinary derivatives.
The generic stationary-series lemma from pinned S5.1 applies because `A0_N=0`
and `A0_NN=-2J` is nonzero. It gives lapse orders one through three and
Hamiltonian orders zero through four, retaining the second-order lapse
correction in the quartic Hamiltonian. The same generic constraint/stationarity
identities are replayed, now with all matter monomials included.

Taylor-jet convolution and binomial recurrence derive the coefficients exactly.
An independent symbolic differential identity validates the fractional-power
jets. A fourth-order N jet of Psi suffices: u differentiation at fixed N does
not lower its N-1 order. There is no time series or Fourier truncation here.

Restore units with `F_phys=(M^2/tau^2)*F`, `F2_phys=M^2*F2`,
`Ai_phys=M^2*Ai`, `phi_phys=tau*u` and `chi_phys=M*chi_bar`. Rescaling spacetime
coordinates by tau gives overall action prefactor `(M*tau)^2`. This scaling
alone does not prove an interacting hierarchy.

At each background point choose **fixed** local units `ell=tau*sqrt(d)`.
Set `x=u/sqrt(d)`, `y=1/sqrt(d)`, so `x^2+y^2=1` and both tails lie in the
compact closure `-1<=x<=1`, `0<=y<=1`. In dimensionless original units,

`Hbar=4x`, `lbar=y^11/10`, `R=(1-x^2)^3`,

`R_u,bar=-6x(1-x^2)^3`, `R_uu,bar=6(8x^2-1)(1-x^2)^3`.

Every coefficient of `d*F` is a polynomial in x. Every rescaled lapse jet is
therefore polynomial in x,y, reduced by `y^2=1-x^2`. Symbolic homogeneous
scaling checks and independent substitutions verify the dictionary. This is a
pointwise constant choice of units, not a global time-dependent canonical map.

Exact all-real-time FLINT sign certificates give `1/10<d*J<8`. Consequently
`-16<hbar_NN<-1/5` and `abs(1/hbar_NN)<5`. On the compact square, the sum of
absolute polynomial coefficients bounds each jet. Every stationary coefficient
is a polynomial in those jets divided only by a power of hbar_NN. A second
triangle bound with inverse bound 5 therefore gives explicit uniform rational
majorants for all 27 lapse and 58 Hamiltonian coefficient records. Denominators
with any other variable factor are rejected by the verifier.

These bounds do not solve the matter-sourced spatial constraint, normalize the
physical modes, include their time-dependent canonical generators, calculate
exchange/contact amplitudes, or establish a perturbative frequency window.
