# Nonlinear D-only chart and lapse-reduced interaction generator

All formulas in sections 1--4 use M=tau=1 and u=t, d=1+u^2. Section 5 restores
dimensions. They apply to the exact D witness, not arbitrary DHOST/matter.
The computations retain full spatial tensor contractions; no Fourier-mode
average or homogeneous truncation defines the nonlinear action.

## 1. Removing lapse derivatives exactly

The Ia completion with f=F2=-1/2, A1=0, A3=4/d^3 gives

    A4=-A3+X^2 A3^2/4,    A5=-X A3^2.

In unitary clock phi=t, X=N^-2. Set V=n(1/N), where n is the unit normal;
positive-expansion extrinsic curvature has K^i_j=H delta^i_j on FLRW.
The ADM density modulo the S1 curvature boundary is

    N sqrt(h) [F + (Kij Kij-K^2)/2 + C K V - 3 C^2 V^2/4
                 + R[h]/2 + E |D ln N|_h^2],
    C=A3/N^3,  E=-X^2 A4.

Use an auxiliary spatial metric, leaving the physical lapse and shift fixed:

    h_ij = exp(2 omega) hat_h_ij,
    omega(u,N)=(N^-4-1)/(2d^3).

This has the explicit inverse hat_h=exp(-2 omega)h for N>0. The Jacobian in
(N,h_ij) is triangular with nonzero diagonal; there is no condition on
Theta or Lambda. It is a field chart, not a replacement of the physical
matter metric. At N=1, omega=0, omega_N=Delta=-2/d^3, reproducing zeta=v+Delta n.

Since omega_N=-C/(2N^2),

    K_old^i_j = hat_K^i_j + (C V/2+W) delta^i_j,    W=omega_u/N.

Completing the full kinetic square gives exactly

    (hat_Kij hat_Kij-hat_K^2)/2 - 2 W hat_K - 3 W^2.

The spatial identity in three dimensions is

    R[h]=exp(-2 omega)(R[hat_h]-4 hat_D^2 omega-2|hat_D omega|^2).

The divergence from integration by parts of -2 N exp(omega) hat_D^2 omega
is -2 hat_D_i(N exp(omega) hat_D^i omega), with the spatial volume factor
understood. The remaining |hat_D ln N|^2 coefficient is

    E+(N omega_N)^2+2N omega_N = 0.

Thus both spatial and temporal lapse derivatives disappear **to all orders**:

    L/(N sqrt(hat_h)) = exp(3 omega)
      [F+(hat_Kij hat_Kij-hat_K^2)/2-2W hat_K-3W^2]
      + exp(omega) R[hat_h]/2.

Spatial boundaries vanish for periodic fields or suitable spatial decay.
The original covariant-to-ADM temporal boundary is the one already displayed
and checked in S1; do not discard additional time-dependent canonical
boundaries in a later change to physical perturbation coordinates.

## 2. Full canonical Hamiltonian

Define p^i_j=pi^i_j/sqrt(hat_h) and p=tr p. The six-component momentum map is

    p^i_j = exp(3 omega)/2
            [hat_K^i_j-hat_K delta^i_j-2W delta^i_j].

The Hamiltonian, with the usual spatial momentum constraints retained, is

    H = sqrt(hat_h) [2 N exp(-3 omega)(tr(p^2)-p^2/2)
          -2 omega_u p -N exp(3 omega)F(u,N^-2)
          -N exp(omega) R[hat_h]/2] + N^i H_i,
    H_i=-2 hat_h_ij hat_D_k pi^{jk}.

Off-diagonal components are counted twice in matrix contractions; code
checks the Legendre transform with all six independent symmetric components.
This formula includes scalar and tensor interactions before the momentum
constraints are solved. Lapse has no conjugate velocity and is algebraic.

Write p^i_j=-H delta^i_j+delta p^i_j, sigma=tr(delta p),
S2=tr[(delta p-sigma I/3)^2], and rho=R[hat_h]. Then

    tr(p^2)-p^2/2=-3H^2/2+H sigma-sigma^2/6+S2.

S2 denotes an invariant of degree two, not a second independent scalar
field. rho is the full nonlinear spatial curvature. Denote the square-bracket
Hamiltonian by h(N,sigma,rho,S2).

## 3. Independent bridge to the linear checkpoint

Differentiation at N=1, sigma=rho=S2=0 gives

    h_N=0,
    h_NN=-2 J,   h_Nsigma=2 Theta,   h_Nrho=-Lambda/2,
    J=2(u^8+2u^6+4u^2+1)/d^5,
    Theta=2u/d,   Lambda=1-2/d^3.

These are derived from the nonlinear expression, not inserted from the
linear reduced action. The replay checks their agreement with S0--S4.
Exact rational Sturm tests also establish

    1/(2d) < J <= 2/d.

The upper inequality follows directly because d^4 minus the numerator
polynomial of J/2 is 2u^6+6u^4. In particular h_NN is nonzero at every
finite time, including gamma and both Lambda zeros. The ordinary implicit
function theorem supplies a unique local analytic lapse branch about the
background at each finite time. Compact time intervals admit finite local
covers. This is an algebraic statement in local curvature/momentum inputs,
not a well-posedness theorem for the full nonlinear field PDE, a global
large-amplitude lapse solution, or a uniform-in-time EFT neighbourhood.

The local Dirac count in this chart is also transparent: the lapse's zero
momentum and its algebraic equation form a second-class pair where h_NN!=0.
The three spatial diffeomorphisms retain their usual first-class momentum
constraints. Excluding the shift multiplier pairs, seven configuration
variables (six metric plus lapse) therefore leave
(14-2-2*3)/2=3 local propagating degrees of freedom. This uses spatial
covariance and the usual constraint-counting lemma; it is not a formal
Poisson-bracket/PDE proof or a nonlinear energy-positivity assertion.

At u=0, J=2 and h_NN=-4. The small/zero lapse denominator is **not** Theta.
At |u|->infinity, J~2/u^2; no uniform unweighted inverse-Hessian bound exists.

## 4. Cubic and quartic invariant generator

Expand n=N-1 with weights sigma=rho=1 and S2=2. At fixed spatial metric,

    h = A(n)+L(n)+Q(n),
    L_k=B_k sigma+C_k rho,    Q_k=D_k sigma^2+E_k S2,

where k denotes the k-th N derivative at N=1, not a Taylor coefficient.
A1=0. Define

    n1=-L1/A2,
    f2=A3 n1^2/2+L2 n1+Q1,    n2=-f2/A2,
    n3=-(A3 n1 n2+A4 n1^3/6+L2 n2+L3 n1^2/2+Q2 n1)/A2.

Here A_k are Hamiltonian jets, **not** the covariant DHOST A_i functions.
The stationary Hamiltonian pieces are

    h0=A0,  h1=L0,
    h2=Q0-L1^2/(2A2),
    h3=A3 n1^3/6+L2 n1^2/2+Q1 n1,
    h4=A4 n1^4/24+L3 n1^3/6+Q2 n1^2/2-f2^2/(2A2).

The final term in h4 is the second-order lapse correction. Substituting only
the linear lapse into a quartic Hamiltonian would miss it. Generic symbolic
polynomial substitution checks the constraint through degree three and the
stationary Hamiltonian through degree four. The report stores all model
coefficients, their finite-time denominator proofs and exact bounce values.

As a compact exact regression at the bounce, with Z=sigma^2-6S2,

    n1=rho/8,    n2=-7Z/12+97rho^2/128,
    h2=-Z/3+rho^2/32,
    h3=-7rho Z/24+29rho^3/256,
    h4=49Z^2/72-697rho^2 Z/384+25553rho^4/24576.

An additional test differentiates the literal nonlinear bounce action with
F(0,X)=-19+20X-5X^2, independently of the model's stored lapse-jet array.

These are invariant Taylor coefficients of H/sqrt(hat_h). They are not yet
physical cubic/quartic vertices: sqrt(hat_h), rho and momentum index raising
must still be expanded in the actual perturbations, the spatial momentum
constraints solved, and the physical modes canonically normalized. There
is no spatial collinearity restriction, but no scattering amplitude yet.

## 5. Physical scales and the remaining obstruction to a cutoff claim

Give the clock phi dimensions of length, so X is dimensionless. For
u=t/tau, phi=t, restore

    F=M^2/tau^2 f(phi/tau,X),  F2=-M^2/2,  Ai=M^2 ai(phi/tau,X),
    a=1+(t/tau)^2,  H=2u/(tau d),  varphi=2M asinh(phi/tau).

All Li scale as tau^-2; x^mu=tau y^mu and phi=tau bar_phi give
S=(M tau)^2 S_dimensionless. The classical family therefore has a tunable
dimensionless prefactor. At a fixed normalized frequency, canonical n-point
vertices scale as (M tau)^(2-n), when the relevant kinetic chart is regular.
This scaling alone does not bound frequency growth, tail normalization or
all interaction orders, and does not determine a strong-coupling scale.

For a nonzero reference at the bounce set E_ref=1/(tau sqrt(d)). Exactly,

    E_curv=sqrt(H^2+|Hdot|),    sqrt(2) E_ref <= E_curv < sqrt(6) E_ref.

Also |d log(A3)/dt|<=6 E_ref. This is not a bound on every canonical
interaction's variation; logarithms of coefficients that cross zero should
not be mistaken for physical singularities.

The old regular b-chart finite-q kinetic coefficient is

    K_b(q)=q J/(q Lambda^2-J),   K_b(infinity)=J/Lambda^2,
    K_b(q)/K_b(infinity)-1=J/(q Lambda^2-J),

with dimensionless q=tau^2 k_physical^2. On |u|<=1/4, |Lambda|>1/2, J<=2.
Consequently for q>=64 the kinetic-chart relative error is <=1/7, uniformly
through gamma. This checks a **kinetic** high-frequency limit only; it is
not an adiabatic bound for all terms or a weak-coupling window. Momentum
coordinates can fail at q Lambda^2=J while the canonical Hamiltonian stays
regular. Such a chart failure is not itself a physical instability.

## 6. Closing the algebraic tail-normalization gap

The unweighted J->0 observation does not by itself establish strong coupling.
At each background time u0, choose a **constant local unit of length**
ell=tau sqrt(1+u0^2). Define x=u0/sqrt(1+u0^2) and z=1-x^2. This compactifies
the two infinite tails to x=-1,+1. The dimensionless local invariants are

    sigma_local=sqrt(d) sigma,   rho_local=d rho,   S2_local=d S2,
    h_local=d h,               H_local=sqrt(d) H=2x.

The action prefactor in these units is (M ell)^2=(M tau)^2 d. This is not
a global time-dependent change of canonical variables; no kinetic terms
from derivatives of ell(t) have been silently dropped.

The pinned overlapping **principal** kinetic charts also have uniform
normalization in these units. For |u|>=1/8 the curvature variable has
K_v/M^2=P/(2x^2), hence 1/8<K_v/M^2<=65/2. For |u|<=1/4, beta=b/ell gives
K_beta/M^2=2P/Lambda^2, hence 1/2<K_beta/M^2<=8. Gradients have the same
principal coefficients by the old exact K=G identity. For the gamma chart,
ell^2 k_physical^2>=68 ensures the section-5 finite-q kinetic bound because
d<=17/16. These statements concern quadratic principal normalization, not
nonlinear canonical vertices or all subprincipal/adiabatic terms.

The complete local Hamiltonian has the same form as section 2 with

    omega=z^3(N^-4-1)/2,
    omega_time=-3x z^3(N^-4-1),
    F_local=F0_local+FX_local(X-1)+FXX_local(X-1)^2/2,
    F0_local=-4(1+x^2),
    FX_local=2(2x^2-1)(1-6z^3),
    FXX_local=-2z^3(45x^8+122x^6 z+96x^4 z^2-12x^2 z^3+5z^4).

Here omega_time means ell times the original physical partial time
derivative, not differentiation while artificially holding x constant.
All functions are regular jointly in x in [-1,1] and N near one. Direct
independent differentiation of this compact nonlinear Hamiltonian gives

    h_local,N=0 on the background,
    h_local,NN=-4P(x),    P(x)=1-6x^4+10x^6-4x^8,
    1/4<P(x)<=1,          -4<=h_local,NN<-1.

The P bound follows from the exact weighted-J proof and the endpoint
values P(-1)=P(1)=1. Joint analyticity, the compact x parameter set and
the uniformly nonzero lapse Hessian imply a uniform **local algebraic**
lapse branch for sufficiently small normalized invariant inputs. This
closes the tail gap for lapse solvability, not for physical scattering or
the full nonlinear evolution. It is the parameter-uniform implicit-function
lemma applied to this displayed algebraic function.

For a monomial sigma^p rho^r S2^w, its normalized Hamiltonian coefficient
is c(u)d^(1-p/2-r-w); for the lapse expansion it is c(u)d^(-p/2-r-w).
All 34 coefficients through the recorded orders become rational functions
of x, with denominators only a nonzero constant times P(x)^k. The report
stores each function, both exact endpoint limits, and the conservative bound

    abs(coefficient) <= 4^k sum(abs(numerator power coefficients))
                         / abs(denominator constant).

This bound uses abs(x)<=1 and P>1/4; there is no numerical sampling or
unproved tail extrapolation. In particular a raw lapse-rho coefficient
growing as -u^2/8 instead tends to -1/8 in these local units. The bound is
not a canonical coupling: the spatial momentum constraints, metric volume,
geometry, physical normalization and exchange/contact terms remain as
listed in the S5 contract. Finite algebraic coefficients at both tails do
not establish a numerical cutoff, a scattering vacuum or UV completion.
