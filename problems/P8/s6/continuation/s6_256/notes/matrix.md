# Complete finite-band matrix and mass-adapted energy chart

Enter S6.254's ENTIRE four-mode central Hamiltonian using the rolling
S6.255 map. Only actual classical Euler contacts are removed:
L0,Vvv,vs_M1,vs_H vanish after variation. The identities d_M1=0,
mass_M1,M1=mass_M1,H=0 hold; mass_H,H=-mu with mu=N U n>0.
Keep c_H,w_H,d_H, all Proca source contacts, the full lapse square and
the nonzero heavy mass. No finite off-clock source is set to zero.

The central physical canonical transformation T0 is the complete
S6.254 map, including its -Hhat Qb Pb time term. With L0=0 the
old-from-new perfect-square shear is

    Pb_old=Pb_new-4 a_hat Theta P^2 Qb/(D L2).

Call its matrix R. The full new Hessian is
R^T H R+Omega R^-1 R_dot. Compute R_dot with a_hat_dot=a_hat Hhat
and the complete Theta_dot,D_dot,L2_dot; do not freeze coefficients.
The shear is symplectic and this connection is symmetric.

Reorder to (Qb,sigma_M1,QL,Pb,p_M1,PL,hbar,p_H). In the light
six-phase block use momentum weights

    Mlight=diag(1,1,P^1/2,P^3/2,P^1/2,P).

The entire weighted matrix is split into P^3/2 Afast_light, the
retained massive heavy oscillator, and the full remaining matrix.
There is no limit in P or n in this subtraction. Define

    alpha=L2^2/(8 J a_hat^3),
    g=r sqrt(zeta)/(D a_hat^2),
    h=a_hat^3 Yv/zeta,
    beta=L2 w_M1/(4 J Z a_hat^3),
    charge=a_hat c_M1/D,  d=a_hat Y,
    rho^4=alpha g^2 h.

The evaluated coefficient box gives 10^-4<rho<10^-2,
1/20<alpha<1/10, -10^-8<g<-10^-10, 10^5<h<10^7,
|beta|,|charge|<1 and 1/2<d<2. Its actual first time derivatives
are below 10^50; rho_dot/rho is obtained by differentiating its
positive fourth-power relation.

The complete light growth chart is

    (Qb, -rho QL/g, alpha Pb/rho, rho^2 PL/(g h),
     sigma_M1-beta Qb/alpha,
     [p_M1-(charge+d beta/alpha) QL/g]/80000).

It sends Afast_light EXACTLY to rho times the four-cycle

    C4=[[0,0,1,0],[1,0,0,0],[0,0,0,1],[0,1,0,0]]

plus the retained two-phase slow block [[0,0],[-d/80000,0]].
The polynomial remainder modulo rho^4-alpha g^2 h is checked in
all 36 entries. The slow block is not assumed diagonalizable.

Keep the heavy core

    Aheavy=[[0,Bh],[-Ah,0]],
    Ah=a_hat^3 mu+a_hat Y P^2,  Bh=1/(Z a_hat^3).

Its positive energy chart diag(sqrt(Ah),sqrt(Bh)) turns it into
sqrt(Ah Bh) times the standard skew matrix. This is an exact
identity at the fixed giant mass. Its frequency is not a norm error
and is not deleted. Positive square roots are normalized by exact
positive factorization before structural zero checks.

Combine the light chart and heavy energy chart as W. Bound every
entry of

    W [M^-1 Awhole M-diag(P^3/2 Afast_light,Aheavy)] W^-1 / P^3/2

on the entire coefficient box and 10^64<=P<=2*10^64, with all
physical time jets up to 10^40, c_H,w_H up to 10^-1000, d_H up to
10^-880 and mu in [n/2,2n]. Exact outward rational enclosures give
a sum of absolute entry bounds below 10^-19. This bounds the
operator two-norm of the ENTIRE eight-by-eight remainder. No raw
mass-squared estimate, finite subset of entries or leading-term
guess substitutes for this computation.

Differentiating the complete light growth chart gives
||S_dot S^-1||<10^63, hence its normalized contribution is below
10^-33 since P_min^3/2=10^96. These are actual coefficient/time
bounds, not independent arbitrary parameters silently assumed
to describe a solution.

The two tensor polarizations have positive canonical cores
A_TT=2 a_hat C P^2, B_TT=1/(a_hat^3 D).
For determinant-one tensor perturbations exp(e q), e_ij e_ij=4,
the ADM density is a_hat^3 D q_dot^2/2-a_hat C P^2 q^2.
The trace is unchanged, so homogeneous matter and the aligned
source introduce no additional tensor potential. The literal
canonical Legendre block is checked.

Each transverse Proca polarization has
A_PT=a_hat K Y P^2/(zeta Z)+a_hat^3 Yv/zeta,
B_PT=zeta/(a_hat^3 K).
Independent Maxwell ADM re-entry gives electric coefficient
a_hat R^-1/4/N, magnetic coefficient N R^1/4/a_hat, and mass
coefficient N a_hat R^-1/4/zeta, exactly matching this full core.

Use diag(sqrt(A),sqrt(B)) for each remaining positive oscillator.
Their cores are skew and their logarithmic diagonal chart derivatives
are all below 10^50. Because the full connection is block diagonal,
its operator two-norm is the maximum diagonal magnitude, not the sum
over mode copies. The normalized bound is 10^-46. This includes
the heavy chart and both copies of TT and transverse Proca.
Thus all eight physical modes and all sixteen phases are retained.
