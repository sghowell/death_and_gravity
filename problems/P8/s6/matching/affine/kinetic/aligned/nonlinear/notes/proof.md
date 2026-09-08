# Full ADM reduction and local nonlinear constraint count

All formulas use the existing normalized clock u=phi, repo X=s^2=N^-2,
source x=-X, and positive spatial metric. U=e^(3 omega) and
h=(1+u^2)^3. The physical metric is not changed as a matter coupling.
The established point chart is h_physical=e^(2 omega) h_hat, where
omega=-log((h-1+X)/h)/4. Work with positive N and the original closed
X tube [9/10,11/10]. Surface terms below use periodic data or the usual
compact support/falloff; they do not prescribe a new bulk lapse datum.

## 1. Complete original scalar ADM action

With K positive for expansion and V=n(s), the frozen original action is

    N sqrt(h_physical) [ F+B(K^2-K_ij K^ij)+2s f_phi K
                        +C K V+D V^2+E a_i a^i-f R3 ].

Here B=f because A1=0, f=-(h-1+X)/(2h)=-2p_affine^2,
C=s(4f_X+X A3), D=X(A3+A4+X A5), E=4X f_X-X^2 A4,
and A3=1/(hX). The scalar F is the entire frozen function, not its
on-clock Taylor representative. Metric braiding is exactly zero.

Set Q_time=s omega_phi+2s omega_X V. Then

    K_physical=K_hat+3Q_time,
    (K_ij K^ij)_physical=(K_ij K^ij)_hat+2Q_time K_hat+3Q_time^2.

Direct substitution, checked at general X and u, gives

    B(K_hat^2-K_hat,ij K_hat^ij)+b0 K_hat+F0+L_V V,
    b0=4B s omega_phi+2s f_phi,
    F0=F+6B X omega_phi^2+6X f_phi omega_phi,
    L_V=12X f_phi omega_X.

Both the K_hat*V and V^2 coefficients vanish. L_V is independent of
K_hat, shear and V; it need not vanish off the clock. Use the smooth
primitive I(u,s), I_s=U L_V, I(u,1)=0. Since u is constant on each slice,

    sqrt(h_hat) I_s (dot(s)-N^i D_i s)
      = d_t(sqrt(h_hat) I)-d_i(sqrt(h_hat) N^i I)
        -sqrt(h_hat) I_phi-N sqrt(h_hat) I K_hat.

Thus the exact time part becomes N sqrt(h_hat) times

    U B(K_hat^2-K_hat,ij K_hat^ij)+b K_hat+f0,
    b=U b0-I,    f0=U F0-s I_phi.

This removes the remaining lapse velocity by an actual local boundary
term. In particular I_N=0 and I_NN=-18u/(1+u^2)^7 on the clock.
Omitting I changes the background lapse Jacobian by 18 at u=0.
The Taylor jet of I used later is exact through the required second
N derivative, not an off-clock replacement for this full primitive.

The spatial-curvature transformation must also be retained:

    R3_physical=e^(-2omega)[R_hat-4 Delta_hat omega-2(D_hat omega)^2].

Let B4=-f and d=N omega_N=-2X omega_X. After integrating its Laplacian
term, the lapse-gradient-square coefficient divided by e^omega is

    E+2B4 d^2+4(B4-2X B4_X)d = 0.

This exact general-X identity leaves only e^omega B4 R_hat. No lapse
time derivative or spatial derivative remains in the scalar action.
Checking the pre-boundary expression instead would give a false
nonzero lapse-gradient coefficient.

## 2. Exact retained source and trace Legendre map

S6.41 has already eliminated all 56 complementary quotient directions
and all four projective gauge directions with a nonsingular all64
Euler lift. S6.42's full source is purely normal in the clock chart:

    S_normal=delta K_hat+c,
    delta=(s^2-1)/h=4p_affine^2-1,
    c=-3H s delta+(3/2)s Q_lower.

Q_lower is the actual frozen scalar-coefficient ODE solution, not an
independent external source. No spatial source or lapse derivative
is hidden in this formula. With T=(W0-N^i Wi)/N, the retained mass is

    N sqrt(h_hat) U [ (T-delta K_hat-c)^2/(2 gamma_t)
                      -h_physical^ij Wi Wj/(2 gamma_s) ].

The shear term is -U B sigma_ij sigma^ij. Since B<0 and U>0, its five
velocity pivots are nonzero. For the trace, a=(2/3)UB and
p=2 pi_trace/(3 sqrt(h_hat)); the normalized Lagrangian is

    L_trace=a K^2+b K+f0+U(T-delta K-c)^2/(2 gamma_t).

At fixed T, a_total=a+U delta^2/(2 gamma_t), and

    K=(p-b+U delta(T-c)/gamma_t)/(2 a_total).

The temporal Gauss contribution to H is -T j, j=D_i Pi^i/sqrt(h_hat).
The simultaneous trace/temporal solution and Hamiltonian are exactly

    K_joint=(p-b-delta j)/(2a),
    T_joint=delta K_joint+c-gamma_t j/U,
    H_after=(p-b-delta j)^2/(4a)-f0+gamma_t j^2/(2U)-c j.

Before eliminating T, H_TT=-U/gamma_eff, with

    gamma_eff=gamma_t+U delta^2/(2a)
             =gamma_t-3 delta^2/(8 p_affine^2).

The inherited closed-tube bounds p_affine^2 >=9/40, |delta|<=1/10
and gamma_t>18/19 imply a correction <=1/60 and
gamma_eff>1061/1140>9/10. Therefore both the fixed-T trace Legendre
map and the temporal auxiliary pivot are nonzero on that tube. A
negative unreduced gravitational trace coefficient is not by itself
a physical ghost verdict: the constraints still have to be imposed.

## 3. Maxwell and physical free matter

Keep the coordinate W0 independent until after the full Maxwell
Legendre transform. Let G=h_hat^ij, kappa=zeta e^omega>0, and
E_i=dot(Wi)-D_i W0-N^j F_ji. In units of sqrt(h_hat),

    L_electric=kappa E^T G E/(2N),
    P=kappa G E/N,    E=N G^-1 P/kappa.

The electric Hamiltonian is N P^T G^-1 P/(2kappa)+P^i D_i W0,
plus the shift term. All three electric velocity pivots are positive.
Pi^i=sqrt(h_hat) P^i is now an independent momentum density. Spatial
integration by parts gives Pi^i D_i W0=-W0 D_i Pi^i plus a boundary.
Only now use W0=N T+N^i Wi. The normal term is -N T D_i Pi^i, while
the vector part of the momentum constraint is

    H_i(vector)=Pi^j F_ij-Wi D_j Pi^j.

The exact one-form Lie-derivative generator verifies this last term;
dropping it would not generate spatial diffeomorphisms on Wi.
Treating Pi as its earlier velocity expression would obscure the
algebraic nature of the canonical temporal equation, even in Proca.

The magnetic energy is N sqrt(h_hat) zeta e^(-omega) F_ij F_hat^ij/4;
the spatial mass energy is N sqrt(h_hat) e^omega Wi W_hat^i/(2 gamma_s).
The matter Hamiltonian is N sqrt(h_hat)[p_chi^2/(2U)
+e^omega(D_hat chi)^2/2] plus its shift generator. These contain no
derivatives of N or T. The five shear Hamiltonian terms are likewise
algebraic in N. Consequently the entire Hamiltonian, not just its
homogeneous truncation, has the two auxiliaries N,T without their
spatial derivatives. The map W0=N T+N^i Wi is an invertible point map;
its primary-momentum shifts vanish on p_W0=0.

## 4. Actual rolling Jacobian and persistence

On the original rolling background, N=1, T=0, Wi=Pi=0,
K_hat=3H, p=-2H, p_chi=ell and the spatial curvature and shear vanish.
Use the entire original scalar potential in the Hamiltonian above.
Q_lower=Q_lower,N=0 and Q_lower,NN=-3h_phi/h^3 follow from its
original ODE. Differentiate with canonical momenta held fixed.
The code checks H_N=0 and the full joint auxiliary Hessian divided
by sqrt(h_hat):

    D = d_(N,T) d_(N,T) H /sqrt(h_hat) = diag(-2J,-1).

This checks the lapse and temporal coupling jointly, rather than
assuming the isolated mass determines the secondary rank. No Theta
division occurs at the center. The actual J is P34(u)/[800(1+u^2)^18];
P34 has strictly positive even-power coefficients and constant 1199.
For |u|<=1/2,

    J >= 1199/[800(5/4)^18] > 1/40.

Hence both background singular values of D exceed 1/20. All exact
Hamiltonian coefficients are smooth on an open set about the closed
clock tube. The implicit function theorem and compactness of the
background segment give a nonempty open canonical/spatial-jet
neighborhood where N,T remain locally algebraically solvable and D
is invertible. For example, a matrix perturbation norm <1/40 is a
sufficient Jacobian-space condition; no numerical field/jet radius
is claimed. Spatial curvature and momentum divergence may enter D
as physical jet variables, not as derivatives acting on N or T.

## 5. Dirac closure and count

There are 15 configuration variables: six h_hat, one N, three shifts,
four vector components and one chi. The ten nonauxiliary velocity
pivots (five shear, one trace, three vector and one chi) are nonzero.
The primary constraints are p_N,p_T and the three shift momenta.
Preserving p_N,p_T gives chi_N,chi_T. Their mixed bracket is the local
invertible multiplication block D, up to the nonzero density and a
consistent overall secondary sign. Write the four-constraint block as

    C=[ 0  -D ],      C^-1=[ D^-T E D^-1   D^-T ],
      [ D^T E ]             [ -D^-1          0  ].

The secondary-secondary bracket E may be a spatial differential
operator. The displayed inverse is an operator-composition identity:
it does not assume E commutes with D, and does not invert E. Thus
these four constraints are second class and their consistency fixes
the remaining multipliers, rather than adding a tertiary branch in
this neighborhood. No extra lapse boundary condition is inferred.

Spatial covariance is exact with the clock fixed. Its three momentum
constraints, extended by p_N D_i N and p_T D_i T as necessary before
restriction, and the three shift momenta are six first-class constraints.
The vector and chi shift generators above retain all matter terms.
The physical local count is therefore

    [30-2*6-4]/2=7 = 2 tensor+1 clock scalar+3 vector+1 chi.

Removing free chi gives six. Gauge and auxiliary affine directions
already removed by the full-rank S6.41 lift must not be counted again.
This is a local nonlinear count about the compact timelike rolling
branch; it is not a nonlinear energy/stability, global well-posedness,
other-background cone, X=0 vacuum, cutoff, quantum or V/G/B theorem.
Original P8 remains open.
