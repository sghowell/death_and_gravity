# S6.277: full homogeneous fold dynamics in the unchanged local parent

Keep the original S6 contract, complete current coefficient functions, all
three vacuum constants, fixed Proca and both scalar stress profiles, masses,
primitive and canonical source. This is a separate classical comparison
result for that literal local action. No quantum state, UV completion,
regulator prescription or old certificate is changed.

Work in unitary clock gauge u=t with positive lapse N, and in the diagonal
homogeneous flat Bianchi-I invariant sector. The full lapse and temporal
equations are retained. Let

    gamma=e^(2alpha)diag(e^(2beta),e^(-2beta),1),
    p=Palpha/(3 kappa e^(3alpha)),
    m=PM/(kappa e^(3alpha)),
    s=Pbeta^2/(8 kappa^2 e^(6alpha)),
    eta=10^100 h, ph=PH/(kappa e^(3alpha)).

All homogeneous spatial Proca coordinates and momenta are zero; the temporal
vector is reconstructed by its full equation. The heavy scalar is initially
zero but remains live. Spatial momentum constraints vanish identically.
M1 and beta are cyclic coordinates with full quadrature reconstruction.

Write rho,P and rho1,P1,rho2,P2 for the actual fixed energy/pressure jets at
u0. Each is bounded in magnitude by epsilon=2/10^400. They are not free
retuning parameters or means of the new classical solution.

At u0,N1,eta=ph=0, every real p in the family

    s=81/160-2P/3+7rho/12+3p^2/8,
    m^2=1217/200-2P+3rho

satisfies the full C=H_N=0 and C_N=0. Here s>49/100, m^2>6,
C_s=3 and D=C_NN=(11391+1400P-1600rho)/400>28.
Either sign of m satisfies these equations. The future endpoint uses m>0.

The full first consistency numerator is

    K=C_u+{C,kappa e^(3alpha)H}
     =p[-3/25+3(P-rho)]+3P1/2-rho1.

The source derivative includes B_uN=-6. When K=0, its unique p satisfies
|p|<21epsilon. A C1 lapse solution would also require

    a N_u^2+b N_u+c=0,

where the ENTIRE second-preservation coefficients in notes/consistency.md
obey a>28, |b|<1 and c>60. Its discriminant is below-6719.
Thus no C1 clock-time lapse solution passes through this specified family.
This is not an assertion that C=0 is an irregular constraint surface:
its gradient C_s is nonzero even though its auxiliary Poisson rank changes.

For p=-1/10 one has J=K>11/1000. The full desingularized vector field

    u_tau=C_N, x_tau=C_N F(u,N,x), N_tau=-K

is tangent to C=0. Its smooth local curve through the fold obeys

    u=-DJ tau^2/2+O(tau^3), N=1-Jtau+O(tau^2).

For small tau>0, C_N<0 and u<0; inverting u gives regular solutions of
the full original canonical equations on that side, approaching u0 in
finite proper time. This is local existence, not a quantitative lifetime.

The actual physical lapse is N, NOT N*R^(-1/4). The physical spatial
metric is R^(-1/2)gamma and its volume is e^(3alpha)R^(-3/4).
For its normal expansion theta and physical Ricci scalar Rphys,

    lim_(tau->0+) tau theta=3/(2D),
    lim_(tau->0+) tau^3 Rphys=-3/(D^2 J).

The physical curvature therefore diverges while lapse and volume have
positive finite limits. No regular C2 physical-metric continuation through
that endpoint is furnished. The auxiliary parameter tau is not a regular
physical clock change at the endpoint; X=1/N^2 tends to1, not0.

All these initial homogeneous momenta differ from the original prepared
state and from S276's mean-preserving Fourier path. S275's corrected finite
hybrid is unchanged. No generic near-bounce instability, quantum mean,
unlocalized Hamiltonian, cutoff-valid singularity, physical UV matching,
omitted-loop/Regge bound or original V/G/B/P8 closure is established.
Written proofs are supported by exact and independent diagnostics, not
by formal verification of every analytic argument.
