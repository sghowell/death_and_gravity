# Actual time-dependent growth, not a frozen-eigenvalue verdict

Fix a compact time interval I of positive length on which every stated
coefficient is smooth, a,D,Z,J,K,Y,Yv,zeta,gamma have strict positive
margins, and r,L2 have nonzero margins. The preceding exact remainder
gives ydot=P^(3/2)[A0(t)+O(P^-1/2)]y with a uniform C1 remainder.
Set rho=(alpha g^2 h)^(1/4)>0 and rho_min=min_I rho>0.

There is an explicit smooth real change of coordinates, with bounded
inverse and bounded time derivative on this fixed interval:

    u=(Qb,-rho QL/g,alpha Pb/rho,rho^2 PL/(g h)),
    sbar_i=sigma_i-beta_i Qb/alpha,
    pbar_i=ps_i-(charge_i+d beta_i/alpha)QL/g.

The entire8x8 principal matrix in these coordinates is a direct sum of
rho C and a two-pair nilpotent block N:

    C(u1,u2,u3,u4)=(u3,u1,u4,u2),
    N(sbar,pbar)=(0,-d sbar).

This identity is checked modulo the exact defining relation
rho^4=alpha g^2 h. In particular, the two slow matter pairs were not
assumed to have a diagonalizable zero eigenvalue, nor were their full
couplings dropped to isolate a handpicked four-dimensional subsystem.

C is a real orthogonal four-cycle. Its unit positive eigenvector is
eplus=(1,1,1,1)/2 and its negative one is eminus=(1,-1,-1,1)/2.
Exactly, (C+C^T)/2=eplus eplus^T-eminus eminus^T. Hence the fast complement
orthogonal to eplus has nonpositive Euclidean energy growth. Write its
positive component as f and include all other components in z.

For the slow block rescale pbar by the fixed positive number
w=rho_min/(4 max_I d). The largest eigenvalue of its symmetric part is
w d/2<=rho_min/8. The algebraic energy margin is the positive square
w d |sbar+pbar_weighted|^2/2. The complete principal complement therefore
has logarithmic norm at most rho_min/8. The additional four physical
modes have zero fast symbol and enter the same complement; their O(P)
terms are part of the full remainder.

The time derivative of this smooth growth chart is O(1), uniformly on I.
Including it, all full remainder terms and their conjugations gives

    fdot=P^(3/2)rho(t)f+E11 f+E12 z,
    zdot=P^(3/2)B0(t)z+E21 f+E22 z,
    ||E||<=C0 P,  lambda_max[(B0+B0^T)/2]<=rho_min/8,

for P above the finite chart/remainder threshold. Increase C0 to at least1.
For sqrt(P)>=8 C0/rho_min and initial f>0,z=0, the cone ||z||<=f is forward
invariant. Indeed, on its boundary the lower derivative estimate for f is
P^(3/2)rho_min f-2C0 P f, whereas the upper Dini derivative of ||z|| is at
most P^(3/2)rho_min f/8+2C0 P f. Their difference is strictly positive.
Inside the cone, fdot>=rho_min P^(3/2) f/2, so

    f(t1)>=exp[rho_min P^(3/2)(t1-t0)/2] f(t0).

All transformations back to the original physical canonical phase and
their inverses have at most polynomial momentum growth on I. The exact
constraint recovery also uses only rational/polynomial momentum factors
with nonzero pivots. Consequently the original full physical Gaussian
propagator has a lower bound c P^-m exp[c1 P^(3/2)|I|] for some finite
m and c,c1>0. No bound with a finite loss of spatial Sobolev derivatives
can hold for that propagator on these fixed prescribed backgrounds.

For a periodic spatial box the real Fourier parity modes already give
the conclusion. On R^3, restrict to a small cone of wave vectors near a
fixed direction, choose smooth transverse frames there, and use smooth
Fourier packets in bounded-width shells with |k| near P. Isotropy makes
the coefficient estimates uniform in that cone and shell. The same
exponential lower bound beats every polynomial Sobolev weight. No claim
about the special exactly homogeneous P=0 constraint is used.

These are initial-data/operator growth estimates on I, not a lower bound
for a newly computed quantum stress or a proof that the fixed prepared
quantum covariance excites a specified amount of the growing branch.
They also do not assert existence of a corresponding unforced nonlinear
solution. Both distinctions matter for the original bounce problem.
