# Positive Minkowski covariances and Cartesian polarization

In this note the displayed covariances exclude only the common
positive factor hbar/kappa. Put omega=sqrt(m^2+q). For Proca define

    Q=(m^2 I+kk^T)/omega,
    P=(omega/m^2)(I-kk^T/omega^2),
    C_P=[[Q,iI],[-iI,P]]/2.

Let R=(m I+kk^T/(omega+m))/sqrt(omega). Its inverse is
sqrt(omega)[I-kk^T/(omega(omega+m))]/m. Both are regular at k=0;
R^2=Q and R^-2=P. The 6 by 3 matrix G=(R;-i R^-1)/sqrt(2)
obeys C_P=G G*, C_P-C_P^T=i Omega6, and M_P C_P+C_P M_P^T=0
at a=1. These identities fix positivity, commutator and time
orientation of the chosen flat state without an angular frame.

Reconstruct (W0,W_i) with W0=-i k.Pi. With eta=diag(-1,1,1,1)
and xi=(-omega,k), its covariance is

    m^2 [eta+xi xi^T/m^2]/(2 omega)
      =m^2 epsilon epsilon^T/(2 omega),
    epsilon=(-k^T/m; I+kk^T/[m(omega+m)]).

Direct checks give xi^T eta epsilon=0 and epsilon^T eta epsilon=I3.
This positive three-column factor is the massive physical
polarization sum. The metric eta alone would have a negative time
entry and is not the state covariance.

For one tensor polarization at a=1,

    C_T=[[1/k,i/2],[-i/2,k/4]]
       =(1/sqrt(k),-i sqrt(k)/2)
        (1/sqrt(k),-i sqrt(k)/2)*,  k>0.

It is stationary in the correctly normalized M_T and has CCR iJ2.
Its k^-1 singularity is infrared, not a missing oscillator or
permission to introduce a zero-frequency atom.

For q>0 define PL=kk^T/q, PT=I-PL, and

    Pi_ij,kl=(PT_ik PT_jl+PT_il PT_jk-PT_ij PT_kl)/2.

PL and PT have ranks one and two. The Cartesian TT matrix Pi
is a self-adjoint idempotent, trace-free and transverse, of rank two.
It is bounded uniformly in the momentum direction. In an aligned
frame Pi=(E1 E1^T+E2 E2^T)/2, with E_lambda:E_mu=2 delta_lambda,mu.

One global-frame-free realization takes nine independent minimal
scalar fields Phi_ij with momenta P_ij and sets
t=2 Pi Phi, pi_TT=Pi P/2. Then [t,pi_TT]=i Pi and the covariance
has only the two physical TT configuration directions. At each
nonzero k this gives the same two C_T blocks. The duplicated
nine-scalar realization is a convenient representation of the
projected algebra, not nine additional physical modes.
