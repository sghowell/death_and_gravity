# P8-A.20 — nonminimal conformal-scalar QSEI with cosmological calibration

This is a separately named one-real-scalar result. It does not reuse a
photon quantum inequality, photon state, or photon anomaly coefficient.
Original P8 remains open.

## Field and reference

Use four-dimensional smooth positive global flat FLRW,
M=I_s times R^3, g=ds^2-a(s)^2 dX^2, H=a'/a, with the FK curvature
convention R=6(H'+2H^2). The real free field satisfies
(Box_g+R/6)Phi=0: mass zero, nonminimal coupling xi=1/6.
The physical stress is the improved conformal scalar stress.
The actual target may be any Hadamard state; homogeneity and
quasifreeness are not required for the field inequality.

Let omega_conf be the conformal Minkowski vacuum on this actual strip.
Write w_omega=<:Phi^2:conf>_omega for its relative Wick square.
This is a difference of states, independent of an additive local
Wick-square renormalization. It is not an absolute positive operator.
The fixed stress prescription has a new finite real coefficient beta_S,
defined in notes/anomaly.md, separate from all photon and older scalar
parameters. Lambda and Newton's coupling are fixed separately.

For every real compact proper sampler, the effective energy density
E=T_UU-tr(T)/2 obeys

    integral E_omega f^2
      >= -7 hbar/(144 pi^2) ||L_H f||_2^2
         -(1/3) integral w_omega |G_H f|^2
         +integral E_conf f^2,
    L_H f=f''-2Hf'+(3H^2/4-3H'/2)f,
    G_H f=f'-3Hf/2.                                      (1)

Both terms use proper time and the physical field. The reference
retains its trace anomaly and finite beta_S term. The inequality also
extends to compact H2_0 samplers with zero value and derivative traces.
It is a quantitative nonoptimal bound, not a sharp or null QEI.

## Cosmological theorem and its exact hypotheses

Assume the actual SEE

    G_FK+Lambda g_FK=-kappa(T_scalar+T_other), kappa>0,

with E_other>=ell. Fix tau>0 and require the actual past interval
[-tau/100,0] with smooth endpoint neighborhoods. On x=s/tau in
[-1/100,0], for some p in [1/2,2/3], put

    E_history(x)=tau H(tau*x)+p/(p/2-x).

The four suprema of its derivatives through order three are at most

    (1/100,1/2,3,16).                                    (H)

The pinned geometric proof then gives H<=-19/(10tau) and the strict
past caps d=(21/10,9,70,800), in units tau^(-j-1).

The following are SEPARATE conditional extension hypotheses: IF a
comoving normal reaches s=tau, its entire segment 0<=s<tau obeys

    |H^(j)(s)|<=c_j/(tau-s)^(j+1), j=0,1,2,3,
    c=(4,128,16384,1048576),                              (F)

and the one-sided field bound below holds there as well as on the
actually existing past interval:

    w_omega(s)<=Phi_*^2, Phi_*^2>=0.                      (W)

Only the past part of (W) is unconditional. For an inhomogeneous
target, require this bound along the normal used in the argument.
It need not be a bound on |w_omega|. Neither (F) nor the future (W)
follows from (H), the field equation, or an assumed sign of an
effective Newton constant.

Set

    delta=kappa hbar/(8pi^2 tau^2),
    zeta=kappa Phi_*^2/3,
    sigma=tau^2[Lambda_+ + kappa max(-ell,0)].

If

    delta(1+|beta_S|)<=10^-8, zeta<=1/5000, sigma<=5,        (Q)

then the stipulated global spacetime has upper proper-clock endpoint
at most tau and is timelike geodesically incomplete in that direction.
Every timelike curve from s=0 has proper length at most tau.
Reversing the contracting normal gives the past-oriented statement.
This concerns the stipulated global spacetime, not every extension.

The exact scalar cost is C0+Cbeta|beta_S|, with

    C0=373108140471263/75000000 <5000000,
    Cbeta=3345309847373/10500000 <5000000,
    Cphi=1679141/10000.

The worst-case dimensionless margin is

    2009079/350000 -18/5 -1/20
      -(1679141/10000)/5000 -(1313/3500)*5
       =63325013/350000000 >9/50.                         (2)

No pointwise initial SEC is an independent premise. This particular
radiation/dust history tube nevertheless has initial timelike
convergence; it is not advertised as an SEC-violating history example.
Cosmological strength means the displayed macroscopic scaling and
dimensionless budgets, not observational verification or optimality.

## Two distinct concrete controls

A smooth finite-energy coherent family on Minkowski spacetime makes
every fixed nonzero compactly averaged E arbitrarily negative.
It proves that deleting all state-dependent qualifications is
impossible even for this conformally coupled field. It is not an
actual self-consistent SEE counterexample.

A NEW one-scalar occupation state with beta_S=Lambda=T_other=0
realizes an exact SEE low branch and the full required past history,
including the small field budget. Its thermal density coefficient
is hbar*pi^2/(30 beta_T^4), its Wick square is hbar/(12a^2 beta_T^2),
and its dimensionless anomaly coupling is lambda=delta/360.
This certifies short-past state/source compatibility. Its future
field strength grows beyond (W) near its finite earlier endpoint;
it is not a proof of that future hypothesis and not a new QEI-only
endpoint argument.

See the six source-pinned notes for the analytic and scope proofs.
The certificate is exact algebra plus written arguments, not a
proof-assistant formalization or independent peer review.
