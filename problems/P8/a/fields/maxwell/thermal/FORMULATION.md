# P8-A.19 — actual thermal-photon SEE history in the A.18 radiation tube

This optional strengthening constructs an actual state/metric realization
of the previously certified short-history hypotheses. It does not reopen
A.18's scoped incompleteness objective or replace its all-Hadamard theorem
by a thermal-state assumption. All frozen ancestor files remain unchanged.

## Named physical specialization and actual state

Use one physical free Maxwell field and the A.16 FK metric/curvature
conventions. Explicitly specialize to beta_M=0, Lambda=0, no additional
matter source and no independent nonzero gravitational curvature-squared
coupling. The Newton coefficient kappa and hbar are positive. This is the
named zero-type-D member of A.16's photon family, not scalar gamma=0 and
not an assertion of actual SEE existence for arbitrary beta_M.

On the flat conformal comparison, construct the physical two-transverse-
polarization quasifree occupation state with

    n(k)=1/(exp(b_T k)-1), b_T>0.

Here b_T is a length in the conformal-time clock, distinct from beta_M.
The physical temperature obeys k_B T_physical=hbar/(a b_T), in c=1 units.
Conformal transport of the gauge-invariant field strength gives a positive
Hadamard state on every smooth positive flat-FLRW strip. Only the physical
polarizations are used; no gauge-potential index/weight formula is assumed.

The exact thermal coefficient and total physical stress are

    Q=hbar*pi^2/(15 b_T^4), r=31 hbar/(480 pi^2),
    rho=Q/a^4+rH^4,
    pressure=Q/(3a^4)-rH^4-(4r/3)H^2 Hdot,
    EED=Q/a^4-r(H^4+2H^2 Hdot).

The thermal charge Q is determined by this actual state, not a radiation
integration constant assumed without a state realization. The anomaly
term rH^4 is not omitted. This is a conformally thermal state, not global
equilibrium under cosmic-time translations of a time-dependent metric.

## Exact semiclassical branch and clock

Fix tau>0 and

    0<delta=kappa*hbar/(8pi^2 tau^2)<=10^-8,
    b=kappa r/3, lambda=b/tau^2=31delta/180.

Let s be the contracting-normal proper parameter, x=s/tau and
y(x)=-tau H_normal(s)>0. Impose a(0)=1, y(0)=2, and choose the state by

    Q=12(1-4lambda)/(kappa tau^2)>0,
    b_T=[hbar*pi^2/(15Q)]^(1/4).

The scale normalization a(0)=1 is explicitly chosen. It is not inferred
from a Hubble-rate bound. In the physical expansion orientation s=T_0-T.
Both independent SEE components are equivalent on the nonzero-H branch to

    H^2-bH^4=kappa Q/(3a^4),
    Hdot=-2H^2(1-bH^2)/(1-2bH^2).

Consequently the exact dimensionless normal equation and first integral are

    y'=F_lambda(y)=2y^2(1-lambda y^2)/(1-2lambda y^2),
    a^4=4(1-4lambda)/[y^2(1-lambda y^2)].

Use only the low branch 0<y<ycrit=1/sqrt(2lambda). A positive algebraic a^4
on the high branch is not sufficient for admission. The smooth low-branch
solution occupies the open interval -infinity<x<x_end, where

    Phi(y)=1/(2y)+(sqrt(lambda)/2)atanh(sqrt(lambda)y),
    x=Phi(2)-Phi(y),
    x_end=Phi(2)-Phi(ycrit), 1/16<x_end<1/4.

These are exact, not order-reduced or leading-radiation approximations of
the named model. The source's finite prescription has already been fixed;
no derivative term is subsequently deleted from the chosen equation.

## Continuous actual radiation-history embedding

The whole required actual interval x in [-1/100,0] is smooth. For
y_rad=2/(1-4x), the continuous derivative bounds through order three are

    ||d_x^j(y-y_rad)|| <=lambda C_j,
    C=(16/25,1728/25,155136/25,14770176/25).

At delta=10^-8 these four upper bounds are exactly

    (31/28125000000, 93/781250000,
     3131/292968750, 149048/146484375),

strictly below A.18's (1/100,1/2,3,16) allowed jet errors. Smaller positive
delta improves them uniformly. Thus the actual thermal state and exact SEE
history lie inside the ANCHORED p=1/2 radiation slice of the A.18 C3 tube,
with actual |H_normal(0)| tau=2 and zero additional-source budget.
This is a bound on every point of the history, not a finite sample scan or
an unspecified small-error term.

## Endpoint and nonclaims

At the finite branch endpoint,

    a_end^4=16lambda(1-4lambda)>0,
    tau^2 R_FK=-12lambda y^4/(1-2lambda y^2) -> -infinity.

The positive endpoint scale holds for each fixed delta>0; no uniform
positive lower limit as delta tends to zero is asserted. The endpoint is
excluded from the smooth spacetime, not silently treated as regular
initial data for another branch. The actual maximal low-branch spacetime
has a finite contracting-direction comoving proper length tau*x_end,
and is timelike geodesically incomplete in that direction. In the
expansion orientation this is its past endpoint. No larger extension
is manufactured to evaluate A.18's future samplers through this endpoint.

This particular thermal branch has positive EED and satisfies timelike
convergence. Its endpoint follows directly from its specified exact SEE;
it is not a new QEI-only argument, all-state SEE theorem or claim of
fundamental EFT validity at anomaly-scale curvature. A.18's all-Hadamard
conditional theorem remains independent of this example. The present
result supplies actual short-history/state/source compatibility as an
optional concrete strengthening, not an extra closure prerequisite.
