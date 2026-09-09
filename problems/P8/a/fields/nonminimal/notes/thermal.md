# Actual one-scalar thermal state and its short SEE history

## Construct the state and both of its physical moments

Take ONE real massless flat oscillator with field amplitude
sqrt(hbar/(2k)) and occupation

    n(k)=1/[exp(beta_T k)-1], beta_T>0,
    <a_dagger a>=n, <a a_dagger>=1+n, <aa>=0.

Nonnegative occupation gives positivity, and the difference of the
two covariances is one, preserving the vacuum commutator.
The thermal-minus-vacuum two-point function is smooth: the
large-k occupation is exponentially decreasing, while at k=0
n(k)=O(k^-1) and radial measure k^2 times field amplitude squared
k^-1 leaves an integrable integrand. Additional spacetime
derivatives add only nonnegative powers of k. Every derivative
is bounded on compacts by an integrable function.

The conformal map Phi=a^-1 phi transports this actual state onto
the actual smooth positive flat-FLRW strip, preserving Hadamard,
positivity and commutators. It is not an instantaneous thermal
state chosen separately at each proper time.

Expanding the positive Bose factor and applying Tonelli gives

    Q=hbar/(2pi^2) integral k^3 n(k) dk
       =hbar*pi^2/(30 beta_T^4),
    w_flat=hbar/(2pi^2) integral k n(k) dk
       =hbar/(12 beta_T^2).

Here Gamma(4)zeta(4)=pi^4/15 and Gamma(2)zeta(2)=pi^2/6.
Both moments belong to the SAME state, and satisfy

    Q=(24pi^2/(5hbar)) w_flat^2,
    rho_thermal=Q/a^4, w_physical=w_flat/a^2.

There are no two photon polarizations. The physical temperature is
k_B T_physical=hbar/(a beta_T). beta_T is an inverse conformal
frequency, not beta_S or a renormalization parameter.

## Fresh scalar source and independent SEE components

For this example ONLY set beta_S=0, Lambda=0, T_other=0, with no
independent nonzero gravitational R^2 coupling. Put
r_S=hbar/(960pi^2). The actual source is

    rho=Q/a^4+r_S H^4,
    pressure=Q/(3a^4)-r_S H^4-(4r_S/3)H^2 H',
    E=Q/a^4-r_S(H^4+2H^2 H').

Conservation is exact; the anomaly is retained. The actual equations
3H^2=kappa rho and 2H'+3H^2=-kappa pressure reduce to

    H^2-bH^4=kappa Q/(3a^4),
    H'=-2H^2(1-bH^2)/(1-2bH^2),
    b=kappa hbar/(2880pi^2).                              (1)

thermal.py independently checks both components, not just the
density constraint. With x=s/tau,y=-tau H and delta as in the
formulation, the new scalar dictionary is

    lambda=b/tau^2=delta/360,
    y'=2y^2(1-lambda y^2)/(1-2lambda y^2)=F(y).             (2)

At a(0)=1,y(0)=2 the constraint fixes
Q=12(1-4lambda)/(kappa tau^2). Choose the actual state parameter
beta_T=[hbar*pi^2/(30Q)]^(1/4). The state therefore realizes the
source; Q is not a separately assumed classical radiation fluid.

On the low branch 0<y<1/sqrt(2lambda),

    a^4=4(1-4lambda)/[y^2(1-lambda y^2)],
    x=Phi(2)-Phi(y),
    Phi(y)=1/(2y)+sqrt(lambda)atanh(sqrt(lambda)y)/2.

Direct differentiation gives d_x log a=-y and Phi' F=-1.
These are exact, not order-reduced or asymptotic equations.
The pin to A.19 is used only for its GENERIC rational-ODE
comparison and endpoint lemmas, after deriving the scalar
equation and verifying its scalar lambda domain afresh.

## Actual past geometry AND state strength

For lambda<=1/16, the analytic low-branch ODE has on
x in [-1/100,0] the invariant backward box 50/27<=y<=2.
Indeed F>0,F<=4y^2; applying these to v=1/y prevents exit
from 1<=y<=2 and proves continuous existence on a neighborhood
of the closed sample interval.

Relative to the anchored radiation history y_rad=2/(1-4x),

    v'-v_rad'=-2lambda y^2/(1-2lambda y^2),
    0<=y_rad-y<=(16/25)lambda.

The exact successive derivative corrections at a common y are

    F-2y^2=2lambda y^4/D,
    F F_y-8y^3=8lambda y^5(6z^2-8z+3)/D^3,
    F(F F_y)_y-48y^4
      =16lambda y^6(6z^2-5z+2)(14z^2-22z+9)/D^5,

where z=lambda y^2<=1/4,D=1-2z>=1/2. The quadratics are
positive decreasing on that interval, bounded by 3,2,9.
Together with the radiation-polynomial Lipschitz constants
8,96,1536, this gives all-x C3 errors at most

    lambda (16/25,1728/25,155136/25,14770176/25).

For delta<=10^-8, lambda<=1/36000000000; every component is
strictly inside (1/100,1/2,3,16). This proves actual membership
in the anchored p=1/2 slice, not in every p slice.

On this past a>=a(0)=1. The actual field budget is therefore
bounded by its initial value zeta0=kappa w_flat/3, with

    zeta0^2=(20/9)delta(1-4lambda) < (1/5000)^2.

This identity follows from the two moments of the SAME state
and its constraint Q. It is not an independently adjustable
field-amplitude allowance. beta_S=0 and sigma=0 also meet
the named past/parameter gates.

## What this actual example does NOT verify

The ODE covers the open strip -infinity<x<x_end, where
x_end=Phi(2)-Phi(1/sqrt(2lambda)) lies strictly between 1/16
and 1/4. The bound follows by integrating 1/F from y=2:
it is below 1/(2y^2), and on [2,4] above 1/(4y^2) because
lambda<1/64. At the endpoint a_end^4=16lambda(1-4lambda)>0
for each fixed positive delta, while y' diverges.
Every field and SEE assertion is on the smooth open strip,
not on its excluded endpoint.

The actual Wick-square parameter tends instead to

    zeta_endpoint^2=50.

Thus the small field cap cannot be claimed on every shorter
future piece up to this endpoint. The cosmological theorem's
future hypothesis is deliberately conditional on reaching tau;
the example already ends earlier. It supplies genuine short-past
compatibility, not a nonvacuous verification of that future cap.

Moreover kappa tau^2 E=3(F-y^2)=3y^2/(1-2lambda y^2)>0.
This state already has timelike convergence. Its endpoint is
not a new QEI-only argument, and its high-curvature behavior
does not certify fundamental quantum-gravity or EFT validity.
