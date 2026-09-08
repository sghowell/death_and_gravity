# S6.54: C1 vector matching and homogeneous clock source

Status: the unchanged retained Gaussian vector's finite metric
first variations have a quantified first time derivative and
homogeneous clock source on the rolling interval. Original P8
remains OPEN; no state, frozen action or matching scheme changes.

Keep S6.50's exact prepared modes, S6.51's full fourth-order
subtraction and S6.53's scalar-coefficient continuation. Work on
I=[-1/2,1/2], all real comoving momenta and m0*tau>=1000.

## Precisely scoped results

The [differentiated integral proof](notes/derivative.md) bounds
derivatives of the actual rational WKB coefficients. With
C=2,000,000 and C1=100,000,000, the exact reference residual
and its derivative are bounded by C/omega^4 and C1/omega^4.
An exact diagonal-phase transformation followed by one oscillatory
integration by parts, including both endpoints and the feedback
term, improves the prepared mixing estimate to

    |Bcal|<=K/nu^6, K=676,000,000,
    nu^2=m^2+k^2/(25/16)^2.

The derivative of the physical quadratic readout is computed
before momentum integration. Its exact/reference difference is
absolutely integrable and bounded by 40K/m. The differentiated
full reference subtraction tail adds at most D1/(72m^2),
D1=100,000,000. These continuous envelopes give a C1 finite
state energy and pressure term, not only a pointwise value.

The [clock-source proof](notes/clock-source.md) derives the
homogeneous clock variation from the original scalar mass functions
and the same counterterm action. It checks ordinary conservation
at the exact mode, full adiabatic-coefficient and matched local
levels. The actual mass energy instead exchanges with the clock:

    J_clock=-[rho'+3H(rho+p)]
           =-[rho_mass'+3H rho_mass].

Compact variations are supported inside I with initial Gaussian
data fixed. The S6.53 regulator-limit result also controls the
distributional clock-source limit; the present C1 bound identifies
it with the stated continuous function.

Continuous actual-clock local envelopes complete the estimate.
At M*tau=10^12,m0*tau=1000, both density derivatives and the
vector clock source are below 10^-14 of M^2/tau^3. The physical
clock is phi=tau*u. Every local, state-evolution and reference-tail
error term is explicit; no unspecified order-one constant is used.

## Boundaries

This is a vector-only one-loop first-variation result at the fixed
classical background. It does not solve the quantum-corrected light
equations or establish their constrained principal symbol, cone,
nonlinear stability or interacting cutoff. Small unsigned source
values cannot protect a classically saturated matter cone.

Finite C1 quantities do not establish second time derivatives,
all-order Hadamard admissibility, higher quantum functional
variations, all-time control, other fields/higher loops, unknown
UV matching operators, extension to the vacuum, finite-gravity
Regge bounds or V/G/B. Original P8 is not finished or closed.
The remaining estimates are research, not a user-choice blocker.
