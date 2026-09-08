# S6.50: actual vector lapse energy and finite curved mode error

Original P8 remains OPEN. Retain unchanged
S6.42 and all frozen certificates. Derive the vector contribution
to the physical lapse equation before setting the clock-dependent
mass coefficients to their on-clock values. In particular a_N and
b_N do not vanish, even though a=b=1 on the rolling clock.

Use the actual canonical transverse and longitudinal modes at fixed
comoving momentum on u in [-1/2,1/2]. Construct a positive fourth-order
WKB reference and exact normalized modes with matching initial data.
Bound their physical lapse-energy difference uniformly in momentum,
with an absolutely convergent momentum integral. This finite evolution
error must not be conflated with the full renormalized energy.

## Precisely certified component

The quadratic physical lapse and isotropic pressure forms retain the
temporal constraint, canonical map and fixed-comoving momentum. The
zero-derivative lapse pole matches C+C_N, not C, of S6.47. Both actual
canonical mode frequencies enter a positive fourth-order WKB reference
W. For m=m0*tau>=1000 its exact oscillator residual is bounded by
C/omega^4, C=2000000, uniformly on the full momentum half-line and
the closed time interval. Its canonical Wronskian is i.

Prepare exact Gaussian modes to match that reference and its derivative
at u=-1/2. Exact variation of constants preserves the CCR and bounds
their physical energy difference. Integrating over all R^3 momenta,
including two transverse and one longitudinal polarization, gives

    abs(Delta rho)/(M^2/tau^2)<=12C/[(m0*tau)*(M*tau)^2],
    abs(Delta p)/(M^2/tau^2)<=108C/[5*(m0*tau)*(M*tau)^2].

At M*tau=10^12,m0*tau=1000 these are 2.4e-20 and 4.32e-20.
They bound integrals of differences, not either unsubtracted integral.
The state preparation, compact interval, physical units, coefficient
envelopes and analytic momentum integral are part of the claim.
No momentum grid, stationary surrogate, arbitrary order-one constant
or claimed exponentially small production is used.

The written derivations are [energy](notes/energy.md) and
[comparison](notes/comparison.md). The report pins every local source
and fully rebuilds the frozen S6.49 ancestry. It is exact symbolic
verification plus a continuous written proof, not formal verification.

## Explicit boundaries

State any adiabatic subtraction scheme separately. A finite-order
adiabatic Gaussian comparison is not automatically an all-order
Hadamard condition, a full covariant counterterm matching, a corrected
bounce or an original V/G/B completion. These distinctions remain
research tasks; no new user choice is currently required.

In particular the exact Gaussian comparison is not asserted to be
all-order Hadamard. Neither the reference subtraction nor its full
covariant finite matching has yet been supplied. No full renormalized
energy, other field loops, all-time state tail, corrected constraints
or cones, interacting cutoff or original V/G/B result follows.
