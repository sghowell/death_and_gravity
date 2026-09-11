# Exact clock stress transfer with explicit physical units and interval

The literal new canonical action, arbitrary-lapse/scale readouts, metric,
mass1000, phase measure and selected all-order state agree with the
ordinary S6.82 vector calculation. Its previous scalar profiles are not
part of this transfer. The fixed local prescription is the ordinary
one specified in renormalization.md, not a different reference-state
normal ordering.

For each polarization let p=v'-d v. The physical energy is
(|p|^2+omega^2|v|^2)/(2a^3), with the pressure readouts in bridge.md.
The finite stress has the full radial measure k^2 dk/(2pi^2), two
transverse polarizations and one longitudinal polarization. Subtract
the actual fourth-order adiabatic readout, keep the complete exact
state and add the matched local coefficients. No momentum cutoff
or grid replaces this continuum definition.

At mu=m, the local energy/pressure pairs multiplying
m^(4-2n)/(64pi^2), n=0,1,2, are

    rho0=-5/2, P0=5/2,
    rho1=-10H^2, P1=(10/3)(2H'+3H^2),
    rho2=12H^2H'+4HH''-2H'^2,
    P2=-12H^2H'-8HH''-6H'^2-(4/3)H'''.

The code compares every pair to the frozen dimensional calculation
and checks rho_n'+3H(rho_n+P_n)=0 without a singular division by H.
The mode and subtraction Ward identities give conservation of the
whole specified ordinary stress on the clock. Vanishing mass/source
insertions do not erase this nonzero gravitational source.

The S6.82 C5 proof uses an eighth-order reference as a proof device,
not a new state. Its exact-ODE-projected derivatives retain all three
quadratic mode products. All nonintegrable subtraction coefficients
cancel; their remaining rational envelopes are integrated over all
momenta. The initial-state mixing and the low, middle, high preparation
bands plus oscillatory-evolution term are retained. Those endpoints
are proof partitions, not physical cutoffs. Uniform integrable
derivative majorants justify all differentiations through order five.

Here tau=1, M*tau=10^400 and kappa=(M*tau)^2=10^800.
The original derivative normalization is M^2/tau^(2+j).
Consequently each stored normalized coefficient is multiplied by
kappa, not by kappa*zeta, to obtain the canonical physical stress.
Every one of the twelve derivatives satisfies, uniformly on
I=[-1/2,1/2],

    |d^j rho/dt^j|, |d^j P/dt^j| <10^30,
    |d^j rho/dt^j|/kappa, |d^j P/dt^j|/kappa <10^-770,
    0<=j<=5.

The exact rational decompositions are retained in the report.
This is an ABSOLUTE specified-vector bound. The reference denominator
is fixed kappa; it is not the zero classical bounce density. It
bounds this contribution to the prescribed-clock metric residual,
not an exact or nearby self-consistent solution. Outside I, qualitative
Hadamard smoothness on compact intervals supplies no numerical
extension of these constants.
