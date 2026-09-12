# Retarded origin and the entire flat nonlocal representative

The S6.199 normalized full spatial spectral cut gives the retarded
sine kernel

    theta(t) integral K(s,p) sin(sqrt(s+q)t)/sqrt(s+q) ds,
    K=sum_i rho_i(s)Q_i(s,p)/s^2.

At a finite spectral limit this is an ordinary integral. For each
E=sqrt(s+q), y=G_E,ret Gamma solves y''+E^2 y=Gamma with common
zero initial data. If Gamma vanishes near the preparation surface,
all its initial jets vanish. Direct differentiation proves

    y=Gamma/E^2-Gamma''/E^4+Gamma''''/E^6
      -G_E,ret Gamma^(6)/E^6.

The remainder transforms to w^3/[E^6(E^2-w)]. Thus after the
full physical tensor cut has been reconstructed, its sixth-derivative
bulk is exactly B6. This does not drop odd endpoints from S6.198's
curved complex-amplitude formula. There the amplitude and phase vary
with time, so this flat sine-kernel simplification is unavailable.

Using |sin|<=1, the full matrix bulk majorant is

    integral (2rho2+rho0)(1+q/s)^2/(s+q)^(7/2) ds
      <=7/(384pi^2) integral_(4m^2)^infinity (s+q)^(-3/2) ds
      =7/(192pi^2 sqrt(4m^2+q)).

Here the complete tensor norm includes both the congruence and trace
terms. The density bound uses 2/128+1/384=7/384. The result is
uniform in external p and absolutely integrable in internal spectral
mass. For compact real time tests on the unit slab, time L1 is bounded
by L2, and Plancherel plus Cauchy-Schwarz gives the bulk bound
1e-5||D||L2||Gamma^(6)||L2 at m=1000.

For a tail beginning at Lambda the same calculation gives
7/(192pi^2 sqrt(Lambda+q)), bounded by 1/(100sqrt(Lambda)).
No time-frequency restriction or spatial momentum cutoff is imposed.

Combining the exact conversion and these bounds defines the entire
specified F as a retarded distribution on the stated test space.
For this flat representative, finite covariant spectral integrals are
polynomial differential numerators applied to retarded massive
Klein-Gordon kernels. They have causal support; the proved distributional
limit preserves it. The spatially nonlocal instantaneous conversion
coefficients cannot individually be discarded on causality grounds.

These facts concern the specified nonlocal representative. Its remaining
physical finite local polynomial and the actual curved CD matching
are not supplied by this argument.
