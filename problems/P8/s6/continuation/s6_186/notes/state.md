# Actual all-order mode remainder, not a vacuum reset

Keep the unchanged all-order Borel Cauchy state, with
the exact evolved mixing relative to W8. The source-pinned
S6.85 uniform bound is

    |beta(t)|<=B6/nu^6+E/nu^10,
    B6=1813229, E=5347035781757616.

The analytic Wronskian gives |alpha|<2. The exact mixing
transport generator has coefficient magnitude bounded
by the residual times |fref|^2. With
|residual|<=C/omega^8, C=5207948696836,
|fref|^2<=1/omega and |alpha|+|beta|<3, its integral
on the unit interval gives the conservative estimate

    |alpha(t)-alpha_initial|<=6C/nu^9.

The reference readout constraint eliminates derivatives
of alpha,beta from the physical momentum. Accordingly
the actual ten-field mode is exactly

    u_actual=alpha_initial u_ref+e,

where alpha_initial is time independent and
||e||<=4000 R nu^(-11/2),
R=B6+E/m^4+6C/m^3<2e6.
The opposite-frequency readout has the same norm;
its spatial i k factors are retained. No derivative
of an oscillatory beta coefficient is estimated or
silently declared slowly varying.

Since |alpha_initial|<2, the full stress matrix bound
gives the analytic reference pair coefficient
norm below2(2*4000)^2 sqrt(nu mu)<1e9 sqrt(nu mu).

For the difference between the actual and reference
pair, expand all three terms, including e_k e_l.
Because R/m^6<1, their coefficient is bounded by

    2*4000^2 sqrt(nu mu)
      [2R(nu^-6+mu^-6)+R^2 nu^-6 mu^-6]
    <=2*4000^2*3R sqrt(nu mu)(nu^-6+mu^-6)
    <1e15 sqrt(nu mu)(nu^-6+mu^-6).

The exact initial phase and nonzero all-order mixing
are retained. This decomposition is only a proof
device for the actual state; it does not replace it
by a finite-WKB vacuum.
