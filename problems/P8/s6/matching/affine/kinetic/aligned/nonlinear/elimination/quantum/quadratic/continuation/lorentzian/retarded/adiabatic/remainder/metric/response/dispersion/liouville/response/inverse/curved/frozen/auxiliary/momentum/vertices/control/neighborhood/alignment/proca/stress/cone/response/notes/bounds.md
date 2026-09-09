# Continuous local response bounds

For each of the physical energy/pressure and the clock lapse/hat-log-scale Euler outputs, write the actual source-linear response at half-adiabatic order n as a sum of rational functions of u multiplying the two source rows through derivative four. Reject nonlinear sources, unknown coefficients and any derivative above four.

Each rational coefficient is bounded on the entire I=[-1/2,1/2] using the unchanged exact polynomial interval-envelope routine. All 120 coefficient reconstructions are independently retained as zero residuals. The bound is the sum of absolute coefficient envelopes times the maximum of both source functions and all their coordinate-time derivatives through four. It is not a sampling norm.

With L=M tau, m tau=1000 and pi²>9, the normalized component bound is

    sum_(n=0)^2 m^(4-2n) row_envelope_n / (576 L²).

The normalization is physical stress or coordinate Euler response divided by M²/tau². At L=10^400 all four component bounds are strictly below 10^-790. The respective upper bounds multiplied by 10^790 are approximately 2.084e-5 (physical energy), 3.126e-5 (physical pressure), 1.758e-4 (clock lapse) and 1.582e-4 (clock scale); the report retains the exact rational bounds, and the decimal displays are rounded upward.

The clock-chart bound includes only the local matched piece of the fixed profile, by the explicit linear decomposition in the chart note. Its order-zero bound is exactly zero. The remaining state-dependent profile contribution and nonlocal kernel are not estimated here. The L^-2 scaling of these fixed-m coordinate bounds is exact, but it is not a bound on canonical interacting higher vertices or an all-frequency inverse.

