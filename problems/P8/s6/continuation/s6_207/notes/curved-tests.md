# Actual curved coefficient and time-jet tests

The independent numerical implementation evaluates the actual CD scale factor, both mode species, all four W8 coefficients and their total time derivatives. It uses the same prepared-state comparison conventions and both analytic Schwarz signs.

The first route builds the ten physical constrained readouts and contracts the complete stress matrices for all nine pairs. The second evaluates the production four-sector projector/scalar expressions. Noncollinear momenta, complex inverse radius, independent source and detector times, multiple compact times, and noncommuting source/detector tensors are used. Separate real-radius checks compare the normalized implementation directly to the original unscaled W8 modes.

Source-time Cauchy extraction retains the detector at fixed time 1/4 while sampling the source in a radius 10^-5 circle. All five endpoint orders and every source time jet are reconstructed through the full integration-by-parts recurrence. The two independent routes agree.

The retained inverse-radius coefficients are then extracted on a radius 10^-6 circle at 125-digit working precision, using 24 and 32 points. Every endpoint/power slot and all thirty-five source time-jet entries are compared. A separate calculation using only the first two W8 coefficients agrees on those retained jets; finite-momentum calculations continue to use all four coefficients.

For each endpoint, subtracting the retained coefficients from a fresh actual evaluation leaves a nonzero remainder consistent with the stated pre-current k^-4 order. The test compares normalized residuals at two still smaller inverse radii. This is a fixed-P check, not a replacement for the written analytic remainder argument.

A private test initially demanded nonzero ultraviolet coefficients for both odd endpoints. The actual calculation showed that this was too strong: the first odd endpoint has a nonzero UV coefficient, while the third's retained UV coefficients cancel in this fixture. Its full finite-momentum endpoint is nevertheless nonzero. The test was corrected to check both facts before freezing; neither endpoint was removed from the algorithm, and no general cancellation theorem is inferred from this fixture.

Controls also show that differentiating after detector/source coincidence adds an unwanted detector derivative, and that merging the distinct T/L inverse phases changes the actual curved result.
