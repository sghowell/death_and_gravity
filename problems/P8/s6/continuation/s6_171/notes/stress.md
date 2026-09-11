# Complete absolute specified free curved tensor

Use the same m=1e200, Delta<=3e197, tau=1e-100 and m0=.99m.
All42 free Dirac copies are in either the exact S6.170 curved in or
out Hadamard state on a=(1+t^2)^2. Six copies have paired varying
masses; the remaining36 are constant but not geometrically inert.

The complete state/projector remainders after the local subtractions
are bounded, uniformly in real time, by

    R_rho<1e410,  R_P<1e416.

The exact rational sums retain the S6.170 half-line state allowance,
twenty-to-first-frame diagonal part, first-frame Taylor remainder,
and both linear and cubic off-diagonal pressure contributions.
Numerically R_rho<4.712e408 and R_P<6.702e414; rounded diagnostics
are not the exact inequality inputs.

Add the referenced local two-derivative terms (reference.md), the
finite Euler terms (curvature.md), and ONLY the three complete
potential components already bounded in S6.168: the fixed mass
anchor, full quartic, and the convergent entire higher-even tail.
The old flat derivative/projector allowances are not counted again.
Potential pressure equals minus potential energy by covariance.

The resulting exact positive sums in stress.py give

    |rho(t)|<1e789,  |P(t)|<1e789,
    max(|rho|,|P|)/kappa<1e-11,  kappa=1e800.

Isotropy makes the specified comoving orthonormal tensor diagonal.
These are all its components in that frame, not uniform bounds over
arbitrary boosts. The potential dominates the numerical allowance.
The ratio to kappa is a named reference comparison, not relative
error against the bounce density: the latter vanishes at H=0.

This is absolute stress of the specified quadratic/reference
one-loop sector. It includes vacuum polarization as well as state
coherence. It is not the full interacting canonical parent tensor,
a cutoff or higher-loop bound, a self-consistent bounce, or control
of metric/scalar response through the short transition.
