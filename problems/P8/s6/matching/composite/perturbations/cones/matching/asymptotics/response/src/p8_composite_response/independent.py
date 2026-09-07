"""Separate exact coefficient tests for the flux Schur restriction and pulse."""

from fractions import Fraction as Q

from p8_composite_modes.independent import Poly


def checks():
    # The immutable primitive has a fixed alphabet; these distinct names
    # are a ring isomorphism, not an import of its physics conventions.
    eps, g, f, vg, vf, w, p, z, e = (Poly.variable(key) for key in
                                      ("x", "A", "B", "alpha", "beta", "C", "p", "z", "y"))
    # Clearing the positive source denominator (1+e) gives independent
    # polynomial derivatives of the literal metric potential/source action.
    force_g = (-(vg+w)*g+eps*w*f)*(1+e)+p*z
    force_f = (eps*w*g-(vf+eps**2*w)*f)*(1+e)+p
    common = eps*force_g+force_f
    expected = -(eps*vg*g+vf*f)*(1+e)+(1+eps*z)*p
    identities = {
        "locked_potential_cancels_relative_mass": common-expected,
        "source_factor_is_physical_composite_weight": (1+eps*z)*p-(1+e)*p-(eps*z-e)*p,
        "relative_mass_energy_square": w*(g-eps*f)**2-w*g**2+2*eps*w*g*f-eps**2*w*f**2,
    }
    if any(not value.is_zero() for value in identities.values()):
        raise ValueError("An independent flux polynomial identity failed")
    # Integrate u²-2u³+u⁴ and (1-u) times it coefficientwise.
    pulse = {2: Q(1), 3: Q(-2), 4: Q(1)}
    moment = sum((value*(Q(1, degree+1)-Q(1, degree+2))
                  for degree, value in pulse.items()), Q(0))
    integral = sum((value/Q(degree+1) for degree, value in pulse.items()), Q(0))
    margins = {"pulse_integral": integral, "pulse_first_moment": moment,
               "heavy_band_coefficient": Q(1, 4)-Q(140, 19)/64,
               "finite_absolute_mismatch": moment-Q(2, 600),
               "finite_relative_error_margin": moment/3-Q(5, 1800)}
    if tuple(margins.values()) != (Q(1, 30), Q(1, 60), Q(41, 304), Q(1, 75), Q(1, 360)):
        raise ValueError("Independent pulse/band/mismatch arithmetic failed")
    return {"coefficientwise_identities": dict.fromkeys(identities, "0"),
            "Fraction_comparison_margins": {key: str(value) for key, value in margins.items()}}
