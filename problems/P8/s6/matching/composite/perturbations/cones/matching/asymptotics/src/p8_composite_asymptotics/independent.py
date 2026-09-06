"""Separate coefficientwise Fraction checks, without the main symbolic engine."""

from fractions import Fraction as Q

from p8_composite_modes.independent import Poly


def identities():
    v, a, z = (Poly.variable(name) for name in ("v", "accel", "z"))
    vp = 1+a-v**2
    tr_numerator = 3*v**2-1-a
    tr_numerator_prime = 6*v*vp+a*v
    nr_numerator = 2*v*tr_numerator_prime-2*tr_numerator*vp+tr_numerator**2
    expected_nr = 3*(a+1)**2+3*v**4-2*v**2
    mass = v**2-v
    ch, sh = Poly.variable("A"), Poly.variable("B")
    scale, scale_prime = 12*ch+2*sh-11, 12*sh+2*ch
    e, k = Poly.variable("x"), Poly.variable("c")
    return {
        "full_relative_normalization_derivative": nr_numerator-expected_nr,
        "general_mass_led_positive_margin": nr_numerator-3*v**2*mass
             -3*(a+1)**2-v**2*(3*v-2),
        "uniform_heavy_quarter_margin": nr_numerator-4*v**2*mass-v**2
             -3*(a+1)**2-v**2*(v-1)*(3-v),
        "scale_first_integral_mod_hyperbolic_identity": scale_prime**2-scale**2-22*scale+19
             +140*(ch**2-sh**2-1),
        "v_upper_bound_square": 140*z**2-19*(z**2+22*z-19)-(11*z-19)**2,
        "source_projection_cleared": k*e*(1+e)-(1+k*e**2)-(k*e-1),
    }


def checks():
    coefficients = identities()
    if any(not value.is_zero() for value in coefficients.values()):
        raise ValueError("An independent asymptotic polynomial identity failed")
    v, a = Q(2), Q(11)
    vp, ap = 1+a-v*v, -a*v
    vpp = ap-2*v*vp
    mass = v*v-v
    nr = (3*(a+1)**2+3*v**4-2*v*v)/(4*v*v)
    heavy = nr-mass
    projected_first_coefficient = -vp/(2*v)
    projected_second_coefficient = -vpp/(2*v)+vp*vp/(4*v*v)-heavy
    values = {"initial_vprime": vp, "initial_vsecond": vpp,
              "initial_mass_squared": mass, "initial_N_relative": nr,
              "initial_heavy_growth_coefficient": heavy,
              "projected_first_coefficient_of_sqrt2": projected_first_coefficient,
              "projected_second_coefficient_of_sqrt2": projected_second_coefficient,
              "strict_v_squared_upper_below_nine": 9-Q(140, 19),
              "backward_z_min_is_positive_square_margin": Q(140)-11**2,
              "backward_z_min_is_below_one_square_margin": 12**2-Q(140)}
    expected = (8, -54, 2, Q(59, 2), Q(55, 2), -2, -10, Q(31, 19), 19, 4)
    if tuple(values.values()) != expected:
        raise ValueError("An independent exact limiting fixture failed")
    return {"coefficientwise_identities": dict.fromkeys(coefficients, "0"),
            "exact_forward_growth_and_domain_fixtures": {key: str(value) for key, value in values.items()}}
