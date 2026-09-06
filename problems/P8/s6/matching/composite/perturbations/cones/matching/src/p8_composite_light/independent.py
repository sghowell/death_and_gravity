"""Separate Fraction coefficient and exact hierarchy arithmetic."""

from fractions import Fraction

from p8_composite_modes.independent import Poly


def identities():
    ks, w, ud, dd, wd, delta = (Poly.variable(name) for name in ("A", "B", "x", "z", "u", "v"))
    hd, gd = ud-w*dd-wd*delta, ud+(1-w)*dd-wd*delta
    weighted = ks*(1-w)*hd**2+ks*w*gd**2-ks*(ud-wd*delta)**2-ks*w*(1-w)*dd**2
    f, g, fd, gd, omega, omega_d, theta, cross = (Poly.variable(name) for name in
                                                ("alpha", "beta", "x", "z", "c", "y", "p", "rho"))
    bg = 2*omega*(gd-theta*g)+cross*g
    badjf = -2*omega*fd+(cross-2*omega_d-2*omega*theta)*f
    adjoint = f*bg-badjf*g-2*(omega_d*f*g+omega*fd*g+omega*f*gd)
    # Independent exact symmetric second-jet difference and completion.
    a = Poly.variable("accel")
    n = (144*a**2-6*a+67)*Fraction(1, 84)
    routh = Fraction(9, 4)-Fraction(9, 2)*a
    return {"weighted_kinetic_and_moving_ratio": weighted,
            "formal_adjoint_full_time_surface": adjoint,
            "symmetric_Routh_hierarchy_difference": 3*n-routh-(36*a**2+30*a+1)*Fraction(1, 7),
            "symmetric_N_positive_completion": n-Fraction(12, 7)*(a-Fraction(1, 48))**2-Fraction(51, 64)}


def checks():
    values = identities()
    if any(not value.is_zero() for value in values.values()):
        raise ValueError("An independent light-mode coefficient identity failed")
    q = Fraction
    lower, upper = q(10511, 500), q(106, 5)
    small_a_numerator = q(1037)-(376+16*q(22))/4
    small_a_denominator = 16*(60996+3047*q(22))
    margins = {"sqrt442_lower": 442-lower**2, "sqrt442_upper": upper**2-442,
               "small_A_omega_lower_above_third": small_a_numerator**2/small_a_denominator-q(1, 3),
               "mass_upper": (107-5*lower)/42,
               "mass_led_ratio": q(1, 4)/q(9, 200),
               "symmetric_positive_Routh_at_4_over_25": q(9, 4)-q(9, 2)*q(4, 25)}
    if not (all(value > 0 for value in margins.values())
            and margins["mass_upper"] == q(9, 200)
            and margins["mass_led_ratio"] == q(50, 9)
            and margins["symmetric_positive_Routh_at_4_over_25"] == q(153, 100)):
        raise ValueError("An independent exact hierarchy control failed")
    return {"coefficientwise_identities": dict.fromkeys(values, "0"),
            "exact_hierarchy_and_countercontrol": {key: str(value) for key, value in margins.items()}}
