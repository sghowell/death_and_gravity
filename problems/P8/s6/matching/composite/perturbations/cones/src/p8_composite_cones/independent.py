"""Coefficientwise Fraction audit, reusing the pinned independent polynomial ring."""

from fractions import Fraction

from p8_composite_modes.independent import Poly


def identities():
    alpha, beta, y, c, pressure, scale = (Poly.variable(name) for name in
                                        ("alpha", "beta", "y", "c", "p", "m4"))
    b1, b2, b3 = (Poly.variable(f"beta{j}") for j in (1, 2, 3))
    r, s = alpha+beta*y, alpha+beta*c
    potential = scale*(b1+2*b2*y+b3*y**2)
    branch = potential-alpha*beta*r**2*pressure
    mu = y*(scale*(b1+b2*(y+c)+b3*c*y)-alpha*beta*r*s*pressure)
    shear = scale*(b2+b3*y)-alpha*beta**2*r*pressure
    return {"g_speed_squared_cleared": r**2-s**2-beta*(y-c)*(r+s),
            "f_speed_squared_cleared": (c*r)**2-(y*s)**2-alpha*(c-y)*(c*r+y*s),
            "tensor_velocity_convex_sum_cleared": alpha+beta*c-s,
            "positive_weights_sum_cleared": alpha+beta*y-r,
            "pressure_stiffness_coefficient_identity": mu-y*branch-y*(c-y)*shear,
            "squared_speed_average_gap_cleared":
                alpha*y*r+beta*c**2*r-y*s**2-alpha*beta*(c-y)**2}


def checks():
    checks = identities()
    if any(not value.is_zero() for value in checks.values()):
        raise ValueError("An independent cone coefficient identity failed")
    a, b, y, c = Fraction(2), Fraction(3), Fraction(4), Fraction(5)
    r, s = a+b*y, a+b*c
    g, f = r/s, c*r/(y*s)
    if not (g < 1 < f and a/r*g+b*y/r*f == 1):
        raise ValueError("Independent exact cone-order fixture failed")
    return {"coefficientwise_identities": dict.fromkeys(checks, "0"),
            "g_speed_fixture": str(g), "f_speed_fixture": str(f),
            "positive_null_monotonicity_fixture": str(-a*r*Fraction(7)/(2*Fraction(11)))}
