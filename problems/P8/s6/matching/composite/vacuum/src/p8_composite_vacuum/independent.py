"""Independent coefficientwise Fraction determinant and matching checks."""

from fractions import Fraction as Q
from itertools import combinations

from p8_composite_modes.independent import Poly


def checks():
    roots = [Poly.variable(name) for name in ("A", "B", "C", "x")]
    # Coefficients of product_i[(1+d_i)+t(1-d_i)], obtained by direct
    # subset multiplication rather than the main characteristic formula.
    coefficients = []
    for degree in range(5):
        coefficient = Poly()
        for chosen in combinations(range(4), degree):
            term = Poly(1)
            for index, root in enumerate(roots):
                term *= 1-root if index in chosen else 1+root
            coefficient += term
        coefficients.append(coefficient)
    e2 = sum((roots[i]*roots[j] for i, j in combinations(range(4), 2)), Poly())
    e4 = roots[0]*roots[1]*roots[2]*roots[3]
    vacuum_shift = sum((coefficient*Q(-3, 8) for coefficient in coefficients), Poly())
    identities = {
        "literal_beta2_determinant": coefficients[2]-6+2*e2-6*e4,
        "physical_volume_source_identity": sum(coefficients, Poly())-16,
        "all_beta_vacuum_shift_is_relative_independent": vacuum_shift+6,
        "vacuum_relative_FP_and_quartic": coefficients[2]+vacuum_shift+2*e2-6*e4,
    }
    if any(not value.is_zero() for value in identities.values()):
        raise ValueError("An independent Cayley/determinant identity failed")
    shifted = (Q(-3, 8), Q(-3, 8), Q(5, 8), Q(-3, 8), Q(-3, 8))
    first, second = shifted[0]+3*shifted[1]+3*shifted[2]+shifted[3], shifted[4]+3*shifted[3]+3*shifted[2]+shifted[1]
    physical_mass = Q(2)*(shifted[1]+2*shifted[2]+shifted[3])/4
    values = {
        "g_vacuum_constraint": first, "f_vacuum_constraint": second,
        "physical_mass_squared_over_m_squared": physical_mass,
        "physical_Planck_squared_over_M_squared": Q(1, 2),
        "old_CD_centre_curvature_X_slope": (Q(-11, 20)-Q(-9, 20))/Q(1, 5),
        "old_CD_centre_A3": Q(1),
        "CD_endpoint_H_times_tau": 4*Q(1, 2)/(1+Q(1, 2)**2),
        "free_M1_null_budget": Q(8)+Q(1, 10)**2,
        "constant_vacuum_potential_gap": Q(3, 8)-Q(1, 10000),
    }
    if tuple(values.values()) != (0, 0, Q(1, 4), Q(1, 2), Q(-1, 2), 1, Q(8, 5), Q(801, 100), Q(3749, 10000)):
        raise ValueError("An independent physical-clock/matching fixture failed")
    return {"coefficientwise_identities": dict.fromkeys(identities, "0"),
            "Fraction_physical_matching_fixtures": {key: str(value) for key, value in values.items()}}
