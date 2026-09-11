"""Full local second-metric contact comparison, not a discarded seagull."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vector_state.comparison import AMAX

FIELD = reference.FIELD_NORM
REMAINDER = s.Integer(2) * 10**6
MODE = s.Rational(15, 2) * FIELD**2 * REMAINDER
DISPLAY = s.Integer(10) ** 8


@cache
def constants():
    A, m = AMAX, reference.MASS
    return {
        "complete_three_mode_contact_numerator": MODE,
        "J5_pi_lower_bound": A**3 / (54 * m * m),
        "contact_L2_coefficient": MODE * A**3 / (54 * m * m),
        "contact_tail_K_squared_numerator": MODE * A**5 / 36,
    }


@cache
def data():
    v = s.Matrix(s.symbols("v0:4", complex=True))
    e = s.Matrix(s.symbols("e0:4", complex=True))
    d = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]])
    g = s.Matrix([[2, -1, 1], [-1, -3, 2], [1, 2, 1]]) / 7
    H = (d * g + g * d) / 2
    eps, eta = s.symbols("epsilon eta", real=True)
    X = eps * d + eta * g
    exp2 = s.eye(3) + X + X * X / 2
    c = constants()
    y = s.Symbol("y", positive=True)
    checks = {
        "full_actual_reference_covariance_difference": (v + e) * (v + e).conjugate().T
        - v * v.conjugate().T
        - v * e.conjugate().T
        - e * v.conjugate().T
        - e * e.conjugate().T,
        "noncommuting_full_metric_contact": exp2.applyfunc(
            lambda x: s.expand(x).coeff(eps, 1).coeff(eta, 1)
        )
        - H,
        "three_modes_half_Hamiltonian_and_three_terms": s.Rational(3, 2)
        * FIELD**2
        * (4 * REMAINDER + REMAINDER)
        - MODE,
        "fifth_complete_massive_radial_integral": s.integrate(
            y * y / (1 + y * y) ** s.Rational(5, 2), (y, 0, s.oo)
        )
        - s.Rational(1, 3),
        "contact_does_not_drop_longitudinal_mode": MODE
        - s.Rational(15, 2) * 4000**2 * (2 * 10**6),
    }
    return {
        "full_local_contact": "Retain -<H_DGamma> from the actual spatial Hamiltonian. The comparison replaces each physical one-mode readout by its same constant-alpha W8 reference, with H=(D Gamma+Gamma D)/2 and the positive contact mass sign.",
        "readout": "After the exact a^(3/2) physical readout rescaling, ||v_ref||<=2*4000 sqrt(nu), ||e||<=4000 R nu^(-11/2), R<2e6. The contact feature norm is <=||D||F||Gamma||F. The temporal constraint mode remains in v, with zero direct unimodular vertex.",
        "three_terms": "Half the sum over all3 polarizations bounds the complete covariance difference by(3/2)4000^2[4R nu^-5+R^2 nu^-11]. Since R/nu^6<1, this is<2.4e14 nu^-5.",
        "result": "The full contact difference is absolutely integrable, using J5=Amax^3/(6pi^2 m^2). Its spacetime pairing is<1e8 ||D||L2||Gamma||L2; its discarded one-momentum tail is<1e14 ||D||L2||Gamma||L2/K^2.",
        "reference_boundary": "The reference contact and reference memory are defined directly from the same finite-regulator mode formula. They are not asserted to be derivatives of a fixed renormalized effective action or the covariance of a new CCR state.",
        "constants": c,
        "checks": checks,
        "gates": {
            "actual_all_order_remainder_strict": 0
            < reference.data()["actual_state_remainder_coefficient"]
            < REMAINDER,
            "relative_error_less_than_one": REMAINDER / reference.MASS**6 < 1,
            "complete_contact_strict_display": c["contact_L2_coefficient"] < DISPLAY,
            "contact_tail_strict_display": c["contact_tail_K_squared_numerator"]
            < 10**14,
            "noncommuting_contact_fixture": d * g != g * d,
        },
    }
