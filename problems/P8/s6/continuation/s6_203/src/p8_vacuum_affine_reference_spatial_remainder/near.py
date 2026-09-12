"""Complementary near-momentum estimates without expanding the sharp band."""

from functools import cache

import sympy as s
from p8_vacuum_affine_reference_spatial_analyticity import domain
from p8_vector_state.comparison import AMAX

from . import real


@cache
def constants():
    A = AMAX
    d = domain.DELTA
    m = real.MASS
    E = real.E
    return {
        "internal_near_radius_over_external_P": 2 * A / d,
        "raw_near_P4_upper": 16 * E * A**3 / (27 * d**4),
        "Taylor_orders_zero_to_three_P4_upper": E
        * (4 * A**4 + s.Rational(10, 3) * A**3)
        / (9 * d**4),
        "complete_near_P4_upper": E
        * (4 * A**4 + s.Rational(26, 3) * A**3)
        / (9 * d**4),
        "complete_near_P5_upper": E * A**4 / (9 * d**5 * m),
    }


@cache
def data():
    k, L, m, P, A, d = s.symbols("k L m P A d", positive=True)
    c = constants()
    checks = {
        "zeroth_radial_power": s.integrate(2 * k**3, (k, 0, L)) - L**4 / 2,
        "first_radial_power": s.integrate(k * k, (k, 0, L)) - L**3 / 3,
        "second_radial_power": s.integrate(A * k, (k, 0, L)) - A * L**2 / 2,
        "third_radial_power": s.integrate(A * A, (k, 0, L)) - A * A * L,
        "fourth_radial_logarithm": s.expand_log(
            s.integrate(A**3 / k, (k, m, L)) - A**3 * s.log(L / m)
        ),
        "raw_near_volume_power": s.factor(
            2 * real.E * (2 * P / d) * (2 * A * P / d) ** 3 / 6
            - s.Rational(16, 3) * real.E * A**3 * P**4 / d**4
        ),
        "complete_P4_coefficient": c["raw_near_P4_upper"]
        + c["Taylor_orders_zero_to_three_P4_upper"]
        - c["complete_near_P4_upper"],
        "logarithm_to_fifth_external_power": s.factor(
            real.E * A**3 * P**4 / (2 * d**4) * (2 * A * P / (d * m))
            - real.E * A**4 * P**5 / (d**5 * m)
        ),
    }
    return {
        "near_domain": "NEAR is real |k|>=m and |P|>delta nu(k)/2. Hence nu<2|P|/delta and |k|<L=2Amax|P|/delta. If L<=m, this region is empty.",
        "raw_piece": "Use the all-real full row4e27nu and the near volume bound. Its internal integral is bounded by the displayed raw coefficient times|P|^4, at every external P.",
        "Taylor_coefficients": "For each fixed |k|>=m, the local joint Cauchy theorem gives ||T_q E(P,k)||<=E delta^-q |P|^q nu^(1-q), q0..4, even when the evaluation point P lies outside that Cauchy ball. This is a bound on the polynomial coefficients, not analytic continuation of E to a large complex P.",
        "radial_bounds": "On m<=|k|<=L, use nu<=2|k| for q0 and nu>=|k|/Amax for q2..4. The five integrals after the1/(2pi^2) measure are L^4/(4pi^2),L^3/(6pi^2),Amax L^2/(4pi^2),Amax^2 L/(2pi^2),Amax^3 log(L/m)/(2pi^2).",
        "full_near": "Combine the raw row and all five Taylor orders. Bound log(L/m)<=L/m on L>m and pi^2>9. The complete near integral is belowC4|P|^4+C5|P|^5, with the displayed rational constants. The original both-leg regulator only decreases this absolute majorant.",
        "boundary": "Only the Taylor remainder is bounded; its subtracted polynomial-integrand sector is kept separately, not deleted or matched to local counterterms.",
        "constants": c,
        "checks": checks,
        "gates": {
            "positive_near_radius": c["internal_near_radius_over_external_P"] > 1,
            "C4_display": c["complete_near_P4_upper"] < 2 * 10**52,
            "C5_display": c["complete_near_P5_upper"] < 2 * 10**54,
            "complete_near_below_two_e54": c["complete_near_P4_upper"]
            + c["complete_near_P5_upper"]
            < 2 * 10**54,
        },
    }
