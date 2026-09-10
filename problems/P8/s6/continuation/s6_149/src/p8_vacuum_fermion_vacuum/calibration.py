"""Actual finite vacuum primitive enclosures with both flavor multiplicities."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_ms_mass import calibration as previous
from p8_vacuum_fermion_self_energy_chord.tail import rational

from . import scalar


def gauge_enclosure(m, a, Q):
    m, a, Q = map(rational, (m, a, Q))
    if m < 2 or a < 0 or not 0 < Q <= 144:
        raise ValueError("Need m>=2,a>=0 and 0<Q_lower<=144")
    return 42 * a * s.Rational(4, 3) * 18 * m**4 / Q**2


@cache
def data():
    p = previous.data()["actual_reference_parameters"]
    m, Y, a = p["mF"], p["Y_upper"], p["a_upper"]
    sc = scalar.enclosure(m, Y, 144)
    g = gauge_enclosure(m, a, 144)
    total = sc["complete_scalar_vacuum_absolute_upper"] + g
    leading_lower = 63 * m**4 / 256
    relative = total / leading_lower
    return {
        "actual_reference_parameters": p,
        "scalar_vacuum_enclosure": sc,
        "all_flavor_gauge_vacuum_absolute_upper": g,
        "both_vacuum_primitive_absolute_upper": total,
        "all_fourteen_flavor_one_loop_vacuum_lower": leading_lower,
        "relative_to_one_loop_vacuum_upper": relative,
        "checks": {
            "same_actual_mass": m - 10**200,
            "scalar_and_gauge_added_once": total
            - sc["complete_scalar_vacuum_absolute_upper"]
            - g,
            "gauge_all_flavor_not_active_only": g
            - 7 * (6 * a * s.Rational(4, 3) * 18 * m**4 / 144**2),
            "one_loop_all_flavor_coefficient": s.Rational(3, 2) * 42 - 63,
            "physical_Q_upper_from_pi_upper": 16 * 4**2 - 256,
        },
        "bounds": {
            "actual_mass_in_domain": bool(m >= 2),
            "actual_Y_upper": bool(0 < p["actual_Y"] < Y),
            "actual_a_upper": bool(0 < p["actual_a"] < a),
            "all_flavor_vacuum_upper_below_one_e_595": bool(0 < total < 10**595),
            "vacuum_relative_to_one_loop_below_one_e_minus_203": bool(
                0 < relative < s.Rational(1, 10**203)
            ),
            "actual_scalar_mass_correction_below_one_e_194": bool(
                0 < sc["actual_scalar_vacuum_mass_correction_absolute_upper"] < 10**194
            ),
            "finite_scalar_vacuum_remainder_below_one_e_minus_607": bool(
                0
                < sc["finite_scalar_vacuum_remainder_absolute_upper"]
                < s.Rational(1, 10**607)
            ),
            "actual_gauge_flavor_multiplicity_matters": bool(
                g > sc["complete_scalar_vacuum_absolute_upper"]
            ),
        },
        "scope": "Both fermion vacuum primitive contributions, not the complete vacuum-reference counterterm or a global potential/semiclassical solution. All remaining canonical/matching/truncation/V/G/B obligations stay open.",
    }
