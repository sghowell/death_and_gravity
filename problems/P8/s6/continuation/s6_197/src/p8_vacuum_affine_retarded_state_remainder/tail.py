"""Uniform external-momentum regulator tail for the time-ordered remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_spatial_current.bounds import require_band
from p8_vector_state.comparison import AMAX

from . import contact, memory

DISPLAY = s.Integer(10) ** 27


@cache
def constants():
    m = s.Integer(1000)
    c = memory.constants()
    low = 72 * (memory.CREF * memory.CERR * AMAX**4 / 9 + memory.CERR**2 * 100 / m**6)
    high = s.Rational(36, 100) * (
        memory.CREF * memory.CERR * c["J4_pi_lower_bound"]
        + memory.CERR**2 * c["J10_pi_lower_bound"]
    )
    contact_tail = contact.constants()["contact_tail_K_squared_numerator"]
    return {
        "low_external_memory_tail_numerator": low,
        "high_external_memory_tail_numerator": high,
        "all_external_memory_tail_numerator": low + high,
        "contact_K_squared_tail_numerator": contact_tail,
        "complete_K_tail_numerator": low + high + contact_tail / m,
    }


def error_bound(K):
    return DISPLAY / require_band(K)


@cache
def data():
    K, r = s.symbols("K r", positive=True)
    c = constants()
    checks = {
        "fourth_radial_half_band_tail": s.integrate(r**-2, (r, K / 2, s.oo)) - 2 / K,
        "tenth_radial_half_band_tail": s.integrate(r**-8, (r, K / 2, s.oo))
        - 128 / (7 * K**7),
        "fifth_radial_full_band_tail": s.integrate(r**-3, (r, K, s.oo))
        - 1 / (2 * K * K),
        "same_both_leg_memory_and_single_leg_contact_regulator": error_bound(10**16)
        - 10**11,
    }
    return {
        "low_external": "For |P|<=K/2 the removed two-particle union forces BOTH |k|,|l|>K/2. Keep the complete memory product; its internal tail is36L(Cref Cerr J4_tail(K/2)+Cerr^2 J10_tail(K/2)).",
        "low_external_radials": "J4_tail(K/2)<=Amax^4/(pi^2 K)<Amax^4/(9K), J10_tail(K/2)<100/K^7. After L<2(1+|P|^2), K>=1000 gives the recorded low-external1/K numerator.",
        "high_external": "For |P|>K/2 use the full internal memory integral and L/(1+|P|^2)<1/(100K). This retains all external momenta, including a small internal leg paired with a large one.",
        "contact": "The local contact has one internal momentum. Its removed tail is at most MODE Amax^5/(4pi^2 K^2)<1e14/K^2. It is not the two-leg overlap and is not dropped at large external transfer.",
        "result": "For K>=1000 the full memory-plus-contact actual-minus-reference current has regulator error<1e27 M[D]M[Gamma]/K. No time-order factor has been replaced or differentiated.",
        "constants": c,
        "checks": checks,
        "gates": {
            "low_external_memory_display": c["low_external_memory_tail_numerator"]
            < 10**26,
            "high_external_memory_display": c["high_external_memory_tail_numerator"]
            < 10**21,
            "all_external_memory_display": c["all_external_memory_tail_numerator"]
            < 10**26,
            "full_contact_tail_display": c["contact_K_squared_tail_numerator"] < 10**14,
            "complete_memory_and_contact_display": c["complete_K_tail_numerator"]
            < DISPLAY,
        },
    }
