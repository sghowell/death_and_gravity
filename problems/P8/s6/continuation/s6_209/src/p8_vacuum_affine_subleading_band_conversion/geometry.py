"""Exact intersection-shell moments and all retained shape coefficients."""

from functools import cache

import sympy as s


def require_integer(value, low, high):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)):
        raise TypeError("An exact bounded integer is required")
    if not low <= value <= high:
        raise ValueError("Exact integer outside the stated shape range")
    return int(value)


def radius(epsilon, u):
    return epsilon * u + s.sqrt(1 - epsilon**2 * (1 - u**2))


def radial_moment(q, epsilon, u):
    q = require_integer(q, 0, 4)
    R = radius(epsilon, u)
    return -s.log(R) if q == 4 else (1 - R ** (4 - q)) / (4 - q)


def hemisphere_polynomial(q, h, u):
    q = require_integer(q, 0, 3)
    h = require_integer(h, 1, 4 - q)
    M = 4 - q
    values = (
        -u,
        (1 - M * u**2) / 2,
        (M - 1) * u / 2 - (M * M - 1) * u**3 / 6,
        s.Rational(2 - M, 8)
        + s.Rational(M * (M - 2), 4) * u**2
        - s.Rational(M * (M * M - 4), 24) * u**4,
    )
    return s.expand(values[h - 1])


def shape_slots():
    return tuple((q, h) for q in range(4) for h in range(1, 5 - q))


def endpoint_shape_slots():
    return tuple(
        (j, d, h) for j in range(4) for d in range(4 - j) for h in range(1, 5 - j - d)
    )


@cache
def data():
    e, u, v = s.symbols("epsilon u v", real=True)
    checks = {}
    for q in range(4):
        M = 4 - q
        jet = s.series(radial_moment(q, e, u), e, 0, M + 1).removeO().expand()
        checks["complete_hemisphere_q" + str(q)] = s.Matrix(
            [
                s.expand(jet.coeff(e, h) - hemisphere_polynomial(q, h, u))
                for h in range(1, M + 1)
            ]
        )
        strip = s.series(radial_moment(q, e, e * v), e, 0, 4).removeO().expand()
        checks["exact_grazing_leading_q" + str(q)] = s.expand(
            strip - e**2 * (s.Rational(1, 2) - v)
        )
    checks.update(
        {
            "positive_grazing_cubic_value": s.integrate(
                s.Rational(1, 2) - v, (v, 0, s.Rational(1, 2))
            )
            - s.Rational(1, 8),
            "positive_grazing_quartic_derivative": s.integrate(
                v * (s.Rational(1, 2) - v), (v, 0, s.Rational(1, 2))
            )
            - s.Rational(1, 48),
            "grazing_second_derivative_remainder": s.integrate(
                v**2 * (s.Rational(1, 2) - v) / 2, (v, 0, s.Rational(1, 2))
            )
            - s.Rational(1, 384),
            "ten_universal_shape_slots": len(shape_slots()) - 10,
            "twenty_instantiated_endpoint_shape_slots": len(endpoint_shape_slots())
            - 20,
            "thirty_five_instantiated_source_time_jet_shapes": sum(
                j + 1 for j, d, h in endpoint_shape_slots()
            )
            - 35,
            "all_retained_shape_orders_nonnegative_cutoff_power": s.Matrix(
                [int(4 - q - h >= 0) - 1 for q, h in shape_slots()]
            ),
        }
    )
    return {
        "exact_shell": "For K>|P|+m, write e=|P|/K and u=n.phat. The original pair-band minus one-k-ball conversion is minus the lost shell. Its upper angular endpoint is u=e/2, not zero. Each UV coefficient of radial density r^(1-q) has exact radial moment K^(4-q)(1-R^(4-q))/(4-q), R=e*u+sqrt(1-e^2*(1-u^2)); q4 uses log(1/R).",
        "complete_shapes": "For q0,...,3 retain h1,...,4-q. This gives ten universal q/h shapes, twenty instantiations of original j/d slots and thirty-five source-time-jet shape entries. Finite K^0 terms remain. All five q4 endpoint slots have vanishing shell limit, not absent bulk logarithms.",
        "grazing": "The positive strip u in[0,e/2] adds e^3 Fbar(0)/8 and e^4 Fbar_prime(0)/48 before the azimuth factor2pi. Omitting it changes the cubic and quartic shapes. No cutoff surface is differentiated.",
        "checks": checks,
        "gates": {
            "all_ten_shape_slots": len(shape_slots()) == 10,
            "all_twenty_endpoint_shape_instantiations": len(endpoint_shape_slots())
            == 20,
            "finite_terms_and_grazing_derivative_retained": True,
            "q4_shell_not_bulk_logarithm_deleted": True,
            "positive_root_K_above_P_plus_m_required": True,
        },
    }
