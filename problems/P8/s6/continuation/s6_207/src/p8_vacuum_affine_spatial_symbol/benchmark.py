"""Exact massive flat collinear benchmark for the finite UV extraction."""

from functools import cache

import sympy as s

from . import extraction

DEGREE = 4


def add(a, b):
    return tuple(s.expand(a[j] + b[j]) for j in range(DEGREE + 1))


def scale(a, c):
    return tuple(s.expand(c * v) for v in a)


def multiply(a, b):
    return tuple(
        s.expand(sum(a[k] * b[j - k] for k in range(j + 1))) for j in range(DEGREE + 1)
    )


def reciprocal(a):
    if a[0] == 0:
        raise ValueError("Nonzero constant term required")
    out = [1 / a[0]]
    for j in range(1, DEGREE + 1):
        out.append(s.expand(-sum(a[k] * out[j - k] for k in range(1, j + 1)) / a[0]))
    return tuple(out)


def square_root_one(a):
    if a[0] != 1:
        raise ValueError("Positive unit square-root branch required")
    out = [s.Integer(1)]
    for j in range(1, DEGREE + 1):
        out.append(s.expand((a[j] - sum(out[k] * out[j - k] for k in range(1, j))) / 2))
    return tuple(out)


@cache
def coefficients():
    m, p = s.symbols("m p", real=True)
    zero = (s.Integer(0),) * 5
    one = (s.Integer(1),) + zero[1:]
    q = (s.Integer(0), s.Integer(0), m * m, s.Integer(0), s.Integer(0))
    omega_k = square_root_one(add(one, q))
    omega_l = square_root_one(
        (s.Integer(1), -2 * p, p * p + m * m, s.Integer(0), s.Integer(0))
    )
    product = multiply(omega_k, omega_l)
    summed = add(omega_k, omega_l)
    inv = reciprocal(summed)
    common = scale(multiply(reciprocal(product), inv), s.Rational(1, 4))
    geom = (s.Integer(1), -p, s.Integer(0), s.Integer(0), s.Integer(0))
    A = add(product, q)
    plus = add(A, geom)
    minus = add(A, scale(geom, -1))
    h = {
        "tensor": multiply(multiply(plus, plus), common),
        "vector": multiply(multiply(q, multiply(summed, summed)), common),
        "scalar": multiply(
            add(
                scale(multiply(minus, minus), s.Rational(1, 3)),
                scale(multiply(A, A), s.Rational(2, 3)),
            ),
            common,
        ),
    }
    rows = {}
    for channel, base in h.items():
        for j in range(5):
            value = base
            for _ in range(j):
                value = multiply(value, inv)
            sign = (1, 0, -1, 0, 1)[j]
            rows[channel, j] = tuple(
                s.expand(sign * value[d]) for d in extraction.retained_degrees(j)
            )
    return m, p, omega_k, omega_l, rows


@cache
def data():
    m, p, wk, wl, rows = coefficients()
    one = (s.Integer(1), s.Integer(0), s.Integer(0), s.Integer(0), s.Integer(0))
    checks = {
        "full_k_frequency_square_root_jet": s.Matrix(multiply(wk, wk))
        - s.Matrix([1, 0, m * m, 0, 0]),
        "full_l_frequency_square_root_jet": s.Matrix(multiply(wl, wl))
        - s.Matrix([1, -2 * p, p * p + m * m, 0, 0]),
        "full_inverse_phase_jet": s.Matrix(
            multiply(add(wk, wl), reciprocal(add(wk, wl)))
        )
        - s.Matrix(one),
        "flat_tensor_zero_transfer_nonzero_subleading": s.Matrix(
            [v.subs(p, 0) for v in rows["tensor", 0]]
        )
        - s.Matrix([s.Rational(1, 2), 0, m * m / 4, 0, -(m**4) / 16]),
        "flat_vector_zero_transfer_mixed_pairs_retained": s.Matrix(
            [v.subs(p, 0) for v in rows["vector", 0]]
        )
        - s.Matrix([0, 0, m * m / 2, 0, -(m**4) / 4]),
        "flat_scalar_zero_transfer_longitudinal_retained": s.Matrix(
            [v.subs(p, 0) for v in rows["scalar", 0]]
        )
        - s.Matrix([s.Rational(1, 12), 0, 5 * m * m / 24, 0, 5 * m**4 / 32]),
        "flat_tensor_fourth_time_jet": rows["tensor", 4][0] - s.Rational(1, 32),
        "flat_scalar_fourth_time_jet": rows["scalar", 4][0] - s.Rational(1, 192),
        "full_three_channel_slot_count": sum(len(v) for v in rows.values()) - 45,
        "all_flat_odd_endpoints_vanish_only_in_this_static_fixture": s.Matrix(
            [v for (channel, j), row in rows.items() if j % 2 for v in row]
        ),
    }
    return {
        "benchmark": "This is an exact massive flat benchmark at a=1, k=r e3 and l=(-r+P)e3, using all three physical channels and the source time jet of order j. It is not the curved CD current or a contact-matched renormalized response.",
        "closed_forms": "Let Omega_k=sqrt(1+m^2x^2), Omega_l=sqrt((1-Px)^2+m^2x^2), S=Omega_k+Omega_l, A=Omega_k Omega_l+m^2x^2 and B=1-Px. The normalized j0 densities are (A+B)^2/(4Omega_k Omega_l S), m^2x^2 S/(4Omega_k Omega_l), and[(A-B)^2/3+2A^2/3]/(4Omega_k Omega_l S) for tensor,vector,scalar. Higher endpoint source-jet coefficients multiply these by cos(j*pi/2)/S^j.",
        "actual_nonzero_coefficients": {
            channel: {str(j): list(rows[channel, j]) for j in range(5)}
            for channel in ("tensor", "vector", "scalar")
        },
        "curved_boundary": "Odd endpoints vanish only for this real static fixture. Curved time-dependent amplitudes and independent detector/source differentiation remain in the general four-sector calculation.",
        "checks": checks,
        "gates": {
            "all_forty_five_benchmark_coefficients_retained": sum(
                len(v) for v in rows.values()
            )
            == 45,
            "spatial_degree_respects_inverse_radius_degree": all(
                s.Poly(v, p).degree() <= d
                for row in rows.values()
                for d, v in enumerate(row)
                if v != 0
            ),
            "massive_vector_mixed_sector_not_deleted": rows["vector", 0][2] != 0,
            "massive_scalar_longitudinal_not_deleted": rows["scalar", 0][0] != 0,
            "not_curved_matching_or_contact_result": True,
        },
    }
