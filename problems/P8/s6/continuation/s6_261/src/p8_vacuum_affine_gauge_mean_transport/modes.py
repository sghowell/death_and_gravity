"""Exact finite Fourier polynomials with all product modes retained."""

from functools import cache

import sympy as s

ZERO = (0, 0, 0)
ALPHA = s.Rational(9, 4)


def tidy(field):
    return {k: v for k, c in sorted(field.items()) if (v := s.expand(c)) != 0}


def add(*fields):
    out = {}
    for field in fields:
        for k, c in field.items():
            out[k] = out.get(k, 0) + c
    return tidy(out)


def scale(c, field):
    return tidy({k: c * v for k, v in field.items()})


def mul(first, second):
    out = {}
    for k, a in first.items():
        for p, b in second.items():
            q = tuple(k[j] + p[j] for j in range(3))
            out[q] = out.get(q, 0) + a * b
    return tidy(out)


def derivative(field, j):
    return tidy({k: s.I * k[j] * a for k, a in field.items()})


def wave(k, amplitude=1, sine=False):
    opposite = tuple(-a for a in k)
    if k == ZERO:
        return {} if sine else {ZERO: s.sympify(amplitude)}
    if sine:
        return {k: amplitude / (2 * s.I), opposite: -amplitude / (2 * s.I)}
    return {k: s.sympify(amplitude) / 2, opposite: s.sympify(amplitude) / 2}


def divergence(vector):
    return add(*(derivative(vector[j], j) for j in range(3)))


def dot(first, second):
    return add(*(mul(a, b) for a, b in zip(first, second, strict=True)))


def gradient(field):
    return [derivative(field, j) for j in range(3)]


def mean(field):
    return field.get(ZERO, s.Integer(0))


def inverse(vector):
    modes = sorted(set().union(*(set(f) for f in vector)))
    out = [{}, {}, {}]
    for mode in modes:
        rhs = s.Matrix([f.get(mode, 0) for f in vector])
        if mode == ZERO:
            if rhs != s.zeros(3, 1):
                raise ValueError(
                    "The nonzero-mode ghost inverse cannot solve a constant source"
                )
            continue
        k = s.Matrix(mode)
        k2 = k.dot(k)
        value = (s.eye(3) - k * k.T / (4 * k2)) * rhs / k2
        for i in range(3):
            out[i][mode] = value[i]
    return [tidy(f) for f in out]


def lie_contravariant_density(Q, xi):
    div = divergence(xi)
    return [
        [
            add(
                dot(xi, gradient(Q[i][j])),
                scale(-1, add(*(mul(Q[i][k], derivative(xi[j], k)) for k in range(3)))),
                scale(-1, add(*(mul(Q[j][k], derivative(xi[i], k)) for k in range(3)))),
                scale(s.Rational(2, 3), mul(Q[i][j], div)),
            )
            for j in range(3)
        ]
        for i in range(3)
    ]


def ghost_from_lie(Q, xi):
    lie = lie_contravariant_density(Q, xi)
    return [add(*(derivative(lie[i][j], j) for j in range(3))) for i in range(3)]


def matrix_action(Q, vector):
    return [add(*(mul(Q[i][j], vector[j]) for j in range(3))) for i in range(3)]


def vector_residual(vector):
    modes = sorted(set().union(*(set(f) for f in vector)) | {ZERO})
    return s.Matrix([f.get(k, 0) for k in modes for f in vector])


def evaluate(field, point):
    return sum(
        complex(c)
        * __import__("cmath").exp(1j * sum(k[j] * point[j] for j in range(3)))
        for k, c in field.items()
    )


@cache
def packet():
    v = add(
        wave((1, 0, 0)),
        wave((0, 1, 1), s.Rational(1, 3), True),
        wave((1, 0, 2), s.Rational(1, 5)),
    )
    n = add(
        wave((1, 0, 0), s.Rational(1, 7)),
        wave((0, 1, 1), s.Rational(1, 4), True),
        wave((2, 0, 0), s.Rational(1, 11)),
    )
    h = [[{} for j in range(3)] for i in range(3)]
    h[1][2] = h[2][1] = wave((1, 0, 0), s.Rational(1, 7))
    h[0][2] = h[2][0] = wave((0, 1, 0), s.Rational(1, 9))
    Q1 = [[scale(-1, h[i][j]) for j in range(3)] for i in range(3)]
    unit = [
        [{ZERO: s.Integer(i == j)} if i == j else {} for j in range(3)]
        for i in range(3)
    ]
    xi1 = inverse([scale(-3, f) for f in gradient(v)])
    M1 = ghost_from_lie(Q1, xi1)
    source2 = [
        add(scale(3, z), scale(-1, m))
        for z, m in zip(matrix_action(h, gradient(v)), M1, strict=True)
    ]
    xi2 = inverse(source2)
    dv1 = scale(s.Rational(1, 3), divergence(xi1))
    dv2 = add(dot(xi1, gradient(v)), scale(s.Rational(1, 3), divergence(xi2)))
    dn2 = dot(xi1, gradient(n))
    direct2 = [
        add(m0, m1, scale(3, force))
        for m0, m1, force in zip(
            ghost_from_lie(unit, xi2), M1, matrix_action(Q1, gradient(v)), strict=True
        )
    ]
    shortened = [
        add(
            add(
                *(
                    mul(h[j][k], derivative(derivative(xi1[i], j), k))
                    for j in range(3)
                    for k in range(3)
                )
            ),
            scale(
                s.Rational(1, 3),
                add(*(mul(h[i][j], derivative(divergence(xi1), j)) for j in range(3))),
            ),
        )
        for i in range(3)
    ]
    Cvv, Cvn = mean(mul(v, v)), mean(mul(v, n))
    volume1 = 3 * mean(dv2) - 6 * mean(dn2)
    volume2 = 9 * mean(mul(v, dv1)) - 18 * mean(mul(n, dv1))
    # Off-slice lower terms are differentiated in ghost_from_lie. Only
    # div(Q1)=0 permits their absence from the reduced M1 expression.
    return {
        "whole_three_direction_scalar_fixture": v,
        "whole_nonindependent_lapse_fixture": n,
        "whole_two_noncollinear_TT_metric_modes": h,
        "whole_first_gauge_displacement": xi1,
        "whole_second_gauge_displacement": xi2,
        "whole_first_log_volume_variation": dv1,
        "whole_second_log_volume_variation": dv2,
        "whole_second_lapse_variation": dn2,
        "whole_finite_Cvv": Cvv,
        "whole_finite_Cvn": Cvn,
        "whole_linear_physical_volume_mean_contact": volume1,
        "whole_quadratic_physical_volume_mean_contact": volume2,
        "checks": {
            "whole_TT_trace": mean(add(*(h[i][i] for i in range(3)))),
            "whole_TT_divergence": vector_residual(
                [add(*(derivative(h[i][j], j) for j in range(3))) for i in range(3)]
            ),
            "whole_first_gauge_equation_from_density_Lie_derivative": vector_residual(
                [
                    add(a, scale(3, b))
                    for a, b in zip(ghost_from_lie(unit, xi1), gradient(v), strict=True)
                ]
            ),
            "whole_second_gauge_equation_all_product_modes": vector_residual(direct2),
            "whole_second_operator_including_off_slice_derivation": vector_residual(
                [add(a, scale(-1, b)) for a, b in zip(M1, shortened, strict=True)]
            ),
            "whole_first_trace_transport": vector_residual(
                [add(dv1, scale(-s.Rational(3, 4), v))]
            ),
            "whole_second_mean_trace_transport": s.factor(mean(dv2) + ALPHA * Cvv),
            "whole_second_mean_lapse_transport": s.factor(mean(dn2) + ALPHA * Cvn),
            "whole_second_displacement_divergence_average": mean(divergence(xi2)),
            "whole_physical_volume_contact_cancellation": s.factor(volume1 + volume2),
            "whole_second_gauge_rhs_has_no_constant_mode": s.Matrix(
                [mean(f) for f in source2]
            ),
        },
        "gates": {
            "nonlinear_gauge_displacement_is_not_zero": any(xi2),
            "mode_products_not_projected_into_input_set": bool(
                set().union(*(set(f) for f in xi2))
                - set(v)
                - set().union(*(set(f) for row in h for f in row))
            ),
            "finite_Cvv_strictly_positive": Cvv > 0,
            "finite_lapse_cross_covariance_not_deleted": Cvn != 0,
            "both_physical_volume_contacts_are_nonzero": volume1 != 0 and volume2 != 0,
            "finite_Fourier_polynomials_not_a_closed_gauge_regulator": True,
        },
    }
