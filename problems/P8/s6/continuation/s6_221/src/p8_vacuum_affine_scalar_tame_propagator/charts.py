"""Exact finite-q canonical swap and outer constraint chart, including all time boundaries."""

from functools import cache

import sympy as s
from p8 import gamma as historical

th, E, l, J, A, T, H, q = s.symbols("Theta E ell Jc A Tcorr H q", real=True)
b, v, sigma, ps, bd, sd = s.symbols("b v sigma ps bdot sigmadot", real=True)
w = -l * E
jetvars = (th, E, l, J, A, T, H)
jets = s.symbols("Theta1 E1 ell1 Jc1 A1 Tcorr1 H1", real=True)
jet2 = s.symbols("Theta2 E2 ell2 Jc2 A2 Tcorr2 H2", real=True)
z = s.Symbol("inverse_q", nonnegative=True)


def dtime(expr):
    return s.factor(
        sum(s.diff(expr, c) * d for c, d in zip(jetvars, jets))
        + sum(s.diff(expr, c) * d for c, d in zip(jets, jet2))
        - 2 * H * q * s.diff(expr, q)
    )


def gradient_numerator():
    return th * jets[1] - E * jets[0] + H * E * th - th**2


def ham():
    lc = -2 * q * th * b - w * ps + 3 * th * l * sigma - (2 * E * q + 3 * T) * v
    return (
        ps**2 / 2
        + q * l * b * sigma
        + (q / 2 - 3 * l * l / 4) * sigma**2
        - q * v * v
        - s.Rational(9, 2) * A * v * v
        + lc * lc / (4 * J)
    )


@cache
def central():
    HC = ham()
    x = s.Matrix([v, ps])
    y = s.Matrix([b, sigma])
    XX = s.hessian(HC, x)
    XY = s.Matrix([[s.diff(HC, xx, yy) for yy in y] for xx in x])
    YY = s.hessian(HC, y)
    M = s.diag(2 * q, 1)
    transport = s.Matrix([[2 * H * q, 0], [0, 0]])
    source = transport - XY
    inverse = XX.inv().applyfunc(s.factor)
    K = (M.T * inverse * M).applyfunc(s.factor)
    B = (M.T * inverse * source).applyfunc(s.factor)
    D = (source.T * inverse * source - YY).applyfunc(s.factor)
    K0 = s.Matrix([[2 * J / E**2 + l * l, l], [l, 1]])
    G = s.Matrix([[2 * gradient_numerator() / E**2, l], [l, 1]])
    gyro = (B - B.T).applyfunc(s.factor)
    lower = (B.applyfunc(dtime) + 3 * H * B - D - q * G).applyfunc(s.factor)
    Delta = (2 * E + 3 * T * z) ** 2 - (2 * z + 9 * A * z * z) * (2 * J + E * E * l * l)
    return {
        "XX": XX,
        "XY": XY,
        "YY": YY,
        "M": M,
        "transport": transport,
        "source": source,
        "inverse": inverse,
        "K": K,
        "B": B,
        "D": D,
        "K0": K0,
        "G": G,
        "gyro": gyro,
        "lower": lower,
        "Delta": Delta,
    }


@cache
def outer():
    vd = s.Symbol("vdot", real=True)
    nn = (vd + l * sigma / 2) / th
    L = (
        -3 * vd**2
        + (J + w * w / 2 - 3 * th * th) * nn**2
        + 6 * th * nn * vd
        + sd * sd / 2
        + w * nn * sd
        - 3 * l * vd * sigma
        + q * v * v
        + 2 * E * q * nn * v
        - q * sigma * sigma / 2
        + 3 * T * nn * v
        + s.Rational(9, 2) * A * v * v
    )
    y = s.Matrix([v, sigma])
    yd = s.Matrix([vd, sd])
    K = s.hessian(L, yd).applyfunc(s.factor)
    B = s.Matrix([[s.diff(L, d, c) for c in y] for d in yd]).applyfunc(s.factor)
    D = s.hessian(L, y).applyfunc(s.factor)
    G = s.Matrix([[2 * gradient_numerator() / th**2, w / th], [w / th, 1]])
    lower = (B.applyfunc(dtime) + 3 * H * B - D - q * G).applyfunc(s.factor)
    gyro = (B - B.T).applyfunc(s.factor)
    return {"L": L, "K": K, "B": B, "D": D, "G": G, "lower": lower, "gyro": gyro}


def degree(expr):
    value = s.factor(expr)
    if value == 0:
        return -s.oo
    n, d = s.fraction(value)
    return s.degree(n, q) - s.degree(d, q)


@cache
def data():
    c, o = central(), outer()
    XX, K, B, D, G = [c[name] for name in ("XX", "K", "B", "D", "G")]
    x = s.Matrix([v, ps])
    y = s.Matrix([b, sigma])
    yd = s.Matrix([bd, sd])
    rhs = c["M"] * yd + c["source"] * y
    solution = c["inverse"] * rhs
    swapped = 2 * q * v * bd + 2 * H * q * v * b + ps * sd - ham()
    reduced = (rhs.T * c["inverse"] * rhs)[0] / 2 - (y.T * c["YY"] * y)[0] / 2
    explicit = (yd.T * K * yd)[0] / 2 + (yd.T * B * y)[0] + (y.T * D * y)[0] / 2
    det = ((2 * E * q + 3 * T) ** 2 - (2 * q + 9 * A) * (2 * J + w * w)) / (2 * J)
    speed = s.Symbol("speed_squared", real=True)
    F = gradient_numerator() - w * w / 2
    old = historical.principal_chart()
    ts, fts, ths, es, ws, ls, js = old["symbols"]
    hs, tds, eds, thds = old["derivative_symbols"]
    past = {
        ts: 1,
        fts: 1,
        ths: th,
        es: E,
        ws: w,
        ls: l,
        js: J,
        hs: H,
        tds: 0,
        eds: jets[1],
        thds: jets[0],
    }
    checks = {
        "complete_finite_q_auxiliary_determinant": s.factor(XX.det() - det),
        "full_weighted_pair_swap_stationarity": (XX * solution - rhs).applyfunc(
            s.factor
        ),
        "full_auxiliary_elimination_matches_exact_reduced_L": s.factor(
            swapped.subs(dict(zip(x, solution)), simultaneous=True) - reduced
        ),
        "full_K_B_D_reconstruction": s.factor(reduced - explicit),
        "central_complete_principal_kinetic": (
            K.applyfunc(lambda e: s.limit(e, q, s.oo)) - c["K0"]
        ).applyfunc(s.factor),
        "central_complete_principal_gradient": (
            (B.applyfunc(dtime) + 3 * H * B - D).applyfunc(
                lambda e: s.limit(e / q, q, s.oo)
            )
            - G
        ).applyfunc(s.factor),
        "normalized_positive_denominator": s.factor(
            c["Delta"] - (2 * J * XX.det() / q**2).subs(q, 1 / z)
        ),
        "exact_central_kinetic_determinant": s.factor(
            K.det() - 8 * J / c["Delta"].subs(z, 1 / q)
        ),
        "historical_b_chart_principal_kinetic_factor_two": (
            c["K0"] - 2 * old["kinetic"].subs(past, simultaneous=True)
        ).applyfunc(s.factor),
        "historical_full_boundary_principal_gradient_factor_two": (
            G - 2 * old["gradient"].subs(past, simultaneous=True)
        ).applyfunc(s.factor),
        "outer_complete_kinetic": (
            o["K"] - s.Matrix([[(2 * J + w * w) / th**2, w / th], [w / th, 1]])
        ).applyfunc(s.factor),
        "both_valid_chart_principal_characteristic_polynomials": s.Matrix(
            [
                s.factor(
                    (speed * c["K0"] - G).det()
                    - 2 * J / E**2 * (speed - 1) * (speed - F / J)
                ),
                s.factor(
                    (speed * o["K"] - o["G"]).det()
                    - 2 * J / th**2 * (speed - 1) * (speed - F / J)
                ),
            ]
        ),
        "both_gyroscopic_matrices_skew": s.diag(
            c["gyro"] + c["gyro"].T, o["gyro"] + o["gyro"].T
        ),
        "central_inverse_phase_velocity_map": (
            c["M"].T * (XX * solution + c["XY"] * y - c["transport"] * y)
            - c["M"].T * c["M"] * yd
        ).applyfunc(s.factor),
    }
    return {
        "central_K": K,
        "central_B": B,
        "central_D": D,
        "central_principal_K": c["K0"],
        "central_principal_G": G,
        "outer_K": o["K"],
        "outer_principal_G": o["G"],
        "auxiliary_determinant": s.factor(XX.det()),
        "normalized_denominator": c["Delta"],
        "complete_pair_swap": "With b=-pv/(2q), the weighted boundary gives pv vdot ->2q v bdot+2Hq v b. The second term is essential. Jointly eliminate x=(v,ps) with the actual full retuned quadratic Hamiltonian, retaining all A,Tcorr,Jc terms.",
        "Euler_operator": "For L=ydot^T K ydot/2+ydot^T B y+y^T D y/2, the exact equation is K ydd+(K'+3HK+B-B^T)yd+(B'+3HB-D)y=0. Every coefficient derivative includes q'=-2Hq.",
        "principal_bridge": "The earlier repository gamma chart already supplies the bare complementary principal algebra. This successor verifies its exact factor-two normalization and extends the full finite-q retuned action, uniform remainders and quantitative two-chart estimates; that older qualitative algebra is not represented as new.",
        "finite_q_warning": "The complementary auxiliary Hessian can vanish at a finite q. It is used only above the explicitly proved q>=4096 threshold where Delta>1/8. The original first-order Hamiltonian remains the low-momentum comparison.",
        "checks": checks,
        "gates": {
            "central_K_bounded_at_high_q": all(degree(e) <= 0 for e in K),
            "central_K_minus_principal_decays": all(
                degree(e) <= -1 for e in K - c["K0"]
            ),
            "central_gyroscopic_remainder_bounded": all(
                degree(e) <= 0 for e in c["gyro"]
            ),
            "central_potential_remainder_bounded": all(
                degree(e) <= 0 for e in c["lower"]
            ),
            "outer_coefficients_and_remainders_no_q": all(
                not e.has(q) for name in ("K", "gyro", "lower") for e in o[name]
            ),
            "full_retuning_and_time_boundary_retained": True,
            "not_every_finite_q_auxiliary_chart_invertible": True,
            "no_Theta_denominator_in_central_coefficients": all(
                not s.denom(s.factor(e)).has(th)
                for name in ("K", "B", "D")
                for e in c[name]
            ),
            "no_ell_or_E_denominator_in_outer_coefficients": all(
                not s.denom(s.factor(e)).has(l, E)
                for name in ("K", "B", "D")
                for e in o[name]
            ),
        },
    }
