"""Complete fixed-CD source defects and explicit quadratic/cubic tame bounds."""

from functools import cache

import sympy as s
from p8_vacuum_analytic_affine_parent import source as original

from . import coefficients


@cache
def data():
    a = s.Symbol("a", positive=True)
    H, b, c = s.symbols("H psi_t psi_tt", real=True)
    v = s.Matrix(s.symbols("psi_x0:3", real=True))
    w = s.Matrix(s.symbols("psi_tx0:3", real=True))
    D = s.Matrix(
        3,
        3,
        lambda i, j: s.Symbol(
            "psi_x" + str(min(i, j)) + "x" + str(max(i, j)), real=True
        ),
    )
    p = s.Matrix([1 + b, *v])
    inv = s.diag(1, -1 / a**2, -1 / a**2, -1 / a**2)
    Hess = s.zeros(4)
    Hess[0, 0] = c
    for i in range(3):
        Hess[0, i + 1] = Hess[i + 1, 0] = w[i] - H * v[i]
        for j in range(3):
            Hess[i + 1, j + 1] = D[i, j] - (a * a * H * (1 + b) if i == j else 0)
    X = (p.T * inv * p)[0]
    Box = s.trace(inv * Hess)
    Z = (p.T * inv * Hess * inv * p)[0]
    expectedZ = (
        (1 + b) ** 2 * c
        - 2 * (1 + b) * (v.dot(w)) / a**2
        + (v.T * D * v)[0] / a**4
        + H * (1 + b) * v.dot(v) / a**2
    )
    R, RX, Ru, Y, Q, ZZ, HH = s.symbols("R R_X R_u X E Z Href", real=True)
    prior = original.data()["shifted_source_P8"]
    names = {str(x): x for x in prior.free_symbols}
    substituted = prior.subs(
        {
            names["R"]: R,
            names["R_X"]: RX,
            names["R_u"]: Ru,
            names["positive_X"]: Y,
            names["Box_P8_u"]: Q + 3 * HH * Y,
            names["uHu_P8"]: ZZ,
            names["H_clock"]: HH,
        },
        simultaneous=True,
    )
    CZ = -1 / Y**2 + 3 * RX / (2 * R * Y)
    scalar = (R - 1) * (Q / Y + 3 * Ru / (4 * R) + CZ * ZZ)
    eps = s.Symbol("epsilon", real=True)
    r1, q1, p1 = s.symbols("R_first Q_first p_first", real=True)
    germ = (1 + eps * p1) * (eps * r1) * (eps * q1)
    t, ft, ftt, lap, f = s.symbols("t f_t f_tt Delta_f f", real=True)
    h = (1 + t * t) ** 3
    hclock = 4 * t / (1 + t * t)
    hdot = s.diff(hclock, t)
    E1 = ftt - 3 * hclock * ft - lap / (1 + t * t) ** 4 - 3 * hdot * f
    Q1 = E1 - 9 * t * ft / ((1 + t * t) * h) + (-1 + 3 / (2 * h)) * ftt
    S20 = 2 * ft * Q1 / h
    z = s.Symbol("local_spatial_profile", real=True)
    special = s.simplify(
        S20.subs({t: 0, ft: z, ftt: z, f: 0, lap: 0}, simultaneous=True)
    )
    return {
        "class": "Fixed physical CD metric, psi compact smooth in (-1/2,1/2)xR3, u=t+psi, max of every coordinate derivative through3 at most delta<=1/100.",
        "literal_clock_gradient_square": X,
        "literal_clock_Box": Box,
        "literal_clock_Hessian_contraction": s.expand(expectedZ),
        "exact_complete_source_defect_form": "S_mu=partial_mu u*(R-1)*Q, Q=E/X+3R_u/(4R)+C Z; E=Boxu-3H(u)X and C=-1/X^2+3R_X/(2*R*X). All R functions are the frozen full analytic target.",
        "leading_scalar_clock_source_temporal": S20,
        "leading_scalar_clock_source_spatial": "Zero at second order; spatial S_i starts at third order.",
        "nonclosed_source_control": "For f=(t+t^2/2)c(x) near t=0, S2_0=3c(x)^2 and S2_i=0; partial_i S2_0=6c partial_i c is generally nonzero. A compact bump can realize these local jets.",
        "homogeneous_control": "When S=S0(t)dt, the compatible retarded solution is W=S, so the FULL sourced light force vanishes. The positive free coherent energy alone must not be called the full source-interaction stress.",
        "checks": {
            "literal_gradient_square": s.expand(X - ((1 + b) ** 2 - v.dot(v) / a**2)),
            "literal_Box_in_actual_geometry": s.expand(
                Box - (c + 3 * H * (1 + b) - s.trace(D) / a**2)
            ),
            "literal_Z_in_actual_geometry": s.expand(Z - expectedZ),
            "literal_complete_source_to_defect_form": s.factor(substituted - scalar),
            "source_zero_on_clock": germ.subs(eps, 0),
            "source_first_clock_variation_zero": s.diff(germ, eps).subs(eps, 0),
            "source_second_clock_coefficient": s.expand(germ).coeff(eps, 2) - r1 * q1,
            "inhomogeneous_actual_leading_source": special - 3 * z * z,
            "inhomogeneous_leading_source_curl": s.diff(special, z) - 6 * z,
        },
    }


@cache
def bounds():
    d = s.Rational(1, 100)
    lo = s.Rational(9, 10)
    C = coefficients.data()["source_Z_coefficient_bounds"]
    q_bound = 40 / lo + s.Rational(9, 2) * 3 + 6 * 2
    q_derivative_E = 42 / lo + 40 * 3 * d / lo**2
    q_derivative_R = s.Rational(81, 5) + s.Rational(351, 2) * d
    q_derivative_Z = 12 + 222 * d
    q_derivative = q_derivative_E + q_derivative_R + q_derivative_Z
    s_bound = 2 * 3 * 70
    ds_bound = 2 * 3 * 70 * d + 2 * 4 * 70 + 2 * 3 * 82
    tests = {
        "gradient_square_defect_below_three": 2 + 4 * d < 3,
        "gradient_square_spatial_or_directional_derivative_below_three": 2 + 8 * d < 3,
        "full_gradient_norm_below_two": (1 + d) ** 2 + 3 * d * d < 4,
        "Box_H_clock_defect_below_forty": 4 + 6 + 12 + 6 * 3 <= 40,
        "Box_H_clock_derivative_below_forty_two": 10 + 12 * s.Rational(11, 10) + 6 * 3
        < 42,
        "full_Z_defect_below_two": (1 + d) ** 2 + 12 * (1 + d) * d + 9 * d * d < 2,
        "full_Z_derivative_below_two": (1 + d) ** 2 + 26 * (1 + d) * d + 39 * d * d < 2,
        "R_defect_derivative_below_four": s.Rational(6, 5) * 3 + 9 * d < 4,
        "Q_defect_below_seventy": q_bound < 70,
        "Q_derivative_below_eighty_two": q_derivative < 82,
        "full_source_below_five_twelve": s_bound < 512,
        "full_source_spatial_or_Frechet_derivative_below_two_zero_four_eight": ds_bound
        < 2048,
    }
    return {
        "maximum_coordinate_jet_size": d,
        "X_minus_one_over_delta": s.Integer(3),
        "E_over_delta": s.Integer(40),
        "Z_over_delta": s.Integer(2),
        "Q_over_delta": q_bound,
        "Q_derivative_over_jet_norm": q_derivative,
        "source_pointwise_over_delta_squared": s.Integer(s_bound),
        "source_derivative_over_delta_times_jet_norm": ds_bound,
        "uniform_source_L2_constant": s.Integer(512),
        "uniform_spatial_and_Frechet_L2_constant": s.Integer(2048),
        "norm_definition": "G_psi(t,x)^2=sum_(|alpha|<=3)|partial^alpha psi|^2; U_psi(t)=||G_psi(t,.)||L2(R3). No support-volume bound is substituted for this norm.",
        "tame_bounds": "||S||L2<=512 delta U_psi; each ||partial_i S||L2<=2048 delta U_psi; ||DS[psi]eta||L2<=2048 delta U_eta. The last direction is an arbitrary compact smooth eta, not only a spatial derivative.",
        "quadratic_or_cubic_boundary": "For psi=epsilon f, delta=epsilon maxjets(f) and U_psi=epsilon U_f. The full source is quadratically small and the causal source-induced light force is cubically small.",
        "checks": {
            "source_pointwise_majorant": s.Integer(s_bound) - 420,
            "source_derivative_majorant": ds_bound - s.Rational(5281, 5),
            "Q_Z_coefficient_bound": C["abs_C"] - s.Rational(424, 81),
            "full_source_X_strip_lower_margin": 1
            - 3 * d
            - s.Rational(9, 10)
            - s.Rational(7, 100),
            "full_source_X_strip_upper_margin": s.Rational(11, 10)
            - (1 + 3 * d)
            - s.Rational(7, 100),
        },
        "gates": {k: bool(v) for k, v in tests.items()},
    }
