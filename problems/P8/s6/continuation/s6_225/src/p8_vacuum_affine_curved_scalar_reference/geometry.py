"""Literal curved scalar geometry in the original exponential spatial chart."""

from functools import cache

import sympy as s
from p8_vacuum_affine_local_tensor_response import operator as tensor_input
from p8_vacuum_affine_proca_gaussian import bridge

eta, x, y, z = s.symbols("conformal_time x y z", real=True)
COORDS = (eta, x, y, z)
ed, eg = s.symbols("detector_variation source_variation")
k = s.Symbol("spatial_transfer", positive=True)
q = k * k
a = s.Function("scale", positive=True)(eta)
wd, cd, wg, cg = [s.Function(name)(eta) for name in ("wD", "cD", "wG", "cG")]
pd, pg = s.exp(-s.I * k * x), s.exp(s.I * k * x)
h = s.diff(a, eta) / a
U = s.diff(h, eta) + h * h
L = s.Matrix([[1, -s.Rational(1, 3)], [0, -1]])
Q = q * s.Matrix([[s.Rational(2, 3), 0], [-1, 0]])
C = s.Matrix([[3 * h, -h], [0, 0]])


def jet(value):
    value = s.expand(value)
    return sum(
        value.coeff(ed, i).coeff(eg, j) * ed**i * eg**j for i in (0, 1) for j in (0, 1)
    )


def exponential(value):
    return jet(1 + value + value * value / 2)


def channels(w, c):
    return s.Matrix(
        [
            s.diff(w, eta, 2)
            + 3 * h * s.diff(w, eta)
            + 2 * q * w / 3
            - s.diff(c, eta, 2) / 3
            - h * s.diff(c, eta),
            -q * w - s.diff(c, eta, 2),
        ]
    )


@cache
def literal():
    w = ed * wd * pd + eg * wg * pg
    c = ed * cd * pd + eg * cg * pg
    metric = s.diag(
        a * a,
        -a * a * exponential(2 * (w - c)),
        -a * a * exponential(2 * w),
        -a * a * exponential(2 * w),
    )
    inverse = s.diag(
        1 / (a * a),
        -exponential(-2 * (w - c)) / (a * a),
        -exponential(-2 * w) / (a * a),
        -exponential(-2 * w) / (a * a),
    )
    Gamma = {}
    for aa in range(4):
        for bb in range(4):
            for cc in range(4):
                Gamma[aa, bb, cc] = jet(
                    inverse[aa, aa]
                    * (
                        s.diff(metric[aa, cc], COORDS[bb])
                        + s.diff(metric[aa, bb], COORDS[cc])
                        - s.diff(metric[bb, cc], COORDS[aa])
                    )
                    / 2
                )
    Ric = s.Matrix(
        4,
        4,
        lambda bb, dd: jet(
            sum(
                s.diff(Gamma[aa, dd, bb], COORDS[aa])
                - s.diff(Gamma[aa, aa, bb], COORDS[dd])
                + sum(
                    Gamma[aa, aa, rr] * Gamma[rr, dd, bb]
                    - Gamma[aa, dd, rr] * Gamma[rr, aa, bb]
                    for rr in range(4)
                )
                for aa in range(4)
            )
        ),
    )
    Rold = -jet(sum(inverse[aa, aa] * Ric[aa, aa] for aa in range(4)))
    volume = jet(a**4 * exponential(3 * w - c))
    return {
        "metric": metric,
        "inverse": inverse,
        "Gamma": Gamma,
        "Ric": Ric,
        "Rold": Rold,
        "volume": volume,
    }


def Weyl_component(aa, bb, cc, dd):
    raw = literal()
    metric, Gamma, Ric = raw["metric"], raw["Gamma"], raw["Ric"]
    riem = jet(
        metric[aa, aa]
        * (
            s.diff(Gamma[aa, dd, bb], COORDS[cc])
            - s.diff(Gamma[aa, cc, bb], COORDS[dd])
            + sum(
                Gamma[aa, cc, rr] * Gamma[rr, dd, bb]
                - Gamma[aa, dd, rr] * Gamma[rr, cc, bb]
                for rr in range(4)
            )
        )
    )
    return jet(
        riem
        - (
            metric[aa, cc] * Ric[bb, dd]
            - metric[aa, dd] * Ric[bb, cc]
            - metric[bb, cc] * Ric[aa, dd]
            + metric[bb, dd] * Ric[aa, cc]
        )
        / 2
        - raw["Rold"]
        * (metric[aa, cc] * metric[bb, dd] - metric[aa, dd] * metric[bb, cc])
        / 6
    )


@cache
def profile():
    old = bridge.clock()
    t = bridge.wkb.u
    aa = old["a"]
    hh = s.factor(aa * old["H"])
    dh = s.factor(aa * s.diff(hh, t))
    uu = s.factor(dh + hh * hh)
    clock = t / (2 * (1 + t * t)) + s.atan(t) / 2
    return {
        "t": t,
        "a": aa,
        "h": hh,
        "hprime": dh,
        "U": uu,
        "conformal_primitive": clock,
    }


@cache
def data():
    raw = literal()
    SD, WD = channels(wd, cd)
    SG, WG = channels(wg, cg)
    bg = profile()
    t = bg["t"]
    test = s.Matrix([wg, cg])
    Weyl = Weyl_component(0, 1, 0, 1)
    checks = {
        "literal_exponential_metric_inverse": s.Matrix(
            (raw["metric"] * raw["inverse"] - s.eye(4)).applyfunc(jet)
        ),
        "actual_background_scalar_curvature": s.simplify(
            raw["Rold"].coeff(ed, 0).coeff(eg, 0) - 6 * U / a**2
        ),
        "complete_source_first_scalar_curvature": s.simplify(
            raw["Rold"].coeff(ed, 0).coeff(eg, 1) - 6 * pg * SG / a**2
        ),
        "complete_detector_first_scalar_curvature": s.simplify(
            raw["Rold"].coeff(ed, 1).coeff(eg, 0) - 6 * pd * SD / a**2
        ),
        "curved_source_Weyl_component": s.simplify(
            Weyl.coeff(ed, 0).coeff(eg, 1) - a * a * pg * WG / 3
        ),
        "curved_detector_Weyl_component": s.simplify(
            Weyl.coeff(ed, 1).coeff(eg, 0) - a * a * pd * WD / 3
        ),
        "zero_background_Weyl_component": s.simplify(Weyl.coeff(ed, 0).coeff(eg, 0)),
        "curvature_coordinate_operator_with_first_order_coefficient": channels(wg, cg)
        - L * test.diff(eta, 2)
        - Q * test
        - C * test.diff(eta),
        "actual_scale_factor": bg["a"] - (1 + t * t) ** 2,
        "actual_conformal_Hubble": bg["h"] - 4 * t * (1 + t * t),
        "actual_conformal_Hubble_derivative": bg["hprime"]
        - 4 * (1 + t * t) ** 2 * (1 + 3 * t * t),
        "actual_rescaled_background_curvature": bg["U"]
        - 4 * (1 + t * t) ** 2 * (1 + 7 * t * t),
        "proper_to_conformal_clock": s.factor(
            s.diff(bg["conformal_primitive"], t) - 1 / bg["a"]
        ),
        "original_S189_background_scalar_curvature_bridge": s.factor(
            6 * bg["U"] / bg["a"] ** 2 - tensor_input.R0.subs(tensor_input.t, t)
        ),
        "coefficient_matrix_Frobenius": s.expand(sum(v * v for v in C) - 10 * h * h),
        "coefficient_derivative_Frobenius": s.expand(
            sum(v * v for v in C.diff(eta)) - 10 * s.diff(h, eta) ** 2
        ),
    }
    return {
        "first_curvature_channels": {
            "detector": s.Matrix([SD, WD]),
            "source": s.Matrix([SG, WG]),
        },
        "coordinate_leading_matrix": L,
        "spatial_matrix": Q,
        "first_order_coefficient": C,
        "actual_profile": bg,
        "continuous_coefficient_bounds": {
            "abs_h": s.Rational(5, 2),
            "abs_hprime": s.Rational(175, 16),
            "abs_U": s.Rational(275, 16),
            "C_norm": s.Integer(8),
            "Cprime_norm": s.Integer(35),
            "a4_upper": s.Rational(390625, 65536),
            "conformal_window_upper": s.Integer(1),
        },
        "chart": "Q=2wI-2cPi in the original exponential synchronous spatial chart, with fixed background conformal time. S=a^2 deltaR_old/6 and C0101/(a^2 phase)=W/3. This is not the full lapse/shift/clock reconstruction.",
        "checks": checks,
        "gates": {
            "coefficient_norm_strictly_below8": 10 * s.Rational(5, 2) ** 2 < 8**2,
            "coefficient_derivative_norm_strictly_below35": 10
            * s.Rational(175, 16) ** 2
            < 35**2,
            "background_curvature_scale_below18": s.Rational(275, 16) < 18,
            "physical_conformal_density_below6": s.Rational(390625, 65536) < 6,
            "leading_second_time_derivative_matrix_invertible": L.det() != 0,
            "actual_scale_not_below_one": s.Poly(bg["a"], t).all_coeffs()
            == [1, 0, 2, 0, 1],
        },
    }
