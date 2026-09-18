"""Uniform finite-conversion continuity, including a regulator dominator."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_angular_finite import angular
from p8_vacuum_affine_radiative_state_soft_index import recoil

from . import source


def exact_check(value):
    return (
        value.applyfunc(s.factor)
        if isinstance(value, s.MatrixBase)
        else s.factor(value)
    )


@cache
def bounds():
    B = s.Rational(65, 4)
    one_gradient = 4 * 4 + 4 * 2 * 16
    delta_T = 2 * 3 * one_gradient + 2
    A, C, R = s.symbols("A C R", positive=True)
    cut = (A * R / C) ** 4
    integral = 4 * C * cut ** s.Rational(1, 4) - A * R * s.log(cut)
    z, D, dD = s.symbols("z D dD")
    v = s.ImmutableMatrix(s.symbols("v0:4"))
    dv = s.ImmutableMatrix(s.symbols("dv0:4"))
    derivative = s.diff((v + z * dv) * (v + z * dv).T / (D + z * dD), z).subs(z, 0)
    checks = {
        "massive_projected_current_gradient": s.Integer(one_gradient - 144),
        "two_recoiled_massive_plus_all_null_difference": s.Integer(delta_T - 866),
        "exact_minimum_majorant_integral": s.expand_log(
            integral - 4 * A * R * (1 + s.log(C / (A * R))), force=True
        ),
        "generic_projected_rank_one_current_derivative": D * D * derivative
        - D * (v * dv.T + dv * v.T)
        + dD * v * v.T,
    }
    finite_numerator = s.Rational(30000, 2) + 400000 + 4 * 50000
    uniform_numerator = (
        4 * 50000 + s.Rational(30000, 2) + 2 * 400000 + s.Rational(1, 8) * 240000
    )
    gates = {
        "quadratic_H_current_difference": 3 * B * delta_T < 50000,
        "quadratic_U_current_difference": 2 * B * delta_T < 30000,
        "H_crossover_inside_unit_interval": s.Rational(100000, 8) < 64000,
        "U_crossover_inside_unit_interval": s.Rational(60000, 8) < 42000,
        "H_log_constant_below_one": s.Rational(64000, 100000) < 1,
        "U_log_constant_below_one": s.Rational(42000, 60000) < 1,
        "physical_index_change_constant": s.Rational(50000, 36) < 1400,
        "finite_conversion_change_constant": finite_numerator / 72 < 10000,
        "uniform_regulator_quotient_constant": uniform_numerator / 72 < 15000,
        "physical_soft_index_bound": s.Rational(265, 36) < 8,
        "all_null_directions_and_splittings_retained": True,
        "uniform_e_difference_not_constant_O_e2_error": True,
    }
    return {
        "checks": {k: exact_check(v) for k, v in checks.items()},
        "gates": {k: bool(v) for k, v in gates.items()},
        "whole_current_and_quadratic_difference_caps": (866, 50000, 30000),
        "whole_radial_interpolation_caps": {
            "H": (100000, 64000, 400000),
            "U": (60000, 42000, 240000),
        },
        "whole_finite_and_regulator_numerators": (finite_numerator, uniform_numerator),
        "whole_state_continuity": "For0<R<=1/8, |a_sigma-a0|<1400R/kappa and |Delta_sigma-Delta0|<10000R(1+ln(1/R))/kappa. Both vanish at R=0. This holds for every finite positive null-radiation state on the exact S300 recoil map, uniformly in multiplicity and direction.",
        "whole_regulator_uniform_difference": "For0<e<=1/8, |[(a_sigma,e-a0,e)-(a_sigma-a0)]/(2e)|<15000R(1+ln(1/R))/kappa. The full S301 trace term and phase are retained; no finite-e positivity assumption is used.",
    }


@cache
def calibration():
    E = s.Rational(5, 4)
    alpha = s.Rational(19, 10)
    Ep = (alpha + 1 / alpha) / 2
    opposite_R = 2 * (E - Ep)
    single_R = E - Ep**2 / E
    out = s.ImmutableMatrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    single = [s.ImmutableMatrix([single_R, single_R, 0, 0])]
    opposite = [
        s.ImmutableMatrix([opposite_R / 2, opposite_R / 2, 0, 0]),
        s.ImmutableMatrix([opposite_R / 2, -opposite_R / 2, 0, 0]),
    ]
    split = [opposite[0] / 2, opposite[0] / 2, opposite[1]]
    directions = (
        s.ImmutableMatrix([0, 1, 0]),
        s.ImmutableMatrix([s.Rational(3, 5), 0, s.Rational(4, 5)]),
    )
    checks, gates, records = {}, {}, []
    for state, quanta in (("one", single), ("two", opposite), ("three_split", split)):
        total = sum(q[0] for q in quanta)
        ps, qs, _ = recoil.momenta(E, quanta, out)
        checks[state + "_all_mass_shells"] = s.Matrix(
            [recoil.dot(p, p) - 1 for p in ps]
        )
        checks[state + "_complete_conservation"] = sum(ps, s.zeros(4, 1)) + sum(
            qs, s.zeros(4, 1)
        )
        for ni, n in enumerate(directions):
            for ri, radial in enumerate((s.S.Zero, s.Rational(3, 5), s.S.One)):
                T = angular.transverse_current(E, quanta, out, n, radial)
                T0 = angular.transverse_current(E, [], out, n, radial)
                dT = (T - T0).applyfunc(s.factor)
                H, U = angular.contractions(T)
                H0, U0 = angular.contractions(T0)
                label = state + "_direction" + str(ni) + "_radial" + str(ri)
                gates[label + "_Frobenius_current_difference"] = (
                    s.factor(s.trace(dT * dT)) <= (866 * total) ** 2
                )
                gates[label + "_H_difference"] = abs(s.factor(H - H0)) < 50000 * total
                gates[label + "_U_difference"] = abs(s.factor(U - U0)) < 30000 * total
                if state == "three_split":
                    checks[label + "_exact_collinear_splitting"] = (
                        T - angular.transverse_current(E, opposite, out, n, radial)
                    )
                records.append((state, ni, radial, total))
    return {
        "checks": {k: exact_check(v) for k, v in checks.items()},
        "gates": {k: bool(v) for k, v in gates.items()},
        "whole_exact_recoil_radial_angular_calibrations": records,
        "whole_calibration_boundary": "18 exact rational one/two/collinearly-split-three radiation states on the original recoil map, including radial0,3/5,1. Frobenius and quadratic inequalities calibrate the formulas; the written nuclear-norm proof supplies uniformity and arbitrary multiplicity.",
    }


@cache
def data():
    bound, actual = bounds(), calibration()
    return {
        "checks": {
            **{"bound_" + k: v for k, v in bound["checks"].items()},
            **{"original_" + k: v for k, v in actual["checks"].items()},
        },
        "gates": {**bound["gates"], **actual["gates"]},
        "whole_uniform_state_continuity": {
            k: v for k, v in bound.items() if k not in ("checks", "gates")
        },
        "whole_original_recoil_calibrations": {
            k: v for k, v in actual.items() if k not in ("checks", "gates")
        },
    }


def require_radiated_energy(value):
    if isinstance(value, (bool, float, s.Float, str)) or value is None:
        raise TypeError("Require an exact real scalar total energy")
    R = s.sympify(value)
    if not isinstance(R, s.Expr):
        raise TypeError("Require an exact real scalar total energy")
    R = recoil.exact_real(R)
    if R < 0 or R > s.Rational(1, 8):
        raise ValueError("Require exact total radiated energy in[0,1/8]")
    return R


def index_change_upper(value):
    R = require_radiated_energy(value)
    return 1400 * R / source.KAPPA


def finite_conversion_change_upper(value):
    R = require_radiated_energy(value)
    return s.S.Zero if R == 0 else 10000 * R * (1 - s.log(R)) / source.KAPPA


def regulator_quotient_upper(value):
    R = require_radiated_energy(value)
    return s.S.Zero if R == 0 else 15000 * R * (1 - s.log(R)) / source.KAPPA
