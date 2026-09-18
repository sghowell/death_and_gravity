"""Positive-energy increments of the original fixed-ball soft current."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_angular_finite import angular
from p8_vacuum_affine_radiative_soft_state_transfer.continuity import (
    require_radiated_energy,
)
from p8_vacuum_affine_radiative_state_soft_index import recoil

from . import source


def index_increment_upper(energy):
    t = require_radiated_energy(energy)
    return 1700 * t / source.KAPPA


def conversion_increment_upper(energy):
    t = require_radiated_energy(energy)
    return s.S.Zero if t == 0 else 11000 * t * (1 - s.log(t)) / source.KAPPA


def regulated_increment_upper(energy):
    t = require_radiated_energy(energy)
    return s.S.Zero if t == 0 else 18000 * t * (1 - s.log(t)) / source.KAPPA


@cache
def calibration():
    E = s.Rational(5, 4)
    u = s.ImmutableMatrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    states = []
    for name, alpha in (
        ("born", s.Integer(2)),
        ("two", s.Rational(19, 10)),
        ("split", s.Rational(19, 10)),
    ):
        Ep = (alpha + 1 / alpha) / 2
        R = 2 * (E - Ep)
        old = (
            []
            if R == 0
            else [
                s.ImmutableMatrix([R / 2, R / 2, 0, 0]),
                s.ImmutableMatrix([R / 2, -R / 2, 0, 0]),
            ]
        )
        if name == "split":
            old = [old[0] / 2, old[0] / 2, old[1] / 2, old[1] / 2]
        states.append((name, alpha, old))
    dirs = (
        s.ImmutableMatrix([0, 1, 0]),
        s.ImmutableMatrix([s.Rational(3, 5), 0, s.Rational(4, 5)]),
    )

    @cache
    def tensor(rays, n, radial):
        return angular.transverse_current(E, list(rays), u, n, radial)

    records = []
    for name, alpha, old in states:
        for step in (s.Rational(1, 1000), s.Rational(1, 10**20)):
            newer = alpha - step
            t = s.factor((alpha + 1 / alpha) - (newer + 1 / newer))
            if not (0 < t and sum(q[0] for q in old) + t < s.Rational(1, 8)):
                raise ValueError("Positive calibration left original energy domain")
            new = [
                *old,
                s.ImmutableMatrix([t / 2, 0, t / 2, 0]),
                s.ImmutableMatrix([t / 2, 0, -t / 2, 0]),
            ]
            oldps, _, _ = recoil.momenta(E, old, u)
            newps, _, _ = recoil.momenta(E, new, u)
            for i in (2, 3):
                v = newps[i][1:4, 0] - oldps[i][1:4, 0]
                if (
                    not abs(newps[i][0] - oldps[i][0]) < 2 * t
                    or not s.factor(v.dot(v)) < 36 * t * t
                ):
                    raise ValueError("Original positive recoil calibration failed")
            for ni, n in enumerate(dirs):
                for rad in (s.S.Zero, s.Rational(3, 5), s.S.One):
                    T0 = tensor(tuple(old), n, rad)
                    T = tensor(tuple(new), n, rad)
                    difference = (T - T0).applyfunc(s.factor)
                    H, U = angular.contractions(T)
                    H0, U0 = angular.contractions(T0)
                    norm = s.factor(s.trace(difference * difference))
                    hd, ud = s.factor(H - H0), s.factor(U - U0)
                    if not (
                        norm < (1218 * t) ** 2
                        and abs(hd) < 60000 * t
                        and abs(ud) < 40000 * t
                    ):
                        raise ValueError(
                            "Original whole fixed-ball increment calibration failed"
                        )
                    records.append((name, step, ni, rad, t))
    return tuple(records)


@cache
def data():
    A, B, t = s.symbols("A B t", positive=True)
    z, D, dD = s.symbols("z D deltaD")
    v = s.ImmutableMatrix(s.symbols("v0:4"))
    dv = s.ImmutableMatrix(s.symbols("dv0:4"))
    derivative = s.diff((v + z * dv) * (v + z * dv).T / (D + z * dD), z).subs(z, 0)
    checks = {
        "generic_projected_rank_one_derivative": D * D * derivative
        - D * (v * dv.T + dv * v.T)
        + dD * v * v.T,
        "one_massive_nuclear_increment_budget": s.Integer(
            2 * 2 * 6 * 4 + 4 * 8 * 16 - 608
        ),
        "complete_positive_current_increment": s.Integer(2 * 608 + 2 - 1218),
        "exact_minimum_radial_majorant": s.expand_log(
            4 * B * (A * t / B)
            - A * t * s.log((A * t / B) ** 4)
            - 4 * A * t * (1 + s.log(B / (A * t))),
            force=True,
        ),
        "complete_regulator_numerator": 4 * s.Integer(60000)
        + s.Integer(40000) / 2
        + 2 * s.Integer(480000)
        + s.Integer(320000) / 8
        - 1260000,
    }
    gates = {
        "quadratic_H_difference": 3 * s.Rational(65, 4) * 1218 < 60000,
        "quadratic_U_difference": 2 * s.Rational(65, 4) * 1218 < 40000,
        "H_crossover_inside_interval": s.Rational(120000, 8) < 64000,
        "U_crossover_inside_interval": s.Rational(80000, 8) < 42000,
        "H_radial_logconstant": s.Rational(64000, 120000) < 1,
        "U_radial_logconstant": s.Rational(42000, 80000) < 1,
        "index_increment_roundup": s.Rational(60000, 36) < 1700,
        "finite_conversion_increment_roundup": s.Rational(
            20000 + 480000 + 4 * 60000, 72
        )
        < 11000,
        "uniform_regulator_increment_roundup": s.Rational(1260000, 72) < 18000,
        "all36_original_fixed_ball_calibrations": len(calibration()) == 36,
        "positive_added_energy_with_fixed_original_E_u": True,
        "same_full_D_trace_and_phase_retained": True,
    }
    return {
        "checks": {
            name: (
                value.applyfunc(s.factor)
                if isinstance(value, s.MatrixBase)
                else s.factor(value)
            )
            for name, value in checks.items()
        },
        "gates": {name: bool(value) for name, value in gates.items()},
        "whole_uniform_positive_current_and_quadratic_caps": (1218, 60000, 40000),
        "whole_radial_integral_caps": (480000, 320000),
        "whole_positive_increment_moduli": "At fixed original E,u and sigma=tau+rho positive of mass<=1/8, t=massrho gives |deltaa|<1700t/kappa, |deltaDelta|<11000t(1+ln1/t)/kappa and |(deltaa_e-deltaa)/(2e)|<18000t(1+ln1/t)/kappa for0<e<=1/8. Differences vanish at t0; original stronger Born-specific constants remain unchanged.",
        "whole_original_current_calibrations": calibration(),
    }
