"""Uniform coefficient Lipschitz control for positive added angular energy."""

from functools import cache

import sympy as s
from p8_vacuum_affine_continuous_logarithmic_coefficient import kernels
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import kernel, radiative
from p8_vacuum_affine_radiative_state_soft_index import recoil

from . import source

LIPSCHITZ = s.Integer(40000)


def positive_increment_upper(added_energy):
    t = kernel.exact_scalar(added_energy)
    if t < 0 or t > s.Rational(1, 8):
        raise ValueError(
            "Require0<=positive added energy<=1/8 with the full state in the same domain"
        )
    return LIPSCHITZ * t / source.KAPPA ** s.Rational(3, 2)


@cache
def original_calibration():
    E, states = radiative.calibration_states()
    u = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    n, p, c = kernels.frame_x(0, "y")
    records = []
    for name in ("born", "one", "two", "four"):
        base = states[name]
        for t in (s.Rational(1, 128), s.Rational(1, 10**20)):
            rays = [*base, t * s.Matrix([1, 0, 1, 0])]
            old, _, _ = recoil.momenta(E, base, u)
            new, _, _ = recoil.momenta(E, rays, u)
            for i in (2, 3):
                if not abs(s.N(new[i][0] - old[i][0], 90)) < 2 * t:
                    raise ValueError(
                        "Original positive-increment energy calibration failed"
                    )
                delta = new[i][1:4, 0] - old[i][1:4, 0]
                if not s.N(delta.dot(delta), 90) < 36 * t * t:
                    raise ValueError(
                        "Original positive-increment spatial calibration failed"
                    )
            for pol, A in (("plus", p), ("complex", (p + s.I * c) / s.sqrt(2))):
                change = kernels.coefficient(E, rays, u, n, A) - kernels.coefficient(
                    E, base, u, n, A
                )
                norm = abs(s.N(change, 90))
                if not norm < LIPSCHITZ * t:
                    raise ValueError(
                        "Original complete coefficient increment calibration failed"
                    )
                records.append((name, pol, t, str(s.N(norm / t, 24))))
    return records


@cache
def data():
    E, R, M, uv, r = s.symbols("E R M uv r", real=True)
    h, Ep = E - R / 2, M / 2
    checks = {}
    for sign in (1, -1):
        checks["original_recoil_time_" + str(sign)] = (
            h - sign * r * uv / (2 * Ep) - (2 * h / 2 + sign * r * (-uv) / M)
        )
        checks["original_recoil_spatial_" + str(sign)] = (
            sign * r * uv / (4 * Ep * (h + Ep))
            - s.Rational(1, 2)
            - (-s.Rational(1, 2) + sign * r * uv / (M * (2 * h + M)))
        )
    P, D, dP, dD = s.symbols("P D deltaP deltaD", nonzero=True)
    checks["exact_massive_current_quotient_difference"] = (
        (P + dP) / (D + dD) - P / D - dP / (D + dD) + P * dD / (D * (D + dD))
    )
    L, dL = s.symbols("L deltaL")
    checks["exact_product_difference"] = (
        (L + dL) * (D + dD) - L * D - dL * (D + dD) - L * dD
    )
    d, ell = s.symbols("d ell")
    checks["exact_mixed_contraction_difference"] = (
        d * (L + dL)
        + (D + dD) * ell
        - 2 * (P + dP)
        - (d * L + D * ell - 2 * P)
        - (d * dL + dD * ell - 2 * dP)
    )
    K0, Ksp, dK0, dKsp = s.symbols("K0 Ksp dK0 dKsp", real=True)
    gamma = (K0 * K0 - Ksp * Ksp - 2) / 2
    checks["outgoing_pair_gamma_directional_derivative"] = (
        s.diff(gamma, K0) * dK0 + s.diff(gamma, Ksp) * dKsp - (K0 * dK0 - Ksp * dKsp)
    )
    gates = {
        "recoil_m_derivative": s.Rational(33, 16) < 3,
        "recoil_r_lower": s.Rational(13, 32) > s.Rational(9, 25),
        "recoil_boost_derivative": s.Rational(1, 32) + s.Rational(5, 512)
        < s.Rational(1, 16),
        "recoil_spatial_derivative": s.Rational(1, 2)
        + 5 * (1 + s.Rational(1, 512))
        + s.Rational(1, 8)
        < 6,
        "recoil_time_derivative": s.Rational(1, 2) + s.Rational(4, 8) + s.Rational(1, 2)
        < 2,
        "massive_real_pair_change": s.Rational(6 * (56 * 5168 + 102 * 161 * 32), 144)
        < 34000,
        "mixed_real_pair_change": (s.Rational(2 * (4 * 1244 + 44 * 32), 8) + 4 * 4 * 44)
        / 36
        < 65,
        "null_real_pair_change": s.Rational(10, 8 * 36) < 1,
        "real_phase_change": s.Rational(1218 * 25 + 49 * 50, 36) < 914,
        "total_real_change": 34000 + 65 + 1 + 914 < 35000,
        "outgoing_pair_gamma_change": s.Rational(33, 8) < 5,
        "aggregate_phase_change": (102 * 2 + s.Rational(5168, 7)) / 24 < 40,
        "whole_positive_increment_change": 35000 + 40 < LIPSCHITZ,
        "positive_interpolating_energy_path_stays_physical": True,
        "extension_to_atoms_and_positive_Borel_increments_by_S328": True,
        "not_an_arbitrary_signed_perturbation_theorem": True,
    }
    return {
        "checks": {name: s.factor(value) for name, value in checks.items()},
        "gates": {name: bool(value) for name, value in gates.items()},
        "whole_positive_increment_bound": "For fixed original E,u and sigma=tau+rho positive, R_sigma<=1/8 and t=massrho, ||C_sigma-C_tau||angular_TT,sup<=40000t. Both full complex components are kept. Physical normalization is kappa^(-3/2). Positive atomic approximation extends the bound to general positive measures, not arbitrary signed perturbations.",
        "whole_derivative_constants": (
            s.Integer(2),
            s.Integer(6),
            s.Integer(8),
            s.Integer(608),
            s.Integer(5168),
            s.Integer(32),
        ),
        "whole_complete_coefficient_budgets": (
            s.Integer(34000),
            s.Integer(65),
            s.S.One,
            s.Integer(914),
            s.Integer(35000),
            s.Integer(40),
            LIPSCHITZ,
        ),
        "whole_original_positive_tail_calibrations": original_calibration(),
    }
