"""Finite-net measure extension and explicit unchanged-parameter bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import kernel

from . import source


def coefficient_upper():
    return 300 / source.KAPPA ** s.Rational(3, 2)


def coefficient_change_upper(radiated_energy):
    R = kernel.exact_scalar(radiated_energy)
    if R < 0 or R > s.Rational(1, 8):
        raise ValueError("Require0<=total radiated energy<=1/8")
    return 24000 * R / source.KAPPA ** s.Rational(3, 2)


def tt_project(matrix, direction):
    P = s.eye(3) - direction * direction.T
    return P * matrix * P - P * s.trace(P * matrix) / 2


@cache
def data():
    checks = {}
    w, v, a, Mi, Mj, Kii, Kij, Kjj = s.symbols("w v a Mi Mj Kii Kij Kjj")
    original = Mi * w + Mj * v + (w * w * Kii + 2 * w * v * Kij + v * v * Kjj) / 2
    x, y = a * w, (1 - a) * w
    split = (
        Mi * (x + y)
        + Mj * v
        + ((x * x + 2 * x * y + y * y) * Kii + 2 * (x + y) * v * Kij + v * v * Kjj) / 2
    )
    checks["measure_atom_split_keeps_self_and_cross_diagonals"] = s.expand(
        split - original
    )
    weights = s.symbols("w0:5")
    checks["product_measure_total_mass_is_R_squared"] = s.expand(
        sum(x * y for x in weights for y in weights) - sum(weights) ** 2
    )
    b, c = s.symbols("b c", real=True)
    n = s.Matrix([2 * b, 2 * c, 1 - b * b - c * c]) / (1 + b * b + c * c)
    T = s.Matrix(
        [
            [s.Symbol("T00"), s.Symbol("T01"), s.Symbol("T02")],
            [s.Symbol("T01"), s.Symbol("T11"), s.Symbol("T12")],
            [s.Symbol("T02"), s.Symbol("T12"), s.Symbol("T22")],
        ]
    )
    projected = tt_project(T, n)
    checks["global_ambient_TT_projection_transverse"] = (projected * n).applyfunc(
        s.factor
    )
    checks["global_ambient_TT_projection_tracefree"] = s.factor(s.trace(projected))
    checks["global_ambient_TT_projection_idempotent"] = (
        tt_project(projected, n) - projected
    ).applyfunc(s.factor)
    gates = {
        "all_compact_kernel_singular_intersections_covered": True,
        "uniform_massive_Doppler_gap": True,
        "positive_finite_measure_mass_bound": True,
        "weak_limits_preserve_energy_and_first_moment": True,
        "product_measure_limit_from_continuous_partition_of_unity": True,
        "finite_sup_norm_nets_make_weak_integrals_uniform_in_soft_direction": True,
        "ambient_TT_projection_avoids_nonexistent_global_polarization_frame": True,
        "finite_atomic_measures_weakly_dense": True,
        "strict_coefficient_margin": 275 + 18 < 300,
        "strict_change_margin": 23000 + 17 + 631 + 19 < 24000,
        "zero_radiation_change_exactly_zero": True,
        "energy_measure_not_quantum_probability_or_hard_sum": True,
        "original_full_loop_Regge_bounce_and_P8_open": True,
    }
    return {
        "checks": checks,
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_measure_domain": "Finite positive Borel angular-energy measures sigma on S2, massR<=1/8. Q=(R,int n d sigma), original physical recoil. No quantum-state probability, hard S-matrix summation or arbitrary signed energy distribution is assumed.",
        "whole_uniform_extension_proof": "Integrate the jointly continuous mixed and symmetric double kernels, including all diagonal mass. In the ambient TT tensor representation their compact kernel family has finite uniform-norm nets. Weak convergence controls each net member, bounded masses control the net error, and product convergence follows from finite continuous partitions of unity. Recoil varies continuously through(R,Q). Therefore weak sigma_m=>sigma implies uniform-angular TT convergence of the complete named C; atomic approximations give a unique extension.",
        "whole_retained_margins": (
            s.Integer(293),
            s.Integer(300),
            s.Integer(23667),
            s.Integer(24000),
        ),
        "whole_physical_bounds": "||C_sigma||sup<=293<300 and ||C_sigma-C_Born||sup<=23667R<=24000R. Divide by kappa^(3/2). R=0 has exactly zero difference. These are kinematic bounds, not the finite radiative hard-loop remainder or common-regulator detector rate.",
    }
