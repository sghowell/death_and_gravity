"""Literal color/flavor traces and the neutral-background boson blocks."""

from functools import cache

import sympy as sp


@cache
def data():
    I = sp.I
    matrices = [
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / 2,
        sp.Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]]) / 2,
        sp.diag(1, -1, 0) / 2,
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / 2,
        sp.Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]]) / 2,
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / 2,
        sp.Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]]) / 2,
        sp.diag(1, 1, -2) / (2 * sp.sqrt(3)),
    ]
    casimir = sum((T * T for T in matrices), sp.zeros(3))
    y = sp.symbols("reference_Yukawa")
    yukawas = [y, -y] + [sp.Integer(0)] * 12
    color_sum = sum(sp.trace(T * T) for T in matrices)
    checks = {}
    for i, T in enumerate(matrices):
        checks[f"one_current_color_trace_{i}"] = sp.trace(T)
        for j, U in enumerate(matrices):
            checks[f"two_current_color_trace_{i}_{j}"] = sp.simplify(
                sp.trace(T * U) - sp.Rational(i == j, 2)
            )
    for i in range(3):
        for j in range(3):
            checks[f"fundamental_Casimir_{i}_{j}"] = sp.simplify(
                casimir[i, j] - sp.Rational(4, 3) * int(i == j)
            )
    for degree in (1, 3, 5):
        checks[f"opposite_flavor_odd_degree_{degree}"] = sum(v**degree for v in yukawas)
    for degree in (2, 4, 6):
        checks[f"active_even_flavor_degree_{degree}"] = (
            sum(v**degree for v in yukawas) - 2 * y**degree
        )
    checks.update(
        {
            "scalar_active_color_flavor_factor": 2 * sp.trace(sp.eye(3)) - 6,
            "one_flavor_gauge_color_contraction": color_sum - 4,
            "all_fourteen_gauge_vacuum_color_factor": 14 * color_sum - 56,
            "active_gauge_Phi_dependent_color_factor": 2 * color_sum - 8,
            "active_gauge_factor_not_scalar_six": 2 * color_sum - sp.Rational(4, 3) * 6,
        }
    )
    return {
        "fundamental_generators": matrices,
        "fundamental_Casimir_matrix": casimir,
        "Dirac_flavor_Yukawas": yukawas,
        "scalar_active_color_flavor_factor": sp.Integer(6),
        "gauge_active_color_flavor_factor": 2 * color_sum,
        "gauge_all_flavor_vacuum_color_factor": 14 * color_sum,
        "bosonic_block_dictionary": "At A=0 and color-singlet Phi background, the tree scalar and gauge Hessians are block diagonal. The gauge propagator is Phi-independent. The exact H Gaussian is already eliminated; H has no direct Yukawa vertex.",
        "mixed_fermion_block_zero": "The Phi-A Hessian of the fermion determinant has a single color generator and vanishes by Tr T^a=0 for arbitrary color-singlet Phi insertions.",
        "gauge_primitive": "One half Tr[D_A F_AA] includes the fermion-fermion-vector sunset and all its Phi background placements. Gauge vertices carry their original Dirac matrices, factors of i and gauge coupling; the overall trace coefficient is not the complete scalar sign or color coefficient.",
        "ghost_boundary": "There is no direct ghost-fermion vertex and no Phi dependence in the tree ghost block. Pure gauge/ghost vacuum graphs belong to the separate bosonic sector.",
        "regulator_boundary": "Retain the vectorlike gauge-preserving dimensional regulator and the whole closed-loop Hessian, including every ordering. No four-dimensional Lorentz contraction, independent-diagram gauge invariance, or completed gauge/IR estimate is inferred from these finite color identities.",
        "scope": "Nongravitational vacuum matter-sector ownership. Finite-gravity fermion graphs, exact gauge spectrum and the G Regge/IR problem are not covered.",
        "checks": checks,
    }
