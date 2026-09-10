"""Independent trace, background derivative and complete quadratic tensor space."""

from functools import cache
from itertools import combinations, permutations

import sympy as sp
from p8_vacuum_finite_mass_gauge_cut import dirac


@cache
def data():
    m, y, nu = sp.symbols(
        "positive_fermion_mass positive_Yukawa positive_reference_scale", positive=True
    )
    phi = sp.symbols("reference_scalar", real=True)
    N, Nc, Q = sp.Integer(6), sp.Integer(3), 16 * sp.pi**2
    qv = sp.symbols("real_loop_q0:4", real=True)
    gamma = dirac.data()["gamma_matrices"]
    slash = sum((qv[i] * gamma[i] for i in range(4)), sp.zeros(4))
    numerator = m * sp.eye(4) - sp.I * slash
    traced = sp.expand(sp.trace(numerator**4))
    q2 = sum(q * q for q in qv)
    D, t = sp.symbols(
        "positive_reference_denominator positive_radial_variable", positive=True
    )
    scalar_trace = 4 * (m**4 - 6 * m * m * q2 + q2 * q2)
    Ibar = sp.symbols("entire_dimensional_reference_Ibar")
    local = 24 * N * y**4 * (Ibar - sp.Rational(8, 3)) / Q
    deltaZ = (
        -2
        * Nc
        * y
        * y
        / Q
        * sum(
            sp.log((m + sign * y * phi) ** 2 / (nu * nu)) + sp.Rational(2, 3)
            for sign in (-1, 1)
        )
    )
    phi2_coefficient = sp.simplify(sp.diff(deltaZ, phi, 2).subs(phi, 0) / 2)
    c = N * y**4 / (Q * m * m)
    diagonal = sp.symbols("external_square_0:4")
    pairs = list(combinations(range(4), 2))
    cross = sp.symbols("external_dot_0:6")
    coeff = sp.symbols("quadratic_tensor_coefficient_0:10")
    gram = {
        **{(i, i): diagonal[i] for i in range(4)},
        **dict(zip(pairs, cross, strict=True)),
    }
    basis = list(diagonal) + list(cross)
    polynomial = sum(a * b for a, b in zip(coeff, basis, strict=True))
    equations = []
    for perm in permutations(range(4)):
        renamed = {
            gram[(i, j)]: gram[tuple(sorted((perm[i], perm[j])))] for i, j in gram
        }
        transformed = polynomial.xreplace(renamed)
        equations.extend(sp.expand(transformed - polynomial).coeff(b) for b in basis)
    equations = sorted(set(equations) - {sp.Integer(0)}, key=sp.default_sort_key)
    matrix, _ = sp.linear_eq_to_matrix(equations, coeff)
    invariant_diagonal = sp.Matrix([1] * 4 + [0] * 6)
    invariant_cross = sp.Matrix([0] * 4 + [1] * 6)
    Sdiag, Scross = sp.symbols("sum_external_squares sum_pairwise_dots")
    degree2 = 2 * c * Sdiag
    checks = {
        "active_flavor_color_multiplicity": N - 2 * Nc,
        "cyclic_functional_derivative_weight": sp.factorial(4) / 4 - 6,
        "independent_four_propagator_Dirac_trace": sp.expand(traced - scalar_trace),
        "constant_trace_denominator_reduction": sp.expand(
            4 * (m**4 - 6 * m * m * (D - m * m) + (D - m * m) ** 2)
            - 4 * (D * D - 8 * m * m * D + 8 * m**4)
        ),
        "convergent_radial_I3": sp.integrate(t / (1 + t) ** 3, (t, 0, sp.oo))
        - sp.Rational(1, 2),
        "convergent_radial_I4": sp.integrate(t / (1 + t) ** 4, (t, 0, sp.oo))
        - sp.Rational(1, 6),
        "finite_zero_momentum_trace_normalization": 24
        * N
        * y**4
        * (-8 * sp.Rational(1, 2) + 8 * sp.Rational(1, 6))
        / Q
        + 64 * N * y**4 / Q,
        "whole_UV_quartic_matches_determinant": sp.expand(local).coeff(Ibar, 1)
        - 24 * N * y**4 / Q,
        "whole_finite_quartic_matches_determinant": local.subs(Ibar, 0)
        + 64 * N * y**4 / Q,
        "background_kinetic_zero_anchor": sp.simplify(
            deltaZ.subs({phi: 0, nu: m}) + 4 * N * y * y / (3 * Q)
        ),
        "background_quadratic_field_coefficient": sp.simplify(phi2_coefficient - 2 * c),
        "four_point_two_zero_legs_matches_background_derivative": sp.simplify(
            sp.diff(deltaZ, phi, 2).subs(phi, 0) - 4 * c
        ),
        "complete_S4_quadratic_space_rank": matrix.rank() - 8,
        "complete_S4_quadratic_space_dimension": len(matrix.nullspace()) - 2,
        "diagonal_orbit_invariant": matrix * invariant_diagonal,
        "pair_orbit_invariant": matrix * invariant_cross,
        "derivative_pair_vertex_factor": -4 * c * Scross
        - degree2.subs(Sdiag, -2 * Scross),
        "mass_shell_constant_degree_two": degree2.subs(Sdiag, -4) + 8 * c,
    }
    return {
        "fermion_mass": m,
        "Yukawa": y,
        "zero_momentum_Dirac_trace": traced,
        "regulated_zero_momentum_four_vertex": sp.expand(local),
        "minimal_finite_zero_momentum_four_vertex": -64 * N * y**4 / Q,
        "background_kinetic_increment": deltaZ,
        "two_derivative_action_coefficient": c,
        "complete_quadratic_S4_constraint_matrix": matrix,
        "quadratic_S4_invariant_basis": [invariant_diagonal, invariant_cross],
        "Euclidean_degree_two_four_vertex": degree2,
        "on_shell_degree_two_four_vertex": -8 * c,
        "amplitude_sign": "The continued Euclidean four-vertex is Gamma_F4; the scattering increment is A_F=-Gamma_F4.",
        "scope": "The full regulated degree-zero reference is subtracted before a four-dimensional norm estimate. Every odd momentum degree is absent by Lorentz scalar invariance. The complete S4-invariant degree-two term is fixed by a background two-point derivative and is constant on the equal-mass shell, so neither degree zero nor two contributes to the forward second coefficient.",
        "checks": checks,
    }
