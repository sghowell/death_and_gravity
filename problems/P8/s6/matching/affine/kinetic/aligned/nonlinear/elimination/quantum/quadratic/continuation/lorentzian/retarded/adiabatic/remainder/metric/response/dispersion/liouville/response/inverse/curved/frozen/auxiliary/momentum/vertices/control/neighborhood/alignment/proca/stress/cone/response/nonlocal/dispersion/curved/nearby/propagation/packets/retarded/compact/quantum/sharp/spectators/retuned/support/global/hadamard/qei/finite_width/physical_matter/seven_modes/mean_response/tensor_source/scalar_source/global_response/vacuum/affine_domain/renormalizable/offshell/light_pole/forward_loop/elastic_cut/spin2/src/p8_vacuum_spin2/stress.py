"""Literal canonical stress insertion and projected local/loop tensor terms."""

from functools import cache

import sympy as sp


@cache
def data():
    eta = sp.diag(1, -1, -1, -1)
    k = sp.Matrix(sp.symbols("internal_k0:4", real=True))
    Q = sp.Symbol("transfer_spatial", real=True)
    q = sp.Matrix([0, 0, 0, Q])
    kp = k + q
    e = sp.Matrix([1, 1, 0, 0])
    mass = sp.Symbol("line_mass_squared", positive=True)
    local, xi = sp.symbols(
        "local_potential_vertex curvature_improvement_vertex", real=True
    )
    dot = lambda a, b: (a.T * eta * b)[0]
    stress = k * kp.T + kp * k.T - eta * (dot(k, kp) - mass)
    lower_e = eta * e
    project = lambda T: (lower_e.T * T * lower_e)[0]
    eq = dot(e, q)
    e2 = dot(e, e)
    ep, el, chi, z, L2, n = sp.symbols(
        "e_dot_external e_dot_loop single_parameter pair_parameter radial_loop_square dimension",
        real=True,
    )
    shifted_numerator = 2 * (el + chi * ep - z * eq) ** 2
    averaged = sp.expand(shifted_numerator).subs({el**2: L2 * e2 / n, el: 0})
    # A local quartic interaction supplies no individual external P vector
    # to its stress bubble. Only loop metric and q tensors survive averaging.
    quartic_numerator = 2 * (el - z * eq) ** 2
    quartic_averaged = sp.expand(quartic_numerator).subs({el**2: L2 * e2 / n, el: 0})
    return {
        "canonical_two_scalar_stress_vertex": stress,
        "null_projection_vector": e,
        "projected_stress_vertex": project(stress),
        "projected_shifted_triangle_numerator": averaged,
        "projected_local_quartic_stress_bubble": quartic_averaged,
        "diagram_boundary": "At this matter one-loop order, the two kinetic stress triangles are the only spin-two contributions. Potential, mass, curved R times scalar-polynomial improvements and the local quartic stress bubble have only metric/transfer tensor structures after integration. The existing light kinetic counterterm supplies the constant Ward normalization.",
        "checks": {
            "explicit_projection_vector_is_null": e2,
            "explicit_projection_vector_is_transfer_transverse": eq,
            "literal_massive_kinetic_stress_projection": sp.expand(
                project(stress) - 2 * dot(e, k) ** 2
            ),
            "metric_potential_and_mass_terms_project_to_zero": project(local * eta),
            "general_scalar_curvature_improvement_projects_to_zero": project(
                xi * (q * q.T - eta * dot(q, q))
            ),
            "actual_triangle_single_parameter_squared_numerator": sp.expand(
                averaged - 2 * chi * chi * ep * ep
            ),
            "local_quartic_stress_bubble_has_no_spin_two_external_weight": quartic_averaged,
        },
    }
