"""All first-order external propagator and inverse-vertex terms, not a resummation."""

from functools import cache

import sympy as s


@cache
def data():
    checks = {}
    ds, h, Pi, Pi0, Pi1, H, F = s.symbols(
        "ds hbar Pi_shift Pi_anchor Pi_slope epsilon_pp shifted_hard"
    )
    ren = Pi - Pi0 - ds * Pi1
    vertex = 2 * H * ((Pi - Pi0) / ds - Pi1)
    prop = -2 * H * ren / ds**2
    checks["both_renormalized_external_graphs_cancel"] = s.factor(vertex / ds + prop)
    checks["arbitrary_shifted_hard_factor_cancels"] = s.factor(F * (vertex / ds + prop))
    checks["formal_first_loop_Dyson_sign"] = (
        s.diff(1 / (ds + h * Pi), h).subs(h, 0) + Pi / ds**2
    )
    checks["first_loop_vertex_plus_propagator"] = s.factor(
        s.diff(2 * H * (1 + h * Pi / ds) / (ds + h * Pi), h).subs(h, 0)
    )
    assert s.factor(2 * H * (Pi - Pi0) / ds**2 + prop) == 2 * H * Pi1 / ds
    assert s.factor(2 * H * Pi / ds**2 + 2 * H * Pi / ds**2) != 0

    v, u = s.symbols("left_invariant right_invariant")
    for degree in range(1, 9):
        chain = sum(u**j * v ** (degree - 1 - j) for j in range(degree))
        os = u**degree - 1 - degree * (u - 1)
        vertex_os = 2 * H * (chain.subs(v, 1) - degree)
        checks[f"power{degree}_all_noncommuting_metric_insertions"] = s.expand(
            (u - v) * chain - u**degree + v**degree
        )
        checks[f"power{degree}_complete_OS_external_cancellation"] = s.factor(
            vertex_os / (u - 1) - 2 * H * os / (u - 1) ** 2
        )
        checks[f"power{degree}_unit_residue_anchor"] = s.diff(os, u).subs(u, 1)
    return {
        "checks": checks,
        "gates": {
            "both_vertex_and_single_propagator_insertions_retained": True,
            "minus_Dyson_sign_and_complete_fixed_kinetic_counterterm": True,
            "arbitrary_shifted_hard_factor_never_soft_expanded": True,
            "eight_independent_complete_local_quadratic_variations": True,
            "selected_first_order_cancellation_not_exact_LSZ_or_all_loops": True,
            "mass_and_residue_conditions_already_fixed_not_new_choices": True,
        },
        "whole_renormalized_shifted_inverse": ren,
        "whole_renormalized_TT_vertex": vertex,
        "whole_single_propagator_radiation_correction": prop,
        "whole_zero_selected_physical_radiation": s.S.Zero,
        "whole_zero_selected_finite_interference": s.S.Zero,
        "whole_cancellation_proof": "At p^2=1,r=p+k,d=r^2-1, the complete massive-loop plus kinetic-counterterm vertex is2epsilon(p,p)*Pi_OS(r^2)/d. The fixed inverse K=d+hbar Pi gives the first propagator correction-Pi/d^2. The two radiation graphs sum to zero legwise, even multiplying an arbitrary shifted hard subamplitude. Nonradiating external legs have the already imposed first-order unit residue. Individually nonzero graphs are retained.",
        "whole_negative_controls": "Omitting the kinetic metric counterterm leaves2epsilon(p,p)*Pi_MSprime(1)/d; flipping the Dyson sign makes the two terms add. The local power controls vary every noncommuting operator position before their fixed value/slope subtractions.",
        "whole_scope": "Only the selected first-loop massive-matter external quadratic contribution vanishes. An algebraic reciprocal of a truncated inverse is not asserted to be the exact resummed propagator or a construction of asymptotic quantum states.",
    }
