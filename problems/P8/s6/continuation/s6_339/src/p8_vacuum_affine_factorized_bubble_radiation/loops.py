"""Both exact-D light insertions and the complete covariant A-B-A chain."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_four_point_loop import amplitude
from p8_vacuum_affine_mixed_source_radiation import triangles

from . import source


def loop_prefactor():
    return source.CUBIC**4 / (32 * s.pi**2 * source.HEAVY_MASS2**2)


@cache
def data():
    z, y, a, e, v = s.symbols("z y pk epsilon p_squared", positive=True)
    m = s.Symbol("heavy_mass_squared", positive=True)
    M = 1 - z * (1 - z) * v
    ML = M - 2 * z * y * a
    MH = M - 2 * (1 - z) * y * a
    MR = M - 2 * z * (1 - z) * a
    prior = triangles.data()
    FL = -z * ML ** (-e) / (a * e)
    FH = -(1 - z) * MH ** (-e) / (a * e)
    checks = {
        "literal_equal_mass_light_triangle": s.expand(
            prior["whole_light_denominator"].subs(m, 1) - ML
        ),
        "literal_equal_mass_second_triangle": s.expand(
            prior["whole_heavy_denominator"].subs(m, 1) - MH
        ),
        "coincident_whole_D_sum_of_parameter_widths": s.expand(
            z * z * (1 - z) + (1 - z) ** 2 * z - z * (1 - z)
        ),
        "first_whole_D_antiderivative": s.simplify(
            s.diff(FL, y) + 2 * z * z * ML ** (-1 - e)
        ),
        "second_whole_D_antiderivative": s.simplify(
            s.diff(FH, y) + 2 * (1 - z) ** 2 * MH ** (-1 - e)
        ),
        "complete_whole_D_bubble_divided_difference": s.simplify(
            FL.subs(y, 1 - z)
            - FL.subs(y, 0)
            + FH.subs(y, z)
            - FH.subs(y, 0)
            + (MR ** (-e) - M ** (-e)) / (a * e)
        ),
        "coincident_exact_D_bubble_derivative": s.simplify(
            s.diff(s.gamma(e) * M ** (-e), v)
            - s.gamma(e + 1) * z * (1 - z) * M ** (-e - 1)
        ),
    }
    n, g, C = amplitude.n, amplitude.g, amplitude.C
    v, u = s.symbols("v u")
    Bl, Br, Delta, Bp = s.symbols("B_left B_right Delta B_prime")
    Al = C + g * g / (n - v)
    Ar = C + g * g / (n - u)
    DA = g * g / ((n - v) * (n - u))
    DB = (Br - Bl) / (u - v)
    checks.update(
        {
            "both_outer_heavy_resolvent_insertions": s.factor((Ar - Al) / (u - v) - DA),
            "complete_ordered_A_B_A_covariant_variation": s.factor(
                DA * Br * Ar
                + Al * DB * Ar
                + Al * Bl * DA
                - (Ar * Ar * Br - Al * Al * Bl) / (u - v)
            ),
            "coincident_full_kernel_derivative": s.factor(
                s.cancel((Ar * Ar * (Bl + (u - v) * Bp) - Al * Al * Bl) / (u - v)).subs(
                    u, v
                )
                - 2 * Al * g * g * Bl / (n - v) ** 2
                - Al * Al * Bp
            ),
        }
    )
    frozen = amplitude.complete_loop().replace(
        lambda x: x.func in (amplitude.Cfun, amplitude.Dfun), lambda x: s.S.Zero
    )
    expected = sum(
        amplitude.fixed_vertex(q) ** 2 * amplitude.Bfun(q)
        for q in (amplitude.svar, amplitude.tvar, amplitude.uvar)
    ) / (32 * s.pi**2)
    checks["entire_frozen_three_channel_bubble_coefficient"] = s.factor(
        frozen - expected
    )
    dc = -3 * C * C * Delta / (32 * s.pi**2)
    dg = -C * g * Delta / (32 * s.pi**2)
    dm = g * g * Delta / (32 * s.pi**2)
    counter = dc / 3 + 2 * g * dg / (n - v) - g * g * dm / (n - v) ** 2
    checks["all_three_original_UV_counterterms"] = s.factor(
        counter + Al * Al * Delta / (32 * s.pi**2)
    )
    checks["full_radiative_UV_counterterm_divided_difference"] = s.factor(
        (counter.subs(v, u) - counter) / (u - v)
        + Delta * (Ar * Ar - Al * Al) / (32 * s.pi**2 * (u - v))
    )
    B0 = s.Symbol("symmetric_bubble")
    symmetric_vertex = C + g * g / (n - s.Rational(4, 3))
    return {
        "checks": checks,
        "gates": {
            "both_internal_light_lines_and_both_outer_heavy_insertions": True,
            "covariant_operator_order_retained_before_metric_variation": True,
            "exact_D_before_subtraction_and_physical_continuation": True,
            "all_contact_trilinear_heavy_mass_UV_counterterms": True,
            "existing_full_parent_OS4_condition_used_by_linearity": True,
            "other_Hessian_classes_and_quadratic_blocks_separate": True,
        },
        "whole_inner_light_denominator": ML,
        "whole_inner_second_denominator": MH,
        "whole_inner_antiderivatives": (FL, FH),
        "whole_generic_endpoint_vertex": Al,
        "whole_flat_selected_bubble_coefficient": expected,
        "whole_generic_UV_counterterms": {"delta_C": dc, "delta_g": dg, "delta_n": dm},
        "whole_existing_linear_finite_contact": -3
        * symmetric_vertex**2
        * B0
        / (32 * s.pi**2),
        "whole_original_normalized_prefactor": loop_prefactor(),
        "whole_covariant_graph_proof": "The local/local Hessian term is the covariant A-B-A chain acting on Phi^2. Both light-line graviton insertions give D_B in exact D, both heavy-line insertions give D_A, and the ordered product identity gives D_(A^2B). All external/source contacts remain. The heavy denominator is never expanded in a loop integral. The pole is cancelled by all three existing local counterterms after the same covariant variation.",
        "whole_branch_boundary": "Specializing the two-mass triangle algebra to equal mass does not import the earlier large-heavy-mass positive gap. Prove the identity on a subthreshold analytic domain and continue to the same Feynman cut boundary, using the separate physical branch estimates. The continuous coincidence limit and the full imaginary part are retained.",
        "whole_scope": "This complete factorized bubble class is not the local/bilocal triangles, six ordered boxes, unresolved mixed quadratic-radiation blocks, internal-graviton loops or full independent curved matching.",
    }
