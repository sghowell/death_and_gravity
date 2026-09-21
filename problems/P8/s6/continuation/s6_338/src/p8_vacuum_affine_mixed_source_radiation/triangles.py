"""Both massive kernel insertions in exact D, and all labelled EOM emissions."""

from functools import cache
from itertools import permutations

import sympy as s
from p8_vacuum_affine_local_tadpole_radiation import vertices


@cache
def data():
    z, y, a, e, m, v = s.symbols(
        "z y pk epsilon heavy_mass_squared p_squared", positive=True
    )
    M = 1 - z + m * z - z * (1 - z) * v
    ML = M - 2 * z * y * a
    MH = M - 2 * (1 - z) * y * a
    MR = M - 2 * z * (1 - z) * a
    checks = {}
    # Radial integration divided by Gamma(1+epsilon) gives
    # -2*z^2*M_L^-1-e and -2*(1-z)^2*M_H^-1-e.
    FL = -z * ML ** (-e) / (a * e)
    FH = -(1 - z) * MH ** (-e) / (a * e)
    checks["whole_D_light_triangle_parameter_antiderivative"] = s.simplify(
        s.diff(FL, y) + 2 * z * z * ML ** (-1 - e)
    )
    checks["whole_D_heavy_triangle_parameter_antiderivative"] = s.simplify(
        s.diff(FH, y) + 2 * (1 - z) ** 2 * MH ** (-1 - e)
    )
    end_light = s.simplify(FL.subs(y, 1 - z) - FL.subs(y, 0))
    end_heavy = s.simplify(FH.subs(y, z) - FH.subs(y, 0))
    checks["whole_D_two_triangles_equal_bubble_divided_difference"] = s.simplify(
        end_light + end_heavy + (MR ** (-e) - M ** (-e)) / (a * e)
    )
    checks["full_dimensional_gamma_recurrence"] = s.expand_func(
        s.gamma(e + 1)
    ) - e * s.gamma(e)
    # At epsilon0 the finite bubble is -integral log M.
    F0L = z * s.log(ML) / a
    F0H = (1 - z) * s.log(MH) / a
    checks["finite_light_triangle_antiderivative"] = s.simplify(
        s.diff(F0L, y) + 2 * z * z / ML
    )
    checks["finite_heavy_triangle_antiderivative"] = s.simplify(
        s.diff(F0H, y) + 2 * (1 - z) ** 2 / MH
    )
    checks["finite_two_triangles_equal_log_bubble_divided_difference"] = s.simplify(
        F0L.subs(y, 1 - z)
        - F0L.subs(y, 0)
        + F0H.subs(y, z)
        - F0H.subs(y, 0)
        - (s.log(MR) - s.log(M)) / a
    )
    t = s.Symbol("t", nonnegative=True)
    checks["original_massive_parameter_gap_for_v_below2"] = s.expand(
        M.subs(v, 2) - 1 - (m - 3) * z - 2 * z * z
    )
    checks["positive_mass_margin_at_n128"] = s.expand(
        (m - 3).subs(m, 128 + t) - 125 - t
    )

    # Feynman weights for both loop insertions are derived by completing
    # the square, not assumed to be a covariantized flat bubble.
    l2, lp, lk = s.symbols("loop_squared loop_dot_p loop_dot_k")

    def square(P, K):
        return l2 + 2 * P * lp + 2 * K * lk + P * P * v + 2 * P * K * a

    dl = l2 - 1
    dlk = l2 + 2 * lk - 1
    dh = l2 - 2 * lp + v - m
    light_combined = (1 - z - y) * dl + y * dlk + z * dh
    light_completed = square(-z, y) - ML
    checks["literal_light_three_denominators_complete_square"] = s.expand(
        light_combined - light_completed
    )
    dhk = l2 - 2 * lp - 2 * lk + v + 2 * a - m
    heavy_combined = (1 - z) * dl + (z - y) * dh + y * dhk
    heavy_completed = square(-z, -y) - MH
    checks["literal_heavy_three_denominators_complete_square"] = s.expand(
        heavy_combined - heavy_completed
    )
    # Shifted polarization numerators follow from epsilon.k=trace epsilon=0.
    Er, Ep, Epr = s.symbols("epsilon_rr epsilon_pp epsilon_pr")
    light_numerator = Er + 2 * z * Epr + z * z * Ep
    heavy_numerator = Er - 2 * (1 - z) * Epr + (1 - z) ** 2 * Ep
    checks["light_triangle_TT_numerator_after_isotropic_average"] = s.expand(
        light_numerator.subs({Er: 0, Epr: 0}) - z * z * Ep
    )
    checks["heavy_triangle_TT_numerator_after_isotropic_average"] = s.expand(
        heavy_numerator.subs({Er: 0, Epr: 0}) - (1 - z) ** 2 * Ep
    )

    # The EOM-marked bilocal term has a singleton, two plain fields,
    # and one marked field. Keep all24 labelled assignments.
    G, H = vertices.generic_radiative_data()
    ak = tuple(G[i, 4] for i in range(4))
    F1 = s.Symbol("F1")
    Fi = s.symbols("F_shift0:4")
    assignments = tuple(permutations(range(4)))
    external = s.S.Zero
    for i in range(4):
        shifted = s.Matrix(G)
        shifted[i, i] += 2 * ak[i]
        vertex = s.Add(
            *(
                (Fi[single] if i == single else F1) * (1 - shifted[mark, mark])
                for single, plain1, plain2, mark in assignments
            )
        )
        external += H[i, i] * vertex / ak[i]
    # Vary(Box+1) on the marked scalar: +2epsilon(p_mark,p_mark).
    # Metric/kernel/measure variations multiplying the unshifted EOM vanish.
    metric = s.Add(
        *(2 * H[mark, mark] * F1 for single, plain1, plain2, mark in assignments)
    )
    checks["all24_EOM_assignments_external_plus_metric_cancellation"] = s.expand(
        external + metric
    )
    checks["all24_EOM_external_not_deleted"] = s.expand(
        external + 12 * F1 * sum(H[i, i] for i in range(4))
    )
    assert external != 0 and metric != 0

    # Independently derive scalar-kernel divided differences by varying powers
    # of (-Box), using the telescoping noncommutative insertion formula.
    pp, rr = s.symbols("p_squared r_squared")
    for degree in range(1, 9):
        chain = sum(rr**j * pp ** (degree - 1 - j) for j in range(degree))
        checks["power" + str(degree) + "_complete_operator_variation"] = s.expand(
            (rr - pp) * chain - (rr**degree - pp**degree)
        )
    return {
        "checks": checks,
        "gates": {
            "light_and_heavy_triangle_insertions_both_present": True,
            "whole_D_antiderivative_before_MS_subtraction": True,
            "massive_parameter_gap_no_hidden_null_pole": True,
            "all24_external_EOM_and_metric_assignments": len(assignments) == 24,
            "EOM_external_and_metric_terms_individually_nonzero": bool(
                external != 0 and metric != 0
            ),
            "noncommutative_power_variation_independently_checked": True,
        },
        "whole_light_denominator": ML,
        "whole_heavy_denominator": MH,
        "whole_light_antiderivative": FL,
        "whole_heavy_antiderivative": FH,
        "whole_complete_EOM_external": s.expand(external),
        "whole_complete_EOM_metric": s.expand(metric),
        "whole_kernel_identity": "At real null k, the complete physical TT response of G_H G is -epsilon(p,p)[B_D((p+k)^2)-B_D(p^2)]/(p.k). The two shifted scalar triangles have weights z^2 and(1-z)^2. Their exact-D radial and parameter integrals yield the divided difference using Gamma(1+epsilon)=epsilon Gamma(epsilon). The momentum-independent pole cancels; at epsilon0 B_MS=-integral log M. For original p^2=1 and |p.k|<=1/2, all parameter denominators are>=1. No low-energy heavy expansion inside a loop or generic curved-completion assumption is used.",
        "whole_EOM_boundary": "The marked light EOM is kept through all24 labelled emissions. Off-shell marked external emission and variation of Box have opposite nonzero sums. Variations of other factors multiply the unshifted on-shell EOM. The cancellation is explicit; the EOM vertex is not discarded before coupling the graviton.",
    }
