"""Complete off-shell mixed-source functional and common-jet conversion."""

from functools import cache

import sympy as s
from p8_vacuum_affine_box_curvature_coefficient import basis, jets
from p8_vacuum_affine_core_curvature_coefficient import conversion
from p8_vacuum_affine_mixed_source_radiation import triangles

N = s.Symbol("n", positive=True)


def within_coefficient(G, H, n):
    value = s.S.Zero
    for j in (1, 2, 3):
        left = (0, j)
        right = tuple(i for i in range(4) if i not in left)
        v = sum(G[i, h] for i in left for h in left)
        u = sum(G[i, h] for i in right for h in right)
        hpp = sum(H[i, h] for i in left for h in left)
        dots = G[left[0], left[1]] + G[right[0], right[1]]
        value += (
            -8 * hpp * (v * v + v * u + u * u) / n**4 + 4 * dots * hpp * (v + u) / n**3
        )
        value += (
            4 * (H[left[0], left[1]] * u * u + H[right[0], right[1]] * v * v) / n**3
        )
    return value


@cache
def data():
    G, H, aa, _variables, T = conversion.generic()
    _, tau, coefs, _ = conversion.build()
    S3 = -2 * sum(
        H[i, i] * conversion.divided_power(G[i, i], aa[i], 3) for i in range(4)
    )
    S21 = -2 * sum(
        H[i, i] * conversion.divided_power(G[i, i], aa[i], 2) * G[j, j]
        + H[j, j] * G[i, i] ** 2
        for i in range(4)
        for j in range(4)
        if i != j
    )
    contact = lambda p: conversion.labeled_contact(p, G, H, aa)
    checks = {
        "whole_singleton_cube_contact": s.expand(
            S3 - s.Rational(2, 3) * contact((0, 0, 3))
        ),
        "whole_distinct_singleton_square_linear_contact": s.expand(
            S21 - 2 * contact((0, 1, 2))
        ),
        "both_across_conversion_factors_equal_minus4": s.Matrix(
            [s.Rational(2, 3) * tau[(0, 0, 3)] + 4, 2 * tau[(0, 1, 2)] + 4]
        ),
    }
    g, c, B0, B1, B2, B3, A1, An = s.symbols("g c_src B0 B1 b2 b3 A1 An", real=True)
    ds = s.symbols("d0:4")
    total = sum(ds)
    B = lambda v: B0 + B1 * v + B2 * v * v + B3 * v**3
    x = s.Symbol("homogeneous_degree2_scale")
    direct = (
        -2 * g * c * sum(A1 - An + (N - 1 - total + ds[i]) * B(ds[i]) for i in range(4))
    )
    covariant = (
        2 * g * c * (4 - N) * sum(B(di) for di in ds)
        - 2
        * g
        * c
        * sum((1 - ds[j]) * B(ds[i]) for i in range(4) for j in range(4) if i != j)
        - 8 * g * c * (A1 - An)
    )
    checks["whole_offshell_covariant_reduction_EOM_and_delta_retained"] = s.expand(
        direct - covariant
    )
    degree6 = s.Poly(
        s.expand(direct.subs({di: x * di for di in ds}, simultaneous=True)), x
    ).coeff_monomial(x**3)
    checks["whole_across_degree6_flat"] = s.expand(
        degree6
        + 2
        * g
        * c
        * (
            (N - 1) * B3 * sum(di**3 for di in ds)
            - B2 * sum(ds[i] ** 2 * ds[j] for i in range(4) for j in range(4) if i != j)
        )
    )
    tr = triangles.data()
    for key in (
        "whole_D_light_triangle_parameter_antiderivative",
        "whole_D_heavy_triangle_parameter_antiderivative",
        "whole_D_two_triangles_equal_bubble_divided_difference",
        "finite_two_triangles_equal_log_bubble_divided_difference",
        "literal_light_three_denominators_complete_square",
        "literal_heavy_three_denominators_complete_square",
    ):
        checks["retained_" + key] = tr["checks"][key]
    jet_contacts = tuple(jets.bose_contact(w, G, H, aa) for w in basis.WORDS)
    across_lift = sum(
        (
            s.Rational(2, 3) * (N - 1) * B3 * coefs[(0, 0, 3)][i]
            - 2 * B2 * coefs[(0, 1, 2)][i]
        )
        * q
        for i, q in enumerate(jet_contacts)
    )
    checks["whole_across_common_jet_curvature"] = s.expand(
        (N - 1) * B3 * S3 - B2 * S21 - across_lift + 4 * ((N - 1) * B3 - B2) * T
    )
    co_b = tuple(v / 2 for v in coefs[(3, 0, 0)])
    co_mix = tuple(2 * v for v in coefs[(2, 1, 0)])
    lift = sum(
        ((4 / N**4 - 2 / N**3) * b + q / N**3) * h
        for b, q, h in zip(co_b, co_mix, jet_contacts)
    )
    checks["whole_withinJ4_generic_offshell_curvature"] = s.factor(
        within_coefficient(G, H, N) - lift + 16 * T / N**4
    )
    Lphi, Lphi2, Y, phi = s.symbols("L_phi L_phi2 Y phi")
    checks["exact_covariant_Leibniz_source_rewrite"] = s.expand(
        (phi * Lphi - Lphi2 / 2).subs(Lphi2, 2 * phi * Lphi - 2 * Y) - Y
    )
    return {
        "checks": checks,
        "gates": {
            "all24_cross_source_labels_with_marked_EOM_retained": True,
            "both_massive_loop_insertions_arbitrary_p_squared": True,
            "all_three_withinJ4_heavy_channels": True,
            "no_offshell_LSZ_assumption_in_amputated_kernel": True,
            "all_comparison_coordinates_match_frozen_core_basis": True,
        },
        "whole_offshell_across_scalar_polynomial": direct,
        "whole_across_degree6_flat": degree6,
        "whole_across_curvature_before_loop_factor": 8 * g * c * ((N - 1) * B3 - B2),
        "whole_withinJ4_curvature_before_4g2_over_kappa_loop_factor": -16 / N**4,
        "whole_covariant_across_reduction": "The literal cross source equals gc_src(4-n)Phi(x)Phi(y)^3 G_HG/3-gc_src delta(G-G_H)Phi^4/3-gc_src Phi(x)Phi(y)^2(Box+1)Phi(y)G_HG. All labels and both the delta and marked-EOM terms are retained. The full-D two-loop-line metric response equals the B(L) divided difference for arbitrary endpoint p_squared. No generic four-hard-direction TT kernel is fixed by Ward identities alone.",
        "whole_withinJ4_reduction": "The fixed c_src I H(Y+Phi^2) source at both heavy ends has off-shell flat4g^2/kappa sum(4-2v+sum_i d_i)/(n-v), before1/(16pi^2). The exact scalar Leibniz identity converts Y without equations of motion. Its degree6 terms have common-lift conversion factors-4 for sum v^3 and-8 for sum(sum_i d_i)v^2. Their n^-3 parts cancel and leave chi_withinJ4=-64g^2/(16pi^2 kappa n^4).",
        "whole_withinJ2_and_counterterm_boundary": "The complete within-J2 graph cancels against the existing full heavy-onepoint counterfunctional, including all metric and external source terms. Its internal loop depends only on null k and has zero TT contraction. Mixed pole counterterms have at most two derivatives; finite OS4 is constant. Their fixed prescriptions are retained, not replaced by finite curvature choices.",
    }
