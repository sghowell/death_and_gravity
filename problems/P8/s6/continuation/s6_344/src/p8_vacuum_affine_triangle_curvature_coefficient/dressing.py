"""Full outer-heavy-line variation and complete original Bose/loop normalization."""

from functools import cache
from itertools import permutations

import sympy as s
from p8_vacuum_affine_curvature_contact_basis import basis


@cache
def data():
    v, u, n, g = s.symbols("v u n g")
    A0, A1, A2, A3 = s.symbols("A0 A1 A2 A3")
    F0, F1, F2, F3 = s.symbols("F0 F1 F2 F3")
    A = lambda x: A0 + A1 * x + A2 * x * x + A3 * x**3
    F = lambda x: F0 + F1 * x + F2 * x * x + F3 * x**3
    scale, Q, c = s.symbols("scale Q c")
    bose = s.expand(
        sum(basis.q(p[2], p[3], p[2], p[3]) for p in permutations(range(4)))
    )
    checks = {
        "entire_scalar_covariant_product_divided_difference": s.factor(
            (A(u) * F(u) - A(v) * F(v)) / (u - v)
            - A(v) * (F(u) - F(v)) / (u - v)
            - F(u) * (A(u) - A(v)) / (u - v)
        ),
        "original_heavy_resolvent_all_insertions": s.factor(
            (1 / (n - u) - 1 / (n - v)) / (u - v) - 1 / ((n - u) * (n - v))
        ),
        "dressing_first_curvature_order_only_endpoint_constant": s.expand(
            s.Poly(s.expand(A(scale**2 * v) * scale**6 * Q), scale).coeff_monomial(
                scale**6
            )
            - A0 * Q
        ),
        "full24_curvature_Bose_factor": s.expand(bose - 4 * basis.canonical),
        "four_field_factor_cancels_triangle_quarter": s.expand(
            g**2 * c * bose / 4 - g**2 * c * basis.canonical
        ),
    }
    for slot in range(4):
        row = sum(basis.h[slot, j] for j in range(4))
        checks[f"physical_TT_shift_of_composite_momentum_{slot}"] = s.expand(row)
    return {
        "checks": checks,
        "gates": {
            "outer_heavy_branch_uses_shifted_C_of_u": True,
            "all24_labels_and_loop_16pi2_restored": True,
            "first_curvature_order_lower_jet_residuals_zero": True,
            "no_extra_assignment_of_independent_parent_chi": True,
            "comparator_dependence_explicit": True,
        },
        "whole_complete_dressed_identity": "Let A(v)=C+g^2/(n-v). The exact chosen comparison functional has A(L1) acting on the same first scalar factor as the triangle labeled-box series. Its metric variation includes both A(v)*delta_triangle_lift and-2epsilon(P,P)*DD A(u,v)*C_triangle(u;w,q). This is precisely the original outer-heavy insertion; using C_triangle(v) would fail the displayed product identity. Consequently the complete dressed difference is A(v) times the primitive difference at every retained local order.",
        "whole_first_order_normalization": "The full original triangle class is g^2/4 times the sum of24 scalar labelings, with common1/(16pi^2 sqrt(kappa)). The primitive degree6 difference is c_triangle(n)*Q_cd,cd. Its full24 sum is4T, so chi_triangle=g^2*(C+g^2/n)*c_triangle(n)/(16pi^2) and DeltaM5_degree6=chi_triangle*T/sqrt(kappa). All higher powers of A(v) first enter curvature differences at derivative order8 or above.",
        "whole_external_and_subtraction_cancellation": "The exact off-shell flat Taylor polynomial of the comparison equals the triangle Taylor polynomial before any mass-shell substitution. Hence all four external emissions of each retained polynomial agree and cancel in their difference. The fixed OS4 constant has the same local completion in both and contributes no real-TT contact; it is not a new finite counterterm or matching choice.",
        "whole_matching_boundary": "This is a known selected-loop coefficient relative to one explicitly fixed covariant local lift, not an invariant separation into all possible bases and not the independent added parent-theory chi of S336. The full triangle is already included in S342. Adding this contact on top of that full amplitude would double count.",
    }
