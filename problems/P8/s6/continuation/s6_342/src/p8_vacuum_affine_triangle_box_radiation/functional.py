"""Full labelled noncommuting covariant trace variation and actual finite controls."""

import itertools
from functools import cache

import sympy as s


@cache
def data():
    checks = {}
    gates = {}
    h = s.Symbol("metric_h")
    L, H, S, P = s.symbols("L H S P", commutative=False)
    dL, dH, dS = s.symbols("dL dH dS", commutative=False)
    triangle = (L + h * dL) * (S + h * dS) * (L + h * dL) * P * (H + h * dH) * P
    triangle_parts = (
        dL * S * L * P * H * P,
        L * dS * L * P * H * P,
        L * S * dL * P * H * P,
        L * S * L * P * dH * P,
    )
    box = (L + h * dL) * P * (H + h * dH) * P * (L + h * dL) * P * (H + h * dH) * P
    box_parts = (
        dL * P * H * P * L * P * H * P,
        L * P * dH * P * L * P * H * P,
        L * P * H * P * dL * P * H * P,
        L * P * H * P * L * P * dH * P,
    )
    checks["generic_triangle_all_noncommuting_metric_positions"] = s.expand(
        s.diff(triangle, h).subs(h, 0) - sum(triangle_parts)
    )
    checks["generic_box_all_noncommuting_metric_positions"] = s.expand(
        s.diff(box, h).subs(h, 0) - sum(box_parts)
    )
    KK, GG, VV = s.symbols("K_inverse G_inverse delta_K_inverse", commutative=False)
    inverse_relation = KK * (-GG * VV * GG) + VV * GG
    checks["inverse_operator_variation_sign"] = (
        s.expand(inverse_relation).subs(KK * GG, s.S.One).expand()
    )
    gates["inverse_operator_wrong_sign_does_not_cancel"] = (
        s.expand(KK * (GG * VV * GG) + VV * GG).subs(KK * GG, s.S.One).expand() != 0
    )
    gates["unreduced_inverse_identity_is_nontrivial"] = s.expand(inverse_relation) != 0
    # Independent finite matrices: direct four-field + metric coefficient
    # versus the entire labelled trace variation. These are diagnostic
    # noncommuting kernels, not a replacement mass or physical propagator.
    G2, C = s.symbols("g_squared contact")
    GL = s.Matrix([[3, 1, 0], [1, 2, 1], [0, 1, 4]])
    GH = s.Matrix([[5, 0, 1], [0, 3, 1], [1, 1, 2]])
    VL = s.Matrix([[1, 2, 0], [2, -1, 1], [0, 1, 3]])
    VH = s.Matrix([[2, -1, 1], [-1, 1, 0], [1, 0, -2]])
    dGL = -GL * VL * GL
    dGH = -GH * VH * GH
    gates["diagnostic_light_and_heavy_noncommute"] = GL * GH != GH * GL
    gates["literal_inverse_derivative_both_sides"] = (
        GL.inv() * dGL + VL * GL == s.zeros(3)
        and dGH * GH.inv() + GH * VH == s.zeros(3)
    )
    fields = (
        s.Matrix([1, 2, 0]),
        s.Matrix([0, 1, 3]),
        s.Matrix([2, -1, 1]),
        s.Matrix([1, 0, -2]),
    )
    ts = s.symbols("field0:4")
    phi = sum((t * f for t, f in zip(ts, fields)), s.zeros(3, 1))
    M = s.diag(*phi)
    LM, HM = GL + h * dGL, GH + h * dGH
    K = C * s.eye(3) + G2 * GH
    dK = G2 * dGH
    KM = C * s.eye(3) + G2 * HM
    mpair = lambda f, g: s.diag(*(K * f.multiply_elementwise(g)))
    dmpair = lambda f, g: s.diag(*(dK * f.multiply_elementwise(g)))
    SM = s.diag(*(KM * phi.multiply_elementwise(phi)))
    directT = -G2 * s.trace(LM * SM * LM * M * HM * M) / 4
    directD = -(G2**2) * s.trace(LM * M * HM * M * LM * M * HM * M) / 4

    def all_fields_coefficient(expr, metric_order):
        return s.Poly(s.expand(expr), *ts, h).coeff_monomial((*([1] * 4), metric_order))

    labelT = s.S.Zero
    labelD = s.S.Zero
    varyT = [s.S.Zero] * 4
    varyD = [s.S.Zero] * 4
    for a, b, c, d in itertools.permutations(range(4)):
        S0 = mpair(fields[a], fields[b])
        dS0 = dmpair(fields[a], fields[b])
        Mc = s.diag(*fields[c])
        Md = s.diag(*fields[d])
        labelT += -G2 * s.trace(GL * S0 * GL * Mc * GH * Md) / 4
        Tterms = (
            dGL * S0 * GL * Mc * GH * Md,
            GL * dS0 * GL * Mc * GH * Md,
            GL * S0 * dGL * Mc * GH * Md,
            GL * S0 * GL * Mc * dGH * Md,
        )
        for j, A in enumerate(Tterms):
            varyT[j] += -G2 * s.trace(A) / 4
        Ms = [s.diag(*fields[j]) for j in (a, b, c, d)]
        labelD += (
            -(G2**2) * s.trace(GL * Ms[0] * GH * Ms[1] * GL * Ms[2] * GH * Ms[3]) / 4
        )
        Ds = (
            dGL * Ms[0] * GH * Ms[1] * GL * Ms[2] * GH * Ms[3],
            GL * Ms[0] * dGH * Ms[1] * GL * Ms[2] * GH * Ms[3],
            GL * Ms[0] * GH * Ms[1] * dGL * Ms[2] * GH * Ms[3],
            GL * Ms[0] * GH * Ms[1] * GL * Ms[2] * dGH * Ms[3],
        )
        for j, A in enumerate(Ds):
            varyD[j] += -(G2**2) * s.trace(A) / 4
    checks["all24_labelled_triangles_from_actual_functional"] = s.expand(
        all_fields_coefficient(directT, 0) - labelT
    )
    checks["all24_labelled_boxes_from_actual_functional"] = s.expand(
        all_fields_coefficient(directD, 0) - labelD
    )
    checks["triangle_all_three_lines_and_outer_K_from_metric_derivative"] = s.expand(
        all_fields_coefficient(directT, 1) - sum(varyT)
    )
    checks["box_all_four_lines_from_metric_derivative"] = s.expand(
        all_fields_coefficient(directD, 1) - sum(varyD)
    )
    gates["each_required_triangle_graph_class_nonzero"] = all(
        s.expand(v) != 0 for v in varyT
    )
    gates["each_required_box_line_class_nonzero"] = all(s.expand(v) != 0 for v in varyD)
    gates["discarding_outer_K_variation_fails"] = (
        s.expand(all_fields_coefficient(directT, 1) - sum(varyT[j] for j in (0, 2, 3)))
        != 0
    )

    return {
        "checks": {key: s.expand(value) for key, value in checks.items()},
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_triangle_metric_positions": triangle_parts,
        "whole_box_metric_positions": box_parts,
        "whole_diagnostic_matrices": {
            "light": GL,
            "heavy": GH,
            "light_vertex": VL,
            "heavy_vertex": VH,
        },
        "whole_diagnostic_four_field_derivatives": {
            "triangle": s.expand(labelT),
            "box": s.expand(labelD),
            "triangle_metric": s.expand(sum(varyT)),
            "box_metric": s.expand(sum(varyD)),
        },
        "whole_covariant_functional": "Wloc=-M_(K Phi^2)/2,Wbi=-g^2 M_Phi G_H M_Phi,K=C I+g^2 G_H. From-Tr(G_L W G_L W)/4, Gamma_T=-g^2/4 Tr[G_L M_(K Phi^2) G_L M_Phi G_H M_Phi] and Gamma_D=-g^4/4 Tr[G_L M_Phi G_H M_Phi G_L M_Phi G_H M_Phi]. Differentiate all24 labelled fields and every inverse, including K, before projection.",
        "whole_metric_terms": "The triangle has both light loop lines, the heavy loop line and the outer K heavy propagator. The box has two light and two heavy loop lines. Potential vertices and volume variations are retained until their pure trace vanishes in real TT projection. No vertex is discarded using its flat on-shell value.",
        "whole_diagnostic_boundary": "The finite three-site matrices are independent noncommuting controls of the literal trace derivatives, not physical propagators or replacement masses. Generic noncommuting Leibniz and inverse identities supply the operator argument.",
    }
