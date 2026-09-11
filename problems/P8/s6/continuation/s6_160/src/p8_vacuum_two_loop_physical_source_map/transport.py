"""Physical-source equality without assuming a small ordinary-Psi residue."""

from functools import cache

import sympy as s


@cache
def data():
    h = s.Symbol("h")
    a0, a1, a2, o1, o2 = s.symbols("A0 A1 A2 o1 o2")
    overlap = 1 + h * o1 + h * h * o2
    amplitude = a0 + h * a1 + h * h * a2
    residue = s.series(overlap**2, h, 0, 3).removeO()
    numerator = s.series(overlap**4 * amplitude, h, 0, 3).removeO()
    amputated = s.series(numerator / residue**4, h, 0, 3).removeO()
    physical = s.series(residue**2 * amputated, h, 0, 3).removeO()
    K, p1, f0, f1, e = s.symbols("K pole1 finite0 epsilon1 epsilon")
    raw = p1 / e + f0 + e * f1
    factor = 1 + h * K + e * h * s.Symbol("Kepsilon")
    return {
        "physical_source_prescription": "Define the transformed regulated generating functional by the exact pullback of the canonical renormalized parent: full S[F(Psi),H,fermions,gauge], transformed fixed parent counterterms, det(F') and J F(Psi)+J_H H. Correlators of this physical source, not an independently minimally subtracted ordinary Psi, are the observables.",
        "continuum_Ward_identity": "The action variation, Jacobian divergence and physical-source variation cancel by the regulated total functional derivative. Iterate through map order three for four physical light sources at two loops. Every local differential Jacobian ghost loop is scaleless in dimensional regularization before finite parts. Heavy elimination commutes because F is independent of H.",
        "formal_ordinary_overlap_for_diagnostic_only": overlap,
        "ordinary_residue_diagnostic": residue,
        "ordinary_amputated_vertex_diagnostic": amputated,
        "proper_LSZ_diagnostic": physical,
        "checks": {
            "ordinary_residue_square_through_two_loops": s.expand(
                residue - 1 - 2 * h * o1 - h * h * (o1 * o1 + 2 * o2)
            ),
            "ordinary_four_vertex_first_coefficient": s.expand(
                amputated.coeff(h, 1) - a1 + 4 * o1 * a0
            ),
            "ordinary_four_vertex_second_coefficient": s.expand(
                amputated.coeff(h, 2) - a2 + 4 * o1 * a1 - (10 * o1**2 - 4 * o2) * a0
            ),
            "all_formal_overlap_factors_cancel": s.expand(physical - amplitude),
            "finite_regulator_pairing_keeps_source_products": s.expand(
                factor * (raw - p1 / e)
            ).coeff(e, 0)
            - (1 + h * K) * f0,
        },
        "not_claimed": "No numerical bound or positivity for o2 or the ordinary Psi two-loop residue. Composite-operator proper subgraphs may mix with fermion operators; their separate MS prescription is not replaced by the physical-source pullback. No global off-shell inverse, bare derivative-truncation equivalence or physical higher-order error is inferred.",
    }
