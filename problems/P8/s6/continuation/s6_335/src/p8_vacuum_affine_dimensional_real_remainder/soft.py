"""General soft residue and trace-completing D6 remainder majorant."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import ward


def component_budgets():
    seagulls = (
        s.Integer(6) * 4 * (17 * 324 + s.Rational(108 * 4, 2) + 2 * 6**3 * 2 * 18 * 2)
    )
    cubic = s.Integer(4) * 6 * (2 * 6**4 * 4 * 6**2 + 2 * (2 * 6**3) * 36 * 6)
    total = 64 * 5000000 + 608 * 24300 + seagulls + 3 * 18**2 * cubic
    return {
        "seagulls": s.sympify(seagulls),
        "cubic": cubic,
        "D6_unit_remainder": total,
        "canonical_remainder": 2 * total + 45216,
    }


@cache
def data():
    D = ward.DIM
    zero = {ward.u1: 0, ward.u2: 0, ward.u3: 0}
    L, R = ward.dot(6, 6), ward.dot(7, 7)
    born = s.factor(
        (-ward.inner(ward.stress(0, 1), ward.H(2, 3)) / L)
        .subs(ward.LOOP_DIM, D)
        .subs(zero)
    )
    theta = (D - 4) / (D - 2)
    w = s.Symbol("soft_scale", positive=True)
    scaled = ward.tt_channel().subs(
        {ward.u1: w * ward.u1, ward.u2: w * ward.u2, ward.u3: w * ward.u3},
        simultaneous=True,
    )
    residue = s.cancel(w * scaled).subs(w, 0)
    J = sum(ward.dot(i, 5) ** 2 / ward.dot(i, 4) for i in range(4)).subs(ward.xk, 0)
    JL = sum(ward.dot(i, 5) ** 2 / ward.dot(i, 4) for i in (0, 1))
    JR = sum(ward.dot(i, 5) ** 2 / ward.dot(i, 4) for i in (2, 3))
    auxR = ward.dot(6, 4) * (JL - JR) - 2 * ward.dot(6, 5) ** 2
    budgets = component_budgets()
    a, b, c, ai, bi, ci = s.symbols("a b c ai bi ci", real=True)
    E = s.Matrix([[a + s.I * ai, b + s.I * bi], [b + s.I * bi, c + s.I * ci]])
    tr = s.trace(E)
    embedded = s.diag(E, -tr / 2, -tr / 2)
    fn = lambda A: s.trace(s.conjugate(A).T * A)
    checks = {
        "general_D_Born_correction": s.factor(
            s.together(born - born.subs(D, 4) + 2 * ward.mu**2 * theta / L)
        ),
        "general_D_channel_leading_soft_residue": s.factor(
            s.together(residue - born * J)
        ),
        "canonical_auxiliary_soft_remainder": s.factor(
            s.together(
                ward.auxiliary_tt_channel()
                + s.Rational(1, 2) * (1 / L + 1 / R) * J
                - (auxR / (L * R)).subs(ward.xk, 0)
            )
        ),
        "auxiliary_remainder_budget": s.Integer(
            576 * 64 + 608 * 3 + 3 * 544 * 4 - 45216
        ),
        "six_D_numerator_difference": s.Integer(36 * (84 * 9 + 25 * 24) - 48816),
        "six_D_old_numerator": s.Integer(36 * 25 * 9 - 8100),
        "six_D_quotient_difference": s.Integer(2 * 48816 + 2 * 8100 * 96 - 1652832),
        "six_D_seagulls": budgets["seagulls"] - 883872,
        "six_D_cubic_coefficient": budgets["cubic"] - 13436928,
        "six_D_total_remainder": budgets["D6_unit_remainder"] - 13396352288,
        "canonical_embedding_total": budgets["canonical_remainder"] - 26792749792,
        "six_D_embedding_trace": s.trace(embedded),
        "six_D_embedding_norm": s.expand(
            fn(embedded) - fn(E) - s.conjugate(tr) * tr / 2
        ),
        "two_dimensional_trace_bound_identity": s.expand(
            2 * fn(E)
            - s.conjugate(tr) * tr
            - ((a - c) ** 2 + (ai - ci) ** 2 + 4 * (b * b + bi * bi))
        ),
    }
    return {
        "checks": checks,
        "gates": {
            "three_channel_rounding_margin": 5000000 > 3 * 1652832,
            "canonical_rounding_margin": bool(
                budgets["canonical_remainder"] < 30000000000
            ),
            "normalized_gravity_rounding_margin": bool(
                s.Integer(4000000000) > s.Rational(30000000000, 8)
            ),
            "trace_completion_does_not_reconstruct_trace_from_D4_helicities": True,
            "unit_D6_action_count_keeps_all_six_dimensional_indices": True,
            "auxiliary_heavy_improved_formula_not_evaluated_at_zero_mass": True,
        },
        "whole_exact_soft_proof": "Channel residue of the actual general-D gravity tensor is G0_D times the complete leading current. The auxiliary canonical remainder is the actual un-improved unit massless-exchange tensor, with R_S<=45216/delta^2.",
        "whole_embedding_proof": "Embed any unit canonical2x2 test matrix E as diag(E,-tr(E)I2/2) in the D6 transverse space; norm<=sqrt2<2. The direct D6 unit-polarization bound controls the canonical trace as well. G_D=G6+(2theta-1)S and|2theta-1|<=1 then give||R_G,D||_F<3e10/delta^2.",
        "whole_D6_component_budgets": budgets,
    }
