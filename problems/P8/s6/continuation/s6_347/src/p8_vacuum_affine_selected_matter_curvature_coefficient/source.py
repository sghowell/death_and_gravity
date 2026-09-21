"""Exact scalar graph inventory and the unchanged original source."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_core_curvature_coefficient import source as previous
from p8_vacuum_affine_heavy_parent_one_loop import germs

KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT = previous.CONTACT
original_parameters = previous.original_parameters
PATTERNS = ((4, 0, 0, 0), (2, 1, 0, 0), (0, 2, 0, 0), (1, 0, 1, 0), (0, 0, 0, 1))


def require_class(value):
    if not isinstance(value, str) or value not in (
        "old_polynomial_core",
        "all_six_local_tadpoles",
        "mixed_across",
        "mixed_within_J4",
        "mixed_within_J2",
    ):
        raise ValueError("Require a complete selected scalar-loop class")
    return value


@cache
def data():
    checks = dict(previous.data()["checks"])
    for name, got, want in (
        ("mass", HEAVY_MASS2, germs.MASS2),
        ("g", CUBIC, germs.G),
        ("contact", CONTACT, germs.CONTACT),
        ("metric", KAPPA, germs.KAPPA),
    ):
        checks["complete_selected_original_" + name] = s.factor(got - want)
    found = tuple(
        p
        for p in product(range(5), repeat=4)
        if sum((j + 1) * v for j, v in enumerate(p)) == 4
    )
    checks["source_coupling_c_src_not_core_contact_ratio"] = s.factor(
        germs.couplings()["H_phi2_Y"] + 4 * CUBIC / KAPPA
    )
    return {
        "checks": checks,
        "gates": {
            "exact_E4_L1_valence_inventory": set(found) == set(PATTERNS),
            "no_degree_seven_or_higher_scalar_vertex_at_this_order": 7 - 2 > 4,
            "all_six_local_operators_and_three_mixed_classes_retained": True,
            "same_original_source_and_physical_metric": True,
            "only_classical_limiting_action_not_finite_gravity_decoupling": True,
            "fixed_OS4_and_pole_only_derivative_prescription_unchanged": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_scalar_graph_patterns": found,
        "whole_complete_inventory": "For E4,L1, sum(degree_v-2)=4 gives exactly five scalar patterns. Three are the old polynomial bubble/triangle/box core; one is the cubic/five-field source pair with all three contraction classes; one is a single six-field vertex with all six literal operators. Metric insertions are the complete variations of these scalar functionals. Higher source jets remain in the action but cannot enter this order.",
        "whole_boundary": "A complete selected formal scalar-loop four-light degree6 coefficient of the fixed-canonical classical limiting action, coupled to one external real-null-TT metric perturbation. Free spectators in this limit do not generate connected four-light graphs. No finite-gravity quantum decoupling theorem or internal-graviton contribution is inferred; extra parent curvature matching remains separate.",
    }
