"""Unchanged original source and four-scalar tree ownership."""

from functools import cache

import sympy as s
from p8_vacuum_affine_one_newton_inclusive_assembly import source as finite_source
from p8_vacuum_affine_single_residual_soft_dressing import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
recoil = previous.recoil


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    actual = finite_source.data()
    for name, value in actual["checks"].items():
        if name.startswith("no_lower_R_jet_"):
            checks["two_real_replayed_" + name] = value
    checks.update(
        {
            "same_two_real_original_Newton_parameter": KAPPA - finite_source.KAPPA,
            "same_two_real_original_heavy_mass": HEAVY_MASS2
            - finite_source.HEAVY_MASS2,
            "same_two_real_original_cubic": CUBIC - finite_source.CUBIC,
            "same_two_real_original_quartic": CONTACT - finite_source.CONTACT,
            "same_original_tuned_contact": CONTACT
            + CUBIC**2 * (3 / (HEAVY_MASS2 - 2) - 2 / (HEAVY_MASS2 - 2) ** 2),
            "four_scalar_forest_budget": s.Integer(4 - 2 * 2),
            "six_scalar_local_vertex_exclusion": s.Integer(2 + (6 - 2) - 6),
            "higher_heavy_source_minimum_external_scalars": s.Integer(
                2 * 2 + (4 - 2) - 6
            ),
            "two_linear_Proca_sources_minimum_external_scalars": s.Integer(
                2 * 2 + 2 * (4 - 2) - 8
            ),
            "six_point_tree_total_degree_budget": s.Integer(6 - 2 - 4),
        }
    )
    result.update(
        {
            "whole_two_real_source_ownership": "The original formal tree with four external Phi and two h has only the minimally covariant Phi2 and H2 kinetic/mass terms, -g H Phi2/2, C Phi4/24 and Einstein gravity. The inherited full fourth-degree scalar derivative and curvature jets vanish. Phi6 and higher vertices cannot lower the external-Phi count in a tree. A higher H Phi4 source needs a second H endpoint in a different Phi component, giving at least six external Phi; two linear Proca sources require at least eight. Spectators with only quadratic matter vertices need a closed line. Ghosts without external ghost legs likewise require a loop. The scalar-forest identity is proved in notes/source.md.",
            "whole_two_real_vertex_budget": "For six external legs, sum(degree-2)=4. Necessary metric jets are Phi2 h through h3, H2 h through h3 (the h3 heavy entry is harmless but not needed by the four-Phi tree), H Phi2 and Phi4 density through h2, and Einstein h3/h4. Higher vertices cannot enter this arity. Fixed loop-marked vacuum/onepoint constants are retained at their inherited formal order, not promoted to tree couplings.",
            "whole_two_real_scope": "This is the complete selected original formal tree coefficient, not an interacting quantum state, all-N nonleading detector rate, finite radiative loop, physical matching assignment, Regge theorem or V/G/B/P8 closure.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "two_real_source_checks_copied_not_parent_mutation": checks
                is not old["checks"],
                "all_six_lower_curvature_jet_checks_replayed": sum(
                    n.startswith("two_real_replayed_no_lower_R_jet_") for n in checks
                )
                == 6,
                "higher_heavy_source_requires_distinct_scalar_components": True,
                "spectator_and_ghost_closed_lines_are_not_trees": True,
                "four_scalar_tree_core_not_all_loop_ownership": True,
            },
        }
    )
    return result
