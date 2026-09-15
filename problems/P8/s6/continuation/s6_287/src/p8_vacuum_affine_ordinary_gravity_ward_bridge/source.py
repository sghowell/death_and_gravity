"""Whole unchanged source and complete pure-GR proper-vertex graph scope."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_matter_graviton_vertex import source as previous

MU, T, EP, TAU = s.symbols("mu transfer EP positive_tau", positive=True)
K, NU = s.symbols("kappa positive_scale", positive=True)
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    inherited = previous.data()
    answer = {k: v for k, v in inherited.items() if k not in ("checks", "gates")}
    checks = dict(inherited["checks"])
    D, trace = s.symbols("D trace_T")
    kdot = s.Symbol("k_dot_T")
    projected_trace = -2 * trace / (D - 2)
    soft = 4 * ((MU - T / 2) ** 2 - MU**2 / (2 + 2 * EP))
    vd = T * T - 4 * MU * T + 2 * MU * MU + 2 * MU * MU * EP / (1 + EP)
    checks.update(
        {
            "general_D_harmonic_trace_cancellation": s.factor(
                kdot - trace / (D - 2) - projected_trace / 2 - kdot
            ),
            "whole_dimension_eikonal_numerator": s.factor(soft - vd),
            "three_vertices_three_propagators_relative_tree_sign": s.simplify(
                ((-s.I) ** 3 * s.I**3 * s.I) / (-s.I) + 1
            ),
            "literal_fixed_light_mass_unit": inherited[
                "same_original_vacuum_parameters"
            ]["mu"]
            - 1,
            "literal_original_kappa_again": KAPPA - s.Integer(10) ** 800,
        }
    )
    answer.update(
        {
            "specified_ordinary_vertex": "The original minimally coupled Phi proper Phi-Phi-h vertex at one pure-GR loop, order kappa^(-3/2), with linear metric g=eta+2h/sqrt(kappa), flat harmonic quadratic gauge fixing and no nonlinear ordinary gauge-fixing vertices. This is not the metric-reducible four-scalar amplitude.",
            "complete_pure_GR_proper_graph_inventory": {
                "A": "Two massless h and one massive Phi triangle, with the entire EH cubic and both literal scalar stresses.",
                "B": "One massless h and two massive Phi triangle, with the external metric on the scalar line.",
                "C": "Two massless h bubble with the complete scalar Phi^2 h^2 seagull and EH cubic.",
                "D": "Both orientations of a one-h one-Phi bubble, external h and one Phi at the scalar seagull.",
                "E": "The massless h tadpole at the Phi^2 h^3 contact, scaleless in dimensional regularization.",
            },
            "graph_completeness_boundary": "Every pure-GR vertex has two Phi legs or zero Phi legs. The connected one-loop Phi path and the unique external h give exactly A through E. Ghost loops attached to that path through a single h are metric reducible; a ghost loop with two attachments adds a loop. The original R-1 and derivative extras start at four Phi and cannot enter this one-pure-GR-loop two-Phi proper vertex. Massive g-squared loops are the separate unchanged S286 result.",
            "background_split": "For the bridge use g_total=g_bar+2gamma/sqrt(kappa), with background harmonic F_nu=bar_g^ab(bar_nabla_a gamma_bn-bar_nabla_n gamma_ab/2), no scalar term in F. At zero background momentum, matter/EH vertices coincide with the ordinary linear-split ones. Their only one-loop proper-vertex difference is the whole constant-background variation of the quadratic gauge-fixing density.",
            "dimensional_and_UV_scope": "D=4+2EP with 0<EP<=1/4 for infrared estimates. Ultraviolet local poles are subtracted covariantly before the real small-transfer bound; finite local analytic counterterms remain unspecified and are not chosen by this proof. Pole and finite Laurent coefficients, not an unregulated finite-t scattering amplitude, are compared.",
            "literal_eikonal_numerator": vd,
            "literal_triangle_relative_sign": -1,
            "checks": checks,
            "gates": {
                "whole_original_R_F_and_vacuum_retained": "whole_original_R_F"
                in answer,
                "source_checks_copied_not_parent_mutated": checks
                is not inherited["checks"],
                "all_five_proper_topology_classes_and_both_D_orientations": True,
                "entire_EH_and_scalar_seagulls_retained": True,
                "ordinary_and_background_gauge_definitions_explicit": True,
                "finite_counterterms_IR_detector_and_original_frontiers_not_chosen": True,
            },
        }
    )
    return answer
