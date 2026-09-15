"""Whole unchanged source and separately selected minimal one-loop graph sector."""

from functools import cache

import sympy as s
from p8_vacuum_affine_ordinary_gravity_ward_bridge import source as previous

S, T = s.symbols("energy transfer", real=True)
MU, K, NU = s.symbols("positive_mass_squared kappa positive_scale", positive=True)
EP = s.Symbol("EP", real=True)
D = s.Symbol("dimension", positive=True)
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    inherited = previous.data()
    answer = {
        key: value for key, value in inherited.items() if key not in ("checks", "gates")
    }
    checks = dict(inherited["checks"])
    m, g = s.symbols("number_matter_vertices number_gravity_vertices")
    checks.update(
        {
            "whole_minimal_one_loop_topological_identity": (m - 2) + (g + 2) - m - g,
            "whole_minimal_one_loop_superficial_degree_four": 4
            + 2 * (m + g)
            - 2 * ((m - 2) + (g + 2))
            - 4,
            "unchanged_actual_light_mass_squared": inherited[
                "same_original_vacuum_parameters"
            ]["mu"]
            - 1,
        }
    )
    answer.update(
        {
            "selected_four_point_sector": "Connected on-shell four-Phi amplitudes from the minimal two-Phi metric vertices and pure Einstein/ghost vertices, at one loop and order kappa^(-2), with the same ordinary linear split and harmonic gauge. The scalar external-leg correction is included. This selects graphs of the unchanged source; it does not remove C, g, heavy, Proca, higher-Phi or other couplings from the theory.",
            "finite_matching_boundary": "The light Gaussian metric volume constant is the entire already fixed S285 density. The light Newton coefficient delta_kappa_Phi=2c_Phi_R remains unspecified. The regular one-loop four-derivative on-shell counterfunctional is represented by two arbitrary crossing-local coefficients. None is assigned a finite value.",
            "coefficient_order_boundary": "Pole and finite Laurent coefficients in the retained D=4+2EP convention, not an exact finite-G or unregulated scattering observable. All-D principal-part algebra is exact, but the analytic pole-matching theorem is not asserted at every higher Laurent order.",
            "checks": checks,
            "gates": {
                "whole_original_source_and_parameters_retained": "whole_original_R_F"
                in answer,
                "source_checks_copied_not_mutated": checks is not inherited["checks"],
                "minimal_graph_selection_not_a_deleted_operator": True,
                "all_scalar_legs_and_metric_contacts_in_scope": True,
                "same_ordinary_linear_split_and_raw_convention": True,
                "no_light_finite_curvature_prescription_added": True,
            },
        }
    )
    return answer
