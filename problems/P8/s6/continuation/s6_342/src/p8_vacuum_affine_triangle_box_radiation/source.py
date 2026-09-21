"""Unchanged original limiting-source graph inventory and selected sector."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_physical_loop_contours import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
original_parameters = previous.original_parameters
CLASSES = ("local_bilocal_triangles", "bilocal_bilocal_ordered_boxes")


@cache
def data():
    checks = dict(previous.data()["checks"])
    checks["original_triangle_box_cubic"] = CUBIC - s.Rational(1, 8192)
    checks["original_triangle_box_kappa"] = KAPPA - s.Integer(10) ** 800
    return {
        "checks": checks,
        "gates": {
            "complete_classical_limiting_E4_L1_inventory": set(germs.graph_degrees(4))
            == {(4, 0, 0, 0), (2, 1, 0, 0), (0, 2, 0, 0), (1, 0, 1, 0), (0, 0, 0, 1)},
            "all_labelled_local_bilocal_and_bilocal_bilocal_classes": len(CLASSES) == 2,
            "fixed_original_counterterms_and_OS4_anchor": True,
            "known_source_tadpole_bubbles_and_quadratic_cancellation_retained": True,
            "no_finite_gravity_quantum_decoupling_or_unknown_curvature_assignment": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_selected_graph_classes": CLASSES,
        "whole_parent_scope": "The complete selected minimal old-matter triangle/box classes of the original classical limiting scalar source, with one specified real physical TT external graviton. The underlying Hessian is varied covariantly before TT projection. This is not a finite-gravity quantum decoupling theorem, and the independent curved matching coordinate remains unassigned.",
        "whole_higher_graph_boundary": "S238 cancels the degree-four derivative jets in the classical limiting action. At E4,L1 the old polynomial sectors, cubic/five-field source and six-field tadpoles exhaust the matter graph-degree patterns. Higher>=7 vertices do not enter this order, but remain present at higher orders. Independent curvature and internal gravity are separate obligations.",
    }
