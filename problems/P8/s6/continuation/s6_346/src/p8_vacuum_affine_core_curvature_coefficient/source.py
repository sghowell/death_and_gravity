"""Original source and pole-only derivative scheme retained without retuning."""

from functools import cache

import sympy as s
from p8_vacuum_affine_box_curvature_coefficient import source as previous
from p8_vacuum_affine_heavy_parent_one_loop import germs

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
original_parameters = previous.original_parameters


@cache
def data():
    checks = dict(previous.data()["checks"])
    for name, got, want in (
        ("mass", HEAVY_MASS2, germs.MASS2),
        ("g", CUBIC, germs.G),
        ("contact", CONTACT, germs.CONTACT),
        ("metric", KAPPA, germs.KAPPA),
    ):
        checks["common_core_original_" + name] = s.factor(got - want)
    return {
        "checks": checks,
        "gates": {
            "same_original_source_and_physical_metric": True,
            "fixed_OS4_value_counterterm_only_not_finite_derivative_retuning": True,
            "old_polynomial_sectors_only_new_source_classes_not_inferred": True,
            "analytic_origin_not_physical_above_threshold_approximation": True,
            "extra_parent_chi_unassigned": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_prescription": "Use H8A420-VAC-OS4 at mu1: the old polynomial bubble poles and on-shell quadratic conditions remain fixed, B_MS(0)=0, and the full finite OS4 adjustment is a constant quartic. No finite derivative subtraction is introduced. The comparison below is a change of operator coordinates, not a change of physical counterterms or metric.",
        "whole_sector_boundary": "Exactly the old polynomial A(v)^2 B(v)/2 bubble, g^2/4 times all24 dressed triangles, and g^4/4 times all24 alternating-mass boxes, with common1/(16pi^2). New local tadpole and mixed-source classes, internal gravitons and the independent extra parent curvature coefficient are not included.",
    }
