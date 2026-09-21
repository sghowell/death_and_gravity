"""Unchanged original source; a known loop coefficient is not parent matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_curvature_contact_basis import source as previous
from p8_vacuum_affine_heavy_parent_one_loop import germs

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
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
        checks["finite_triangle_original_" + name] = s.factor(got - want)
    return {
        "checks": checks,
        "gates": {
            "same_original_masses_couplings_and_physical_metric": True,
            "specified_triangle_class_not_all_curved_matching": True,
            "fixed_finite_OS4_constant_and_counterterms_unchanged": True,
            "off_shell_analytic_jet_not_a_massless_physical_limit": True,
            "no_addition_of_an_already_counted_loop_to_S342": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Compute one finite known selected-matter triangle coefficient in a fully stated local comparison convention. Neither the source nor the independent extra parent curvature coefficient is assigned or changed. The full triangle is already present in S342 and is not counted twice.",
    }
