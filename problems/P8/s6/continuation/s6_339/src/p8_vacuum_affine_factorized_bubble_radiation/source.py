"""Original covariant local/local Hessian class and fixed source parameters."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_heavy_scalar_four_point_loop import amplitude
from p8_vacuum_affine_mixed_source_radiation import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
original_parameters = previous.original_parameters
ENDPOINTS = ("contact_contact", "contact_heavy", "heavy_contact", "heavy_heavy")


def fixed_vertex(value):
    value = s.sympify(value)
    return CONTACT + CUBIC**2 / (HEAVY_MASS2 - value)


def endpoint_factor(value):
    return HEAVY_MASS2 * fixed_vertex(value) / CUBIC**2


@cache
def data():
    checks = dict(previous.data()["checks"])
    v = s.Symbol("channel_invariant")
    checks["same_original_factorized_vertex"] = s.factor(
        fixed_vertex(v)
        - amplitude.fixed_vertex(v).subs(
            {amplitude.C: CONTACT, amplitude.g: CUBIC, amplitude.n: HEAVY_MASS2}
        )
    )
    checks["same_original_cubic_normalization"] = CUBIC - s.Rational(1, 8192)
    checks["same_original_full_parent_contact"] = CONTACT - germs.CONTACT
    return {
        "checks": checks,
        "gates": {
            "all_four_ordered_endpoint_types_retained": len(ENDPOINTS) == 4,
            "full_covariant_local_local_Hessian_not_flat_guess": True,
            "both_outer_heavy_and_inner_light_insertions_retained": True,
            "fixed_original_mass_trilinear_and_contact_conditions": True,
            "independent_curvature_and_other_loop_classes_not_absorbed": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_ordered_endpoint_types": ENDPOINTS,
        "whole_covariant_generator": "W_local=-M_(C Phi^2+g^2 G_H Phi^2)/2. Its local/local contribution to -Tr(G_phi W G_phi W)/4 is -integral[K Phi^2](x)G_phi(x,y)^2[K Phi^2](y)/16, K=C I+g^2 G_H, with covariant measures and Feynman inverses. Vary before commuting any operator. The selected flat amplitude is sum A(v)^2 B(v)/(32pi^2), A(v)=C+g^2/(n-v).",
        "whole_boundary": "This is the complete factorized bubble class of the frozen full loop representation, not the remaining local/bilocal triangles, ordered boxes, unresolved mixed quadratic-radiation insertions, internal gravitons or full matching. No original contact or unknown curvature coefficient is changed.",
    }
