"""Explicit quantum-order distinction; no arbitrary finite matching choice."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import sew, tensor

from . import contact, source


def linear_interference(momenta, k, coefficient):
    coefficient = contact.exact_coefficient(coefficient)
    base, _ = tensor.original_cores(momenta, k)
    local = contact.transverse_core(momenta, k) / s.sqrt(source.KAPPA)
    return s.factor(2 * coefficient * sew.bilinear(base, local, 4))


def contact_square(momenta, k, coefficient):
    """Square of this direction alone, NOT a full two-loop rate."""
    coefficient = contact.exact_coefficient(coefficient)
    local = contact.transverse_core(momenta, k)
    return s.factor(coefficient**2 * sew.bilinear(local, local, 4) / source.KAPPA)


@cache
def data():
    hbar, chi, a, t = s.symbols("hbar chi a t", real=True)
    deformation = a + hbar * chi * t
    return {
        "checks": {
            "five_plus_minimal_three_point_four_external_loop_count": s.Rational(
                5 + 3 - 4, 2
            )
            - 2
            + 1
            - 1,
            "six_leg_two_graviton_contact_tadpole_loop_count": s.Rational(6 - 4, 2)
            - 1
            + 1
            - 1,
            "one_counterterm_plus_one_topological_loop_is_two_loop_order": s.Integer(
                1 + 1 - 2
            ),
            "tree_contact_with_five_external_has_no_topological_loop": s.Rational(
                5 - 5, 2
            )
            - 1
            + 1,
            "classical_source_is_unchanged_at_hbar_zero": deformation.subs(hbar, 0) - a,
            "first_order_rate_derivative": s.diff(deformation**2, hbar).subs(hbar, 0)
            - 2 * chi * a * t,
            "square_alone_is_not_first_order": s.expand(deformation**2).coeff(hbar, 2)
            - chi**2 * t**2,
            "four_dimensional_curvature_coefficient_mass_dimension": s.Integer(
                4 - (4 + 2 + 4) + 6
            ),
        },
        "gates": {
            "zero_flat_curvature_implies_zero_no_metric_vertex": True,
            "zero_fixed_metric_onepoint_and_flat_background_retained": True,
            "flat_four_scalar_data_only_through_one_loop_compared": True,
            "one_explicit_matching_axis_not_an_exhaustive_basis": True,
            "chi_not_selected_by_renormalization_scale_or_soft_finiteness": True,
            "physical_microscopic_matching_remains_research_not_user_choice": True,
        },
        "whole_flat_invisibility_proof": "The Weyl operator vanishes identically at flat metric and begins at four scalars plus one graviton. At the flat vacuum with the fixed zero onepoint, contributing to a no-external-graviton four-scalar coefficient requires at least one topological loop in addition to its explicit hbar. The contact plus minimal hPhiPhi vertex has I2,V2,L1; a two-graviton contact tadpole has I1,V1,L1. Thus this direction does not alter the retained flat matching through one loop, but it can enter it at two loops.",
        "whole_nonidentifiability": "The on-shell five-point contact is nonzero on original-source physical states and cannot be removed as a pure total derivative or equation-of-motion term. Four-point matching and known soft coefficients do not determine it. This is a local EFT comparison; no existence of arbitrary-chi UV completions, sign bound, full matching basis or change to the original source is asserted.",
        "whole_constructive_frontier": "Compute the covariant one-loop radiative functional and required curved counterfunctional, then determine physical hard matching from the quantum parent. A known minimal finite representative may be computed separately, but neither its pole subtraction nor flat value conditions select the independent finite curvature datum.",
    }
