"""Conditional compact low-window rate estimates; no physical coefficient bound."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_real_remainder import born
from p8_vacuum_affine_dimensional_real_remainder import bounds as old_bounds

from . import contact, source

CONTACT_BOUND = s.Integer(1536)
TREE_REMAINDER_BOUND = s.Integer(4000000000)


def tensor_upper(omega):
    omega = old_bounds.resolution(omega)
    return CONTACT_BOUND * omega**2


def low_interference_upper(ss, z, resolution, coefficient):
    coefficient = contact.exact_coefficient(coefficient)
    ss, z = born.domain(ss, z)
    delta = born.transfer_gap(ss, z)
    x = old_bounds.resolution(resolution)
    if x > delta / 192:
        raise ValueError("Require resolution within the actual low recoil window")
    A0 = born.original(ss, z, 0)
    return (
        abs(coefficient)
        * (
            128 * CONTACT_BOUND * x**3 / 3
            + CONTACT_BOUND * TREE_REMAINDER_BOUND * x**4 / (2 * delta)
        )
        / (4 * s.pi**2 * source.KAPPA * A0)
    )


def contact_square_upper(ss, z, resolution, coefficient):
    coefficient = contact.exact_coefficient(coefficient)
    ss, z = born.domain(ss, z)
    x = old_bounds.resolution(resolution)
    A0 = born.original(ss, z, 0)
    return (
        coefficient**2 * CONTACT_BOUND**2 * x**6 / (24 * s.pi**2 * source.KAPPA * A0**2)
    )


@cache
def data():
    w, x, d = s.symbols("omega resolution delta", positive=True)
    C, B = CONTACT_BOUND, TREE_REMAINDER_BOUND
    rows = {
        "pair_spatial_dyad_norm": s.Integer(4 * 2 + 4 * 2 - 16),
        "six_pair_transverse_Frobenius_budget": s.Integer(6 * 16**2) - C,
        "entire_linear_interference_radial_integral": s.integrate(
            w * 2 * C * w**2 * (64 / w + B / d), (w, 0, x)
        )
        - (128 * C * x**3 / 3 + C * B * x**4 / (2 * d)),
        "entire_contact_square_radial_integral": s.integrate(w * C**2 * w**4, (w, 0, x))
        - C**2 * x**6 / 6,
        "linear_cutoff_tail_vanishes": s.limit(
            128 * C * x**3 / 3 + C * B * x**4 / (2 * d), x, 0
        ),
        "quadratic_cutoff_tail_vanishes": s.limit(C**2 * x**6 / 6, x, 0),
    }
    return {
        "checks": {key: s.factor(value) for key, value in rows.items()},
        "gates": {
            "uses_actual_original_recoil_energy_and_spatial_bounds": True,
            "both_polarizations_bounded_by_canonical_Frobenius_norm": True,
            "original_full_positive_Born_and_actual_phase_retained": True,
            "low_tree_estimate_not_used_beyond_delta_over192": True,
            "all_bounds_conditional_on_explicit_chi_not_physical_matching": True,
            "contact_square_not_promoted_to_complete_two_loop_rate": True,
        },
        "whole_uniform_contact_bound": "For5/4<=E<=2,0<=omega<=1/8, every original scalar has |p_spatial|<=2 and |k.p|<=4omega. Each v_ij has spatial norm<=16omega. Projection cannot increase it, so the sum of all six dyads has Frobenius norm<=1536omega^2. This is an analytic all-angle estimate, not inferred from finitely many states.",
        "whole_conditional_rate_bounds": "Let delta=min(1,-t,-u)>0 and A0 be the original full D4 Born. Atx<=delta/192, S335 gives ||Mtree||/A0<=(64/omega+4e9/delta)/sqrt(kappa). Cauchy-Schwarz with both TT modes, actual J<=1 and radial measure omega/(4pi^2) yields |chi|[128*1536*x^3/3+1536*4e9*x^4/(2delta)]/(4pi^2*kappa*A0) for the integrated absolute first-order interference. The contact square alone is <=chi^2*1536^2*x^6/(24pi^2*kappa*A0^2) for all x<=1/8. No value or physical bound on chi follows.",
    }
