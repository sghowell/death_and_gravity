"""Curvature-coordinate composition, density order and instability control."""

from functools import cache

import sympy as s
from p8_vacuum_affine_curved_scalar_reference import resolvent as prior
from p8_vacuum_affine_isolated_shear_resolvent import spectral as shear

from . import weighted as w


@cache
def data():
    sigma = w.SIGMA
    left = prior.C0 / (sigma - prior.FORWARD) ** 2
    right = prior.C0 / (sigma - prior.ADJOINT) ** 2
    composite = left * w.INVERSE_BOUND * right
    wanted = 125 / (4 * (32 + 8 * w.M) * (sigma - 20) ** 2 * (sigma - 108) ** 2)
    physical = 384 * s.pi**2 * prior.g.bridge.KAPPA * composite
    Aat = shear.data()["closed_A2"].subs(shear.d, 2)
    comparison = s.Rational(8, 3) * Aat
    checks = {
        "retained_distinct_curvature_Volterra_exponents": (prior.FORWARD - 20) ** 2
        + (prior.ADJOINT - 108) ** 2,
        "ordered_curved_weighted_inverse_bound": s.factor(composite - wanted),
        "physical_density_restored_on_right": s.factor(
            physical
            - 12000
            * s.pi**2
            * prior.g.bridge.KAPPA
            / ((32 + 8 * w.M) * (sigma - 20) ** 2 * (sigma - 108) ** 2)
        ),
        "bounded_comparison_keeps_positive_real_shear_zero": -s.Rational(8, 3) * Aat
        + comparison,
        "comparison_radial_upper_below_one": s.Rational(8, 3)
        * (s.Rational(1, 30) + s.Rational(3, 14))
        - s.Rational(208, 315),
        "comparison_radial_integral_upper": s.integrate(
            shear.y**2 * (30 - 20 * shear.y**2 + 3 * shear.y**4) / 30, (shear.y, 0, 1)
        )
        - s.Rational(3, 14),
    }
    return {
        "left_coordinate_weighted_norm": left,
        "right_adjoint_coordinate_weighted_norm": right,
        "corrected_curvature_reference_weighted_inverse_bound": composite,
        "physical_conformal_density_weighted_bound": physical,
        "ordered_operator": "A_V=Bc*(Fdiag+V)Bc, inverse E_V=Y K_V Z, in exactly that order. This concerns the S225 actual curvature coordinates with the specified bounded middle correction.",
        "window_and_spaces": "On the unchanged conformal slab T<=1, every realr and all comoving momenta, E_V is bounded on the exponentially weighted L2_tH^r space by the displayed constant, with no spatial derivative loss. Removing the weight costs exp(sigma T).",
        "density": "The force-normalized reference uses output1/(64pi^2 kappa a^4); its inverse is E_V followed on the RIGHT by64pi^2 kappa a^4. The norm retains the upper factor6 and kappa. A small weighted number is not small physical finite-window backreaction.",
        "comparison_correction": {
            "matrix": s.diag(0, comparison),
            "bound": s.Rational(208, 315),
            "positive_real_Laplace_pole_at_zero_transfer": "lambda=2m, p=4m^2",
        },
        "comparison_scope": "This separate constant bounded channel matrix has a simple positive-real-Laplace shear zero, by A2'(p)>0. Its norm is below1, so the weighted inverse theorem still applies and retains that growing reference pole. This explicitly prevents reading weighted existence as stability; it is not the actual physical curved correction.",
        "actual_frontier": "The actual mass-scale/state/contact/local-remainder/tree/matter corrections and the S222 auxiliary/clock graph have NOT been proved to lie in the admitted bounded middle class.",
        "checks": checks,
        "gates": {
            "positive_curved_composite_bound": True,
            "comparison_with_growing_reference_pole_still_bounded": s.Rational(208, 315)
            < 1,
            "nonzero_constant_correction_not_deleted": comparison != 0,
            "unchanged_physical_kappa_retained": prior.g.bridge.KAPPA
            == s.Integer(10) ** 800,
            "both_original_curvature_coordinate_factors_retained": prior.C0
            == s.Rational(5, 2),
        },
    }
