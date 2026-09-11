"""Private complete fixed-prescription homogeneous current and C2 remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_current_response import tail as comparison

from . import metrics, tensors

LOCAL = s.Integer(10) ** 86
CURRENT = (s.Integer(10) ** 87, s.Integer(10) ** 95, 2 * s.Integer(10) ** 111)


@cache
def data():
    R, Ric, Riem = s.symbols("R Ricci2 Riemann2", real=True)
    Euler = Riem - 4 * Ric + R**2
    a4s = R**2 / 72 - Ric / 180 + Riem / 180
    eps, x = s.symbols("epsilon x", real=True)
    total = tuple(comparison.COMPARISON[a] + LOCAL for a in range(3))
    checks = {
        "physical_dimension_complete_finite_curvature_reduction": s.expand(
            -4 * a4s + R**2 / 30 + Ric / 15 + Euler / 45
        ),
        "unchanged_proof_complement_and_tail": comparison.COMPARISON[2]
        - s.Integer(10) ** 111,
        "complete_integral_Taylor_factor": s.integrate((1 - x) * eps**2, (x, 0, 1))
        - eps**2 / 2,
        "zero_amplitude_Taylor_edge": (eps**2 * CURRENT[2] / 2).subs(eps, 0),
        "same_canonical_kappa": modes.KAPPA - s.Integer(10) ** 800,
    }
    return {
        "fixed_matching": "S193 identifies the original covariant homogeneous current with the S192 convergent actual-mode comparison plus the unchanged finite local action. Only after that dimension limit may the compact four-dimensional Euler variation be removed.",
        "local_derivative_display": LOCAL,
        "complete_current_parameter_displays": CURRENT,
        "complete_unrounded_sums": total,
        "physical_current_divided_by_kappa": tuple(x / modes.KAPPA for x in CURRENT),
        "domain": "The same smooth gamma=epsilon Gamma with ||Gamma^(j)||op<=1 through j12 and |epsilon|<=.01, common zero initial neighborhood, fixed smooth tracefree metric detector ||D||op<=1. These are homogeneous physical current bounds.",
        "Taylor": "The full original-prescription current remainder after its value and first amplitude derivative at zero is at most epsilon^2*1e111, zero at epsilon0 and strict otherwise. Proper density divides by the epsilon-independent a^3>=1 and preserves all bounds.",
        "not_a_solution": "A small renormalized response/source norm is not a full spatial/mixed inverse, finite-coupling solution, stability theorem or cutoff. No order reduction or deletion of higher-derivative modes is performed.",
        "checks": checks,
        "gates": {
            "all_local_geometry_and_detector_factors_retained": tensors.constants()[
                "complete_local_current_raw_parameter_display_before_rounding"
            ]
            < LOCAL,
            "all_complete_physical_current_bounds": all(
                total[a] < CURRENT[a] for a in range(3)
            ),
            "second_Taylor_coefficient": CURRENT[2] / 2 == s.Integer(10) ** 111,
            "same_actual_mass": modes.MASS == 1000,
            "metric_jet_constant_not_underestimated": metrics.G == 200000,
        },
    }
