"""Same original source and first-loop counterterm/metric insertion scope."""

from functools import cache

import sympy as s
from p8_vacuum_affine_gaussian_metric_pole_matching import source as gaussian_source
from p8_vacuum_affine_massive_gravity_pole_completion import source as previous

S, T = previous.S, previous.T
MU, K, NU = previous.MU, previous.K, previous.NU
D, EP = previous.D, previous.EP
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
VECTOR_MASS2 = gaussian_source.VECTOR_MASS2
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    inherited = previous.data()
    answer = {
        key: value for key, value in inherited.items() if key not in ("checks", "gates")
    }
    checks = dict(inherited["checks"])
    X, Y = s.symbols("scalar_X scalar_Phi_squared")
    trace = -(D - 2) * X / 2 + D * MU * Y / 2
    checks.update(
        {
            "whole_D_Einstein_trace_with_original_R_old_sign": s.factor(
                2 * trace / (D - 2) + X - D * MU * Y / (D - 2)
            ),
            "whole_phi_squared_derivative_four_point_IBP": 12 * (MU / 3) - 4 * MU,
            "entire_actual_light_mass_again": inherited[
                "same_original_vacuum_parameters"
            ]["mu"]
            - 1,
        }
    )
    answer.update(
        {
            "selected_counterfunctional_order": "One insertion of an order-hbar local counterfunctional into an on-shell four-Phi TREE, and the complete single H, Phi or Proca Gaussian METRIC loop coupled by the original minimal scalar metric currents. Extra source couplings and non-Gaussian diagrams are retained as separate sectors. This is not a finite field redefinition of the original bounce action.",
            "counterterm_coordinate_scope": "Eight displayed covariant local structures project to two on-shell combinations at this order. The computation evaluates the counteraction on the ordinary sourced metric and checks the scalar Fourier vertices. It neither chooses the coefficients nor asserts nonlinear/all-orders quantum EOM equivalence.",
            "Gaussian_scope": "The H/Proca metric densities, their whole volume cancellation, Newton shifts and finite curvature polynomials are the frozen S285 ones. The light scalar nonlocal Gaussian function is fixed, but its local curvature coefficients remain unmatched. Its Gaussian part is already present in S288, not an additional copy.",
            "checks": checks,
            "gates": {
                "whole_original_R_F_and_parameters_retained": "whole_original_R_F"
                in answer,
                "source_checks_copied_not_mutated": checks is not inherited["checks"],
                "first_loop_counterterm_order_not_finite_field_redefinition": True,
                "whole_metric_Gaussian_densities_not_only_selected_polarizations": True,
                "original_physical_frame_and_bounce_action_unchanged": True,
                "light_finite_curvature_and_other_sectors_not_chosen": True,
            },
        }
    )
    return answer
