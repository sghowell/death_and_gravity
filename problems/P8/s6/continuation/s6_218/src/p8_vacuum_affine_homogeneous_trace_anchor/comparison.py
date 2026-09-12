"""Hypothesis-by-hypothesis transfer of full homogeneous Gaussian comparisons."""

from functools import cache

import sympy as s
from p8_vacuum_affine_current_response import adiabatic, low, tail, vertices
from p8_vacuum_affine_matrix_adiabatic import initial, riccati
from p8_vacuum_affine_matrix_response_tail import covariance, reference, variation
from p8_vacuum_affine_proca_gaussian import gaussian

from . import homogeneous


@cache
def data():
    dom = homogeneous.domination_data()
    old = tail.constants()
    xi = s.Symbol("xi", real=True)
    source = gaussian.data()["checks"]
    checks = {
        "same_original_proof_partition": initial.PARTITION - s.Integer(10) ** 16,
        "same_actual_initial_graph_error": initial.GRAPH_ERROR - s.Integer(10) ** 30,
        "same_complete_comparison_bounds": s.Matrix(
            [tail.COMPARISON[j] - s.Integer(10) ** (77 + 17 * j) for j in range(3)]
        ),
        "full_trace_vertex_bound_is_stronger": s.Matrix([2, 25, 650])
        - s.Matrix(vertices.G),
        "same_full_parent_clock_source": source["actual_covariant_clock_source_value"],
        "same_full_parent_first_source": source[
            "actual_covariant_clock_source_first_variation"
        ],
        "same_zero_mean_second_source_insertion": source[
            "linear_source_contact_has_zero_mean"
        ],
        "source_squared_mean_term_no_metric_Hessian": s.diff(xi**4, xi, 2).subs(xi, 0),
        "complete_integral_Taylor_factor": s.integrate(1 - xi, (xi, 0, 1))
        - s.Rational(1, 2),
    }
    return {
        "complete_hypothesis_transfer": "The new scalar-frequency and squeeze mixed jets are dominated entry by entry by S191, R=0, and the complete normalized trace current vertices are dominated through two amplitude derivatives by S192. The same positive Hamiltonian, unchanged pure Gaussian initial graph and zero parameter jets, frequency floor, fixed split, six-quadrature trace factor and even block-current parity apply. Reuse the generic finite Riccati/covariance/contour arguments only with these verified hypotheses; no old tracefree Hamiltonian identity is relabeled as a trace identity.",
        "finite_and_infinite_comparison": "The complete actual-state-minus-J_ad4 trace current, including the finite complement and full infinite tail, has coordinate-density C0,C1,C2 bounds(1e77,1e94,1e111) on this admitted scalar family. The common all-momentum bounds justify C2 amplitude differentiation into C0(I), not merely pointwise finite-k differentiation.",
        "same_reference_not_new_state": "The finite tenth-order comparison and fourth-order marker polynomial are mathematical subtractions. Their actual all-order state mismatch remains. No state reset, finite momentum grid replacing infinity or physical cutoff is introduced.",
        "contact_and_parity": "Every derivative of the full trace Hamiltonian and current includes the temporal-constraint and volume dependence. The physical metric vertex is real block diagonal, so only odd diagonal current marker coefficients vanish; odd QP covariance blocks are not discarded.",
        "parent_boundary": "S176's actual parent has source S=DS=0 at the reference and zero-mean second source insertion, so its Gaussian metric Hessian is the connected Proca response used here. Away from the reference the affine source need not vanish. The finite-amplitude C2 estimate here concerns the Gaussian determinant current component and is used to prove this reference Hessian; it is not a finite-amplitude remainder for the full sourced nonlinear parent.",
        "comparison_unrounded": old["full_actual_minus_fourth_order_comparison_bounds"],
        "checks": checks,
        "gates": {
            "all_new_domination_hypotheses": all(
                bool(v) for v in dom["gates"].values()
            ),
            "same_finite_Riccati_reference_smallness": riccati.data()["gates"][
                "reference_small_at_every_actual_momentum"
            ],
            "same_complete_mixed_reference_bounds": reference.data()["gates"][
                "same_reference_smallness"
            ],
            "same_actual_first_and_second_error": all(
                bool(variation.data()["gates"][k])
                for k in ("first_actual_error_display", "second_actual_error_display")
            ),
            "same_covariance_tail_estimates": all(
                bool(v) for v in covariance.data()["gates"].values()
            ),
            "same_complete_contour_parity_estimates": all(
                bool(v) for v in adiabatic.data()["gates"].values()
            ),
            "same_finite_band_estimates": all(
                bool(v) for v in low.data()["gates"].values()
            ),
            "all_complete_comparison_displays": all(
                old["full_actual_minus_fourth_order_comparison_bounds"][j]
                < tail.COMPARISON[j]
                for j in range(3)
            ),
            "no_finite_amplitude_full_parent_claim": True,
        },
    }
