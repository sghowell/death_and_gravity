"""Existing paired results and exact extent of the unevaluated frontier."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_spectral_bubble import calibration as bubble
from p8_vacuum_fermion_spectral_bubble import ownership as old_bubble
from p8_vacuum_fermion_spectral_pole import calibration as pole
from p8_vacuum_fermion_spectral_pole import quadratic as old_quadratic
from p8_vacuum_light_pole import normalization as original_quadratic

from . import families


def ownership_frontier():
    rows = []
    for record in families.catalog():
        ident = record["id"]
        rows.append(
            {
                "id": ident,
                "status": "BOUNDED_PAIRED_SECTOR"
                if ident == "scalar_Phi2_W1_F0"
                else "LOCAL_QUARTIC_SUBSET_ONLY"
                if ident == "scalar_Phi4_W2_F0"
                else "UNEVALUATED",
            }
        )
    return rows


def validate_frontier(rows):
    if rows != ownership_frontier():
        raise ValueError(
            "The two-loop ownership frontier differs from its exact declared scope"
        )
    return True


def rows_at_degree(degree):
    if type(degree) is not int:
        raise TypeError("Require an integer background degree")
    if degree not in (0, 2, 4):
        raise ValueError("Only background degrees zero, two and four are catalogued")
    return [r for r in families.catalog() if r["external_background_degree"] == degree]


@cache
def data():
    q = old_quadratic.data()
    b = bubble.data()
    p = pole.data()
    actual = b["actual_reference_parameters"]
    frontier = ownership_frontier()
    W = original_quadratic.data()["actual_light_quadratic_Hessian_insertion"]
    checks = {
        "same_actual_complete_nonlocal_quadratic_Hessian": sp.factor(
            q["literal_full_quadratic_insertion"]
            - sp.trace(q["inserted_light_covariance"] * W) / 2
        ),
        "same_mixed_quadratic_single_covariance_position": q["placement_count"] - 1,
        "same_local_quartic_three_channels_two_positions": len(
            old_bubble.data()["placements"]
        )
        - 6,
        "one_completed_paired_quadratic_row": sum(
            r["status"] == "BOUNDED_PAIRED_SECTOR" for r in frontier
        )
        - 1,
        "one_partially_bounded_quartic_row": sum(
            r["status"] == "LOCAL_QUARTIC_SUBSET_ONLY" for r in frontier
        )
        - 1,
        "seven_primitive_rows_still_unevaluated": sum(
            r["status"] == "UNEVALUATED" for r in frontier
        )
        - 7,
        "all_nine_rows_have_exact_status": len(frontier) - 9,
        "same_four_point_subset_bound": b["literal_family_enclosure"][
            "forward_second_coefficient_upper"
        ]
        - 2 * actual["Y_upper"] * actual["L"] ** 2 / (144**2 * actual["mF"] ** 4),
    }
    return {
        "exact_primitive_frontier": frontier,
        "same_actual_full_reduced_quartic_action_Hessian": W,
        "inherited_local_quartic_subset_upper": b["literal_family_enclosure"][
            "forward_second_coefficient_upper"
        ],
        "inherited_complete_paired_quadratic_family": p[
            "quadratic_covariance_insertion_family_enclosure"
        ],
        "separate_uncomputed_counterterm_tasks": [
            "All remaining proper-subgraph counterterm insertions with epsilon dependence",
            "All required second-order local references in the declared interaction scheme",
            "Complete finite interaction and canonical-field conversion through order two",
        ],
        "nonclosure": "A complete list of primitive families is not their evaluation. One paired quadratic row is bounded and one quartic row has only its local-vertex subset bounded; seven primitive rows and the separate counterterm/conversion tasks remain.",
        "scope": "Do not sum the inherited subset estimates as if this were the full enlarged-model two-loop error. Higher-loop, V-contour, G, B and original P8 closure remain separate.",
        "checks": checks,
    }
