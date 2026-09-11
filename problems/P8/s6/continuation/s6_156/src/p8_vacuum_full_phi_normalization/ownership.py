"""Complete quadratic ownership in a direct MS fermion representation."""

from functools import cache

from p8_vacuum_finite_field_covariance import audit as previous
from p8_vacuum_two_loop_light_pole import graphs

QUADRATIC_IDS = ("scalar_Phi2_W0_F2", "scalar_Phi2_W1_F0", "gauge_Phi2")


def rows():
    return [r for r in previous.frontier() if r["id"] in QUADRATIC_IDS]


def validate_quadratic_rows(given):
    if given != rows():
        raise ValueError("The direct-MS quadratic ownership rows changed")
    return True


@cache
def data():
    all_rows = previous.frontier()
    selected = rows()
    return {
        "all_fermionic_quadratic_rows": selected,
        "pure_scalar_raw_refinement_count": graphs.data()["raw_refinement_count"],
        "proper_counterterm_ownership": [
            "scalar Phi2 W0 F2 and gauge Phi2: complete proper fermion mass, kinetic and Yukawa MS insertions; both vertex anchors are MS analytic reference values, not canonical Yukawa parameter conversions",
            "scalar Phi2 W1 F0: entire physical inner fermion OS mass/residue grouping and fixed H source; local mass/source allocation has zero slope and zero outer OS remainder",
            "pure scalar 32: entire scalar inner OS grouping, full proper interaction-MS conversion, nested alpha term and all required scalar local insertions",
            "first fermion Phi4 counterterm inserted in a light tadpole and the whole fermion-cycle contractions are momentum independent, so no additional slope or nonlocal OS term",
            "outer Phi quadratic counterterms impose the same physical mass and residue; their finite slope is t_H before subtraction, not discarded",
        ],
        "direct_MS_coordinate_choice": "Fermion kinetic/mass/Yukawa and gauge parameters remain MS. Frozen proper zero-momentum anchors are analytic decompositions with their MS finite remainders restored, not additional coordinate redefinitions. The only change between this hybrid field and boundary MS Phi is the regulated Phi normalization.",
        "checks": {
            "exact_three_fermionic_quadratic_rows": len(selected) - 3,
            "all_nine_primitive_rows_retained": len(all_rows) - 9,
            "all_thirty_two_scalar_refinements": graphs.data()["raw_refinement_count"]
            - 32,
            "exact_quadratic_ID_partition": len(
                set(QUADRATIC_IDS) - {r["id"] for r in selected}
            ),
            "all_quadratic_rows_have_completed_primitive_bounds": sum(
                r["status"] == "OPEN" for r in selected
            ),
        },
        "scope": "Completeness is for the two-point slope and OS remainder in the direct-MS representation, not an unqualified summation of primitive bounds. Assigned first counterterms are removed from separate insertion lists. Vacuum/source values and the four-point matching assembly remain separate.",
    }
