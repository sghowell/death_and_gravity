"""Read-only whole nonlinear spatial gauge and finite ghost-vertex report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_nonlinear_auxiliary_measure import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-spatial-gauge.json"
PARENT_SHA = "ae63ae3c213a5a5cd8b02174ec9eacfd4d19b3d0fb103c05aa221649ead4e05f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_spatial_gauge/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen full nonlinear auxiliary and boundary parent changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_257_full_retained_nonlinear_canonical_auxiliary_boundary_and_ancestry_rebuilt": PARENT_SHA,
        "same_current_functions_profiles_and_canonical_boundary_prescription": True,
        "local_spatial_gauge_not_a_replacement_fixed_quantum_preparation": True,
        "finite_ghost_vertices_not_completed_covariant_quantum_measure_or_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A whole spatial generator, ghost operator, local slice or scope gate failed"
        )
    packets = audit.packets()
    geometry = tuple(
        name
        for name in packets
        if name
        not in (
            "whole_off_gauge_three_ghost_Fourier_operator",
            "whole_finite_three_ghost_logdet_shape_vertices",
            "whole_second_order_gauge_restoration_and_volume_contact",
        )
    )
    return {
        "schema": 1,
        "claim": "P8-S6.258.COMPLETE_RETAINED_NONLINEAR_SPATIAL_DIRAC_GAUGE_OPERATOR_LOCAL_SLICE_AND_FULL_FINITE_THREE_GHOST_VERTICES_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "EXACT_WHOLE_SPATIAL_GENERATOR_FRAME_OFF_GAUGE_GHOST_AND_SECOND_ORDER_RESTORATION_WITH_WRITTEN_LOCAL_SLICE_AND_FINITE_DETERMINANT_VERTICES; NOT_GLOBAL_GAUGE_BRST_COVARIANT_QUANTUM_MEASURE_FIXED_MEAN_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/spatial.md",
            "notes/gauge.md",
            "notes/local.md",
            "notes/ghost.md",
            "notes/reference.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_spatial_generator_frame_operator_and_local_slice": serialize(
            {name: payload(packets[name]) for name in geometry}
        ),
        "whole_finite_ghost_vertices_and_second_order_reference_gauge_restoration": serialize(
            {name: payload(packets[name]) for name in packets if name not in geometry}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged retained current parent has the complete displayed spatial momentum generator and nonlinear conformal Dirac gauge operator, including all off-gauge lower terms. The full physical/hat conformal map leaves its shape density exactly invariant, and its linearization matches the existing scalar-volume plus TT reference. The entire nonzero-momentum symbol has positive symmetrizer and explicit inverse. A local flat-torus gauge slice exists with translations retained separately; on the exact slice a shape operator-norm distance at most 1/8 gives full weak coercivity at least 3/4 and an explicit volume-dependent inverse gap. All three ghosts remain in the finite determinant calculation: the entire exponential shape family has second logdet variation 3, including a mandatory second shape contact. Two distinct TT waves require a nonlinear gauge-restoring coordinate change, whose retained scalar-volume coefficient is 1/32 in the displayed mode. The mean-zero complement is not a Lie subalgebra, finite summation by parts is not a first-class regulator, and no global gauge, infinite-volume gap, BRST/covariant continuum quantum measure, affine/projective quantization, physical subtraction, interacting state, fixed mean, cutoff, UV/Regge or original P8 closure follows.",
        "not_established": [
            "A globally unique spatial gauge or removal of the residual translation group and its global constraints",
            "A uniform infinite-volume infrared inverse or replacement of the original fixed quantum preparation by the finite torus diagnostic",
            "A diffeomorphism/BRST-preserving finite regulator, completed covariant quantum measure or affine/projective quantum determinant",
            "Continuum or time-regulated ghost determinants, physical counterterms, interacting state or fixed quantum mean",
            "Identification of the nonlinear quantum reference action with an unforced classical reduction after discarding one-point source contacts",
            "A physical nonlinear Cauchy or bounce theorem, controlled Wilsonian matching, UV scattering, gravitational IR/Regge remainder or original V/G/B/P8 closure",
            "Formalization of the written local gauge implicit-function, weak-coercivity, elliptic regularity and finite-determinant arguments",
        ],
        "verification_boundary": "Exact full tensor/canonical identities and independent finite-parameter diagnostics support written local spatial-gauge and finite-ghost claims. No external GR CMC or negative-curvature theorem is imported into this parent. Native/direct/ordinary/CLI use original SymPy; only captured full regression uses the audited exact-GCD adapter. Earlier frozen sources, profiles and fixed preparation remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The whole spatial gauge and ghost report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.258 whole nonlinear spatial gauge and finite three-ghost vertices replay passed; full quantum measure, fixed mean and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
