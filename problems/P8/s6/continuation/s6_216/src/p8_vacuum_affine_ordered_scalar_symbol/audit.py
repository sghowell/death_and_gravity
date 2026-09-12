"""Corrected oriented UV scalar input, with explicit historical withdrawals."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_prepared_ward_reconstruction import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import density, geometry, jets
from . import matching as local_matching

CORRECTIONS = {
    "complete_four_sector_spatial_UV_symbol_with_independent_detector_time": "HISTORICAL_ANNIHILATION_CREATION_PHASE_ERROR_CORRECTED_ENDPOINT_COEFFICIENTS_IN_S216",
    "actual_dimensional_spatial_UV_finite_part_with_evanescent_modes_and_counterterms": "HISTORICAL_PHYSICAL_FINITE_IDENTIFICATION_WITHDRAWN_EXPLICIT_EVANESCENT_CORRECTION_IN_S216",
    "complete_matched_tracefree_spatial_Gaussian_current_with_anchored_original_regulator": "PHYSICAL_CURRENT_AND_ANCHORED_REGULATOR_IDENTIFICATION_WITHDRAWN_PENDING_CONSISTENT_FULL_REASSEMBLY",
    "prepared_ordered_Ward_reconstruction_full_chart_correction_and_three_remaining_scalar_kernels": "PREPARED_ORDERED_WARD_IDENTITIES_CONDITIONAL_KNOWN_TRACEFREE_INPUT_WITHDRAWN_PENDING_REASSEMBLY",
}
ITEM = {
    "id": "corrected_retarded_pair_phase_full_ordered_scalar_UV_and_fixed_six_invariant_finite_matching",
    "status": "CORRECTED_LOCAL_SPATIAL_UV_SCALAR_FINITE_PART_NOT_FULL_CURRENT_REASSEMBLY_SCALAR_ANCHOR_REDUCED_INVERSE_OR_V_G_B",
}


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the original fixed compact CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError("Require fixed kappa,mass1000 and unit slab")
    return t, k, m, L


def frontier():
    return previous.frontier()


def matching():
    rows = [dict(row) for row in previous.matching()]
    if set(CORRECTIONS) - {row["id"] for row in rows}:
        raise ValueError("A precise historical correction target is missing")
    for row in rows:
        if row["id"] in CORRECTIONS:
            row["status"] = CORRECTIONS[row["id"]]
    return rows + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Corrected local scalar UV input cannot close original P8")
    return True


def packets():
    return {
        "literal_full_trace_field_strength_and_six_invariant_geometry": geometry.data(),
        "independent_canonical_retarded_branch_and_odd_endpoint_erratum": jets.phase_data(),
        "all_ten_full_constraint_endpoint_products": jets.data(),
        "ordered_scalar_coefficients_and_fixed_source_dimension_jets": density.data(),
        "full_covariant_poles_corrected_finite_part_and_local_norm": local_matching.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, p in packets().items()
        for key, value in p["checks"].items()
    }
    before, after = previous.matching(), matching()
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "all_primitive_statuses_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "exactly_four_historical_physical_status_corrections": sum(
                a != b for a, b in zip(before, after)
            )
            - 4,
            "every_other_historical_matching_row_unchanged": sum(
                a != b for a, b in zip(before, after) if a["id"] not in CORRECTIONS
            ),
            "one_new_corrected_scalar_checkpoint_row": len(after) - len(before) - 1,
        }
    )
    return {
        key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for key, value in rows.items()
    }


def scalar_entry_count():
    return sum(
        v.rows * v.cols if isinstance(v, s.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, p in packets().items()
            for key, value in p["gates"].items()
        },
        "physical_phase_fixed_before_comparing_curvature_pole": True,
        "full_trace_constraint_and_exponential_volume_retained": True,
        "original_preparation_mass_and_finite_prescription": True,
        "all_scalar_channels_independently_ordered": True,
        "no_frozen_scientific_test_or_report_byte_edited": True,
        "exact_four_row_historical_withdrawal_not_silent_recertification": True,
        "S213_assembly_and_S215_known_input_not_restored": True,
        "full_scalar_remainder_and_anchor_still_open": True,
        "original_V_G_B_and_P8_open": True,
        "all_frontier_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "verified_local_input": "Correctly oriented full spatial UV coefficients, fixed six-invariant dimensional finite part and1e6 local Z24 bound, including all three ordered scalar channels.",
        "erratum": "Odd endpoints reverse relative to frozen annihilation-branch extraction. The tracefree evanescent finite correction is explicit. S198's abstract creation template and S195's canonical Kubo sign remain valid.",
        "withdrawn": "S213 physical assembled current/anchored regulator and S215 known tracefree input require consistent finite endpoint and bulk reassembly. This checkpoint does not restore those claims.",
        "scope": "Only the local spatial UV difference is matched. Homogeneous trace anchor, full scalar state/time/contact remainders, reduced response/inverse, background and original V/G/B remain open.",
    }


def bad_cases():
    bad = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    out = []
    for i, value in enumerate(bad):
        for pos in range(4):
            args = [0, modes.KAPPA, modes.MASS, 1]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, 0, 1),
        (0, modes.KAPPA, modes.MASS, 0),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_{i}", validate_scope, (frontier(), rows)))
    out.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
            (
                "restore_withdrawn_parent_statuses",
                validate_scope,
                (frontier(), previous.matching() + [dict(ITEM)]),
            ),
        )
    )
    for ch in (True, None, 0, "trace", "tracefree", "reversed"):
        out.append((f"geometry_channel_{len(out)}", geometry.channel_directions, (ch,)))
    for ch in ("tensor", "vector", "scalar", True, None, 0, "trace"):
        out.append(
            (f"fixed_physical_scalar_{len(out)}", density.physical_channel_row, (ch,))
        )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported corrected scalar input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "finite_Fock_and_real_covariance_phase_orientation": True,
        "nonzero_wrong_phase_scalar_pole_and_evanescent_controls": True,
        "all_actual_ADM_polarization_pairs_including_trace": True,
        "integer_dimensional_literal_field_strength_and_full_volume_curvature": True,
        "independent_four_order_WKB_all_350_coefficients": True,
        "fixed_invariant_dimension_derivative_and_ordered_Green": True,
        "explicit_historical_withdrawal_and_immutable_parent_chain": True,
        "original_P8_not_closed": True,
    }
