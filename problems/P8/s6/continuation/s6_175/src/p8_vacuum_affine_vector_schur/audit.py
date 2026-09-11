"""Evidence and scope ledger for exact sourced-vector integration."""

from functools import cache

import sympy as s
from p8_vacuum_analytic_affine_parent import audit as previous

from . import gaussian, homogeneous, lorentz, vacuum

ITEM = {
    "id": "same_classical_affine_parent_closed_source_lift_and_scoped_vector_Schur_control",
    "status": "EXACT_CLASSICAL_AND_COMMON_REGULATOR_IDENTITIES_FULL_EUCLIDEAN_AND_LEADING_LORENTZ_BOUNDS_NOT_FULL_QUANTUM_MATCHING",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Sourced-vector identities do not close original P8 gates")
    return True


def observable():
    return {
        "parent": "The SAME CD-REG-AFFINE-ISO of S6.174, with complete S6.109 analytic light target and original M1; no new action or transfer from GY14.",
        "classical": "Closed-source and all homogeneous classical lifts; causal fixed-background response for compact sources with the sourced divergence constraint.",
        "Gaussian": "Exact source contact and Schur identity at a specified common regulator and boundary prescription, not a renormalized determinant or a state-independent in-in stress bound.",
        "quantitative": "Complete nonlinear source L2 bound and an independent Euclidean positive action bound on the small vacuum class; leading S4-only Lorentzian degree-eight functional estimate below the mass pole.",
        "remaining": "Full-source real-time/state/contour control, interacting renormalized parent stress and response, physical cutoff and omitted orders, V/G/B and original P8 remain open.",
    }


@cache
def residuals():
    out = {}
    for name, packet in (
        ("operator", gaussian.data()),
        ("variation", gaussian.variation()),
        ("homogeneous", homogeneous.data()),
        ("vacuum", vacuum.data()),
        ("Lorentz", lorentz.data()),
    ):
        for k, v in packet["checks"].items():
            out[name + "_" + k] = v
    out.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_same_parent_scoped_response_row_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {
        name: value.applyfunc(s.cancel)
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for name, value in out.items()
    }


def scalar_entry_count():
    return sum(
        v.rows * v.cols if isinstance(v, s.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    return {
        **vacuum.data()["gates"],
        **lorentz.data()["gates"],
        "same_full_analytic_parent_source_and_physical_metric": True,
        "all_source_contact_and_divergence_terms_retained": True,
        "closed_source_lift_uses_compatible_vector_data": True,
        "CD_global_hyperbolicity_and_compact_source_hypotheses_checked": True,
        "common_regulator_translation_not_renormalized_stress": True,
        "retarded_equation_not_single_branch_variational_kernel": True,
        "full_nonlinear_source_not_assumed_bandlimited": True,
        "Euclidean_bound_not_Lorentzian_matching_transfer": True,
        "finite_germ_bound_uses_complete_propagator_only_on_its_support": True,
        "generic_operator_unboundedness_not_nonlinear_image_no_go": True,
        "no_qualitative_Green_to_Hadamard_or_stress_inference": True,
        "all_prior_frozen_scientific_bytes_unchanged": True,
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
    calls = []
    for i, value in enumerate(bad):
        for pos in range(3):
            args = [1024, vacuum.KAPPA, vacuum.MAX_ZETA]
            args[pos] = value
            calls.append(
                (f"small_type_{pos}_{i}", vacuum.require_small_domain, tuple(args))
            )
        for pos in range(2):
            args = [4, vacuum.MAX_ZETA]
            args[pos] = value
            calls.append(
                (f"band_type_{pos}_{i}", lorentz.require_finite_band, tuple(args))
            )
    for n, k, z in (
        (2, 1024, vacuum.MAX_ZETA),
        (5, 1024, vacuum.MAX_ZETA),
        (s.Rational(9, 2), 1024, vacuum.MAX_ZETA),
        (1024, 65535, vacuum.MAX_ZETA),
        (1024, vacuum.KAPPA, 0),
        (1024, vacuum.KAPPA, -1),
        (1024, vacuum.KAPPA, s.Rational(1, 1000)),
    ):
        calls.append(
            (f"small_domain_{len(calls)}", vacuum.require_small_domain, (n, k, z))
        )
    for b, z in (
        (-1, vacuum.MAX_ZETA),
        (4, 0),
        (4, -1),
        (4, s.Rational(1, 16)),
        (5, s.Rational(1, 16)),
    ):
        calls.append((f"pole_domain_{len(calls)}", lorentz.require_finite_band, (b, z)))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        calls.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        calls.append((f"matching_{i}", validate_scope, (frontier(), rows)))
    calls.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
        )
    )
    return calls


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported sourced-vector claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "free_Lorenz_constraint_not_imposed_on_sourced_vector": True,
        "source_squared_contact_not_dropped": True,
        "homogeneous_lift_not_asserted_for_arbitrary_vector_initial_data": True,
        "retarded_kernel_not_varied_as_symmetric_single_branch_action": True,
        "Euclidean_norm_not_transferred_to_real_time": True,
        "finite_germ_support_not_assigned_to_full_analytic_source": True,
        "generic_pole_packets_not_an_actual_source_image_exclusion": True,
        "original_P8_not_closed": True,
    }
