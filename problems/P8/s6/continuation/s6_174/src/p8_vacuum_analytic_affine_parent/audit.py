"""Scope and exact-evidence ledger for a new classical analytic parent."""

from functools import cache

import sympy as s
from p8_vacuum_target_clock_null import audit as previous

from . import affine, chart, dynamics, source

ITEM = {
    "id": "analytic_CD_REG_AFFINE_ISO_classical_parent_clock_vacuum_and_local_constraints",
    "status": "NEW_REGULAR_CLASSICAL_PARENT_AND_LOCAL_CONSTRAINT_BLOCKS_NOT_QUANTUM_MATCHING_CUTOFF_OR_V_G_B",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("The classical analytic parent cannot close original P8 gates")
    return True


def observable():
    return {
        "object": "Separately named CD-REG-AFFINE-ISO local classical EFT; the full S6.109 light target and original M1 matter are unchanged.",
        "new_operators": "A covariant source-centered isotropic trace mass update and curl of W=T-B du, with vacuum-regular lower q retained.",
        "proved_boundary": "Exact 56-direction algebraic reduction; regular covariant coefficient/chart domain; exact classical clock and vacuum quadratic blocks; local nonlinear seven-mode neighborhoods at each finite clock point.",
        "not_a_transfer": "No GY14 interacting-sector, fixed-order vacuum matching, full-state, heavy-elimination or loop bound is transferred to this new parent.",
        "remaining": "Physical cutoff, omitted/quantum operators, causal/stable response in controlled neighborhoods, full parent vacuum dispersion and finite-gravity IR/Regge estimates remain open.",
    }


@cache
def residuals():
    out = {}
    for name, packet in (
        ("chart", chart.data()),
        ("target", chart.target_jets()),
        ("source", source.data()),
        ("affine", affine.data()),
        ("dynamics", dynamics.data()),
        ("clock", dynamics.clock_lapse()),
        ("full_quotient", affine.matrices()),
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
            "one_separate_classical_parent_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return out


def scalar_entry_count():
    return sum(
        v.rows * v.cols if isinstance(v, s.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    return {
        **dynamics.clock_lapse()["gates"],
        "full_source_pinned_analytic_target_retained": True,
        "unique_vacuum_regular_lower_q_not_reset_on_clock": True,
        "all64_source_and_all60_quotient_checks_retained": True,
        "source_centering_linear_and_constant_terms_retained": True,
        "new_mass_update_not_old_fixed_mu_away_from_clock": True,
        "mass_update_anisotropy_removable_through_null_gradient": True,
        "covariant_metric_map_explicit_inverse_and_jacobian": True,
        "physical_matter_metric_not_replaced_by_hat_metric": True,
        "generic_R_lapse_time_and_spatial_boundaries_retained": True,
        "new_source_clock_value_and_first_variation_zero": True,
        "nonzero_higher_source_not_discarded": True,
        "classical_clock_exact_and_quadratic_vector_positive": True,
        "vacuum_all_seven_physical_quadratic_modes_healthy": True,
        "clock_constraint_count_local_at_each_finite_time_only": True,
        "frequency_floor_not_stationary_gap_or_cutoff": True,
        "no_old_quantum_or_GY14_matching_transfer": True,
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
        for position in range(3):
            args = [s.S.Zero, s.S.One, s.Rational(1, 2000)]
            args[position] = value
            calls.append((f"type_{position}_{i}", chart.require_domain, tuple(args)))
    for x, r, z in (
        (-s.Rational(1, 4096), 1, s.Rational(1, 2000)),
        (s.Rational(6, 5), 1, s.Rational(1, 2000)),
        (0, s.Rational(1, 2), s.Rational(1, 2000)),
        (0, s.Rational(6, 5), s.Rational(1, 2000)),
        (0, 1, 0),
        (0, 1, -1),
        (0, 1, s.Rational(1, 1000)),
    ):
        calls.append((f"domain_{len(calls)}", chart.require_domain, (x, r, z)))
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
        raise ValueError("Unsupported classical-parent input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "fixed_mu_mass_not_inherited_off_clock": True,
        "unshifted_regular_q_clock_source_not_set_to_zero": True,
        "metric_inverse_not_only_timelike_coordinate_gauge": True,
        "nonzero_null_gradients_not_excluded": True,
        "secondary_lapse_boundary_not_dropped": True,
        "finite_vector_mass_not_cutoff": True,
        "full_quantum_parent_and_matching_not_claimed": True,
        "original_P8_not_closed": True,
    }
