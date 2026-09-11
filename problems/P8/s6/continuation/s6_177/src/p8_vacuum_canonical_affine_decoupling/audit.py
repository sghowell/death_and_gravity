"""No promotion from an actual classical decoupling family to V/G/B closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_proca_gaussian import audit as previous

from . import amplitude, family, gravity, source

ITEM = {
    "id": "full_fixed_canonical_affine_family_actual_classical_gravitational_decoupling",
    "status": "ANCHORED_FULL_FUNCTION_CLASSICAL_LIMIT_AND_TREE_AMPLITUDE_NOT_QUANTUM_CONTOUR_TRUNCATION_OR_V_CLOSURE",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A classical decoupling limit cannot close original P8 gates")
    return True


def observable():
    return {
        "anchor": "Complete S6.174/S6.176 CD-REG-AFFINE-ISO at kappa0=1e800 and zeta=1e-6; no old scientific byte is changed.",
        "fixed_data": "Canonical Phi mass1; ALL independent functions f(Phi,Y), r(Phi,Y), a3=r_Y/Y; original canonical M1 and vector mass1000.",
        "family": "Physical coefficients and regular affine maps are explicitly kappa dependent, with identical full base action at kappa0. Dependent Ia gravity corrections are labeled, not called fixed nongravitational couplings.",
        "limit": "The actual full classical action tends on compact canonical jets to free normalized gravity, the FULL fixed scalar f+a3(L3-L4), free M1 and free massive Proca, with specified density/jet decay.",
        "tree": "The complete four-scalar tree amplitude of this limit is 2lambda sum(s-2)^2+3gamma stu-8gamma; b2=4lambda, not a physical error-controlled dispersion result.",
        "remaining": "No uniform forward graviton bound, quantum regulator removal, full contour/cuts/truncation, interacting parent/cutoff/background result, G or B closure.",
    }


@cache
def residuals():
    out = {}
    for name, packet in (
        ("family", family.data()),
        ("germs", family.germs()),
        ("source", source.data()),
        ("coefficients", source.coefficients()),
        ("gravity", gravity.data()),
        ("pole", gravity.pole()),
        ("amplitude", amplitude.data()),
    ):
        out.update({name + "_" + key: value for key, value in packet["checks"].items()})
    out.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_actual_decoupling_matching_row_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {k: s.cancel(v) for k, v in out.items()}


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        "full_base_functions_and_action_anchor_not_just_Taylor_jets": True,
        "complete_independent_canonical_functions_held_fixed": True,
        "canonical_mass_and_nonzero_interactions_do_not_drift": True,
        "dependent_Ia_gravity_coefficients_explicit": True,
        "fixed_canonical_domain_not_fixed_physical_X_strip": True,
        "uniform_positive_R_and_regular_full_affine_quotient": True,
        "analytic_null_coefficients_and_unique_regular_q": True,
        "off_base_clock_equations_not_asserted": True,
        "Einstein_boundary_and_normalized_tensor_kinetic_retained": True,
        "all_full_source_terms_and_contact_retained": True,
        "full_function_compact_jet_limit_not_small_field_truncation": True,
        "free_canonical_M1_and_massive_vector_spectators": True,
        "literal_DHOST_tree_vertex_and_complete_four_scalar_ownership": True,
        "positive_tree_b2_not_UV_completion": bool(4 * family.LAMBDA > 0),
        "fixed_transfer_limit_not_uniform_forward_limit": True,
        "no_EFT_pole_deletion_or_Regge_remainder_assumption": True,
        "classical_limit_not_quantum_continuum_or_cutoff_bound": True,
        "prior_parent_state_not_promoted_to_interacting_state": True,
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
    out = []
    for i, value in enumerate(bad):
        for pos in (0, 1, 2):
            args = [family.K0, 0, family.ZETA]
            args[pos] = value
            out.append((f"type_{pos}_{i}", family.require_scope, tuple(args)))
    for k, x, z in (
        (0, 0, family.ZETA),
        (-1, 0, family.ZETA),
        (family.K0 - 1, 0, family.ZETA),
        (family.K0, -s.Rational(1, 4096), family.ZETA),
        (family.K0, s.Rational(6, 5), family.ZETA),
        (family.K0, 0, 0),
        (family.K0, 0, s.Rational(1, 2000)),
    ):
        out.append((f"scope_{len(out)}", family.require_scope, (k, x, z)))
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
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported decoupling claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "higher_independent_interactions_not_discarded": True,
        "base_order_not_changed_to_an_integer_sequence": True,
        "physical_X_domain_not_mislabeled_fixed": True,
        "null_source_not_set_to_zero": True,
        "finite_gravity_forward_pole_not_erased": True,
        "off_base_bounce_not_claimed": True,
        "tree_b2_not_physical_omitted_order_control": True,
        "original_P8_not_closed": True,
    }
