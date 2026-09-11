"""Explicit external-theorem applicability without a numerical stress claim."""

from functools import cache

import sympy as s
from p8_vacuum_flat_dirac_production import audit as previous

from . import scattering, symbols

MODULES = (symbols, scattering)
ITEM = {
    "id": "Hadamard_regular_free_flat_SAT8_in_out_state_and_uniform_one_particle_norm_limit",
    "status": "HADAMARD_FOR_SPECIFIED_FREE_EXTERNAL_MASS_OPERATOR_USING_IDENTIFIED_EXTERNAL_THEOREM_NOT_INTERACTING_CURVED_STATE_OR_NUMERICAL_STRESS",
}


def application():
    return {
        "primary_source": "https://arxiv.org/pdf/2108.11955",
        "identified_version": "Gerard and Stoskopf, arXiv:2108.11955v2, 18 October 2021",
        "theorem": "Theorem 1.1; hypotheses H1-H3 in 3.1.3, H4 in 3.3; symbol definition in 4.2",
        "geometry": "The specified four-dimensional flat spacetime with Cauchy surface R^3 has the trivial spin structure and Euclidean bounded geometry.",
        "H1": "h_out=h_in=identity, c_out=c_in=1, including bounded inverse metric and lapse.",
        "H2": "h-h_asym=0, shift=0, c-c_asym=0, so the geometric tails satisfy every required symbol bound.",
        "H3": "Our all-order proof gives the real source mass -M(t) minus its asymptote in S^-8; all spatial derivatives vanish.",
        "H4": "The asymptotic Fourier Hamiltonians have eigenvalues +/-sqrt(p^2+M_asym^2), with M_asym>0. Their spectra omit zero.",
        "state_identification": "Our norm limits of transported asymptotic projectors agree with the theorem's limits; the checked Hamiltonian-sign dictionary exchanges spectral labels.",
        "conclusion": "The theorem therefore applies to these specified free in/out states: they are pure quasifree Hadamard states.",
        "quantitative_boundary": "The external theorem is used, not re-proved or formalized here. It supplies no numerical stress or Hadamard-remainder seminorm bound in this application.",
    }


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("The fixed state and primitive ledger must not change")
    return True


def validate_application(candidate):
    if candidate != application():
        raise ValueError(
            "The identified external theorem and every application hypothesis must be retained"
        )
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(value)
        for mod in MODULES
        for key, value in mod.data()["checks"].items()
    }
    out.update(
        {
            "nine_prior_primitive_rows_retained": len(frontier()) - 9,
            "all_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_prior_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "only_named_free_Hadamard_item_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return out


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **symbols.data()["gates"],
        **scattering.data()["gates"],
        "all_order_symbol_tail_not_finite_derivative_sampling": True,
        "compact_time_derivatives_finite_for_each_order": True,
        "both_active_profile_signs_and_constant_inert_masses": True,
        "bounded_geometry_and_trivial_spin_structure_of_R3": True,
        "external_theorem_and_hypotheses_identified": True,
        "physical_and_source_Dirac_operators_coincide": True,
        "source_spectral_labels_exchanged_explicitly": True,
        "operator_norm_limit_is_first_quantized": True,
        "adjoint_limits_preserve_unitarity_and_projector_purity": True,
        "finite_energy_not_used_as_Hadamard_proof": True,
        "numerical_stress_and_interacting_curved_state_not_inferred": True,
        "all_frozen_scientific_sources_and_reports_unchanged": True,
    }


def bad_cases():
    invalid = (
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
    for i, value in enumerate(invalid + (-1, s.Integer(2))):
        out.append((f"derivative_order_{i}", symbols.tail_coefficient, (value,)))
        out.append((f"polynomial_order_{i}", symbols.polynomial, (value,)))
    for call, name, good in (
        (symbols.mass_derivative_tail, "derivative_tail", [1, 1, 1, 2]),
        (scattering.norm_caps, "norm_tail", [1, 1, 2]),
        (scattering.asymptotic_masses, "gap", [2, 1]),
    ):
        start = 1 if name == "derivative_tail" else 0
        for slot in range(start, len(good)):
            for i, value in enumerate(invalid):
                args = good.copy()
                args[slot] = value
                out.append((f"{name}_type_{slot}_{i}", call, tuple(args)))
    for args in ((1, -1, 1, 2), (1, 1, 0, 2), (1, 1, 1, 1)):
        out.append((f"derivative_domain_{args}", symbols.mass_derivative_tail, args))
    for args in ((-1, 1, 2), (1, 0, 2), (1, 1, 0)):
        out.append((f"norm_domain_{args}", scattering.norm_caps, args))
    for args in ((0, 0), (2, -1), (2, 2), (2, 3)):
        out.append((f"gap_domain_{args}", scattering.asymptotic_masses, args))
    for key in application():
        changed = application()
        changed[key] = "unsupported conclusion"
        out.append((f"application_{key}", validate_application, (changed,)))
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
        raise ValueError("Unsupported Hadamard application input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "finite_derivative_checks_not_all_order_proof": True,
        "negative_source_mass_is_only_convention": True,
        "first_quantized_unitary_not_global_Fock_implementer": True,
        "Hadamard_smoothness_not_numeric_remainder_bound": True,
        "free_flat_state_not_interacting_curved_B_state": True,
        "old_energy_enclosure_and_all_frozen_inputs_unchanged": True,
        "external_theorem_not_claimed_locally_formalized": True,
        "original_P8_not_closed": True,
    }
