"""Full homogeneous shear Hamiltonian, contact and modewise state remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_tensor_noise_response import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import control, hamiltonian, response, vertices

ITEM = {
    "id": "actual_CD_homogeneous_shear_Hamiltonian_contact_and_modewise_remainder",
    "status": "ACTUAL_SHEAR_MODEWISE_FINITE_AMPLITUDE_REMAINDER_NOT_CONTINUUM_SELF_ENERGY_OR_V_G_B",
}


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the fixed unit CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError(
            "Require actual fixed CD parent, mass1000 and unit observation slab"
        )
    return t, k, m, L


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A modewise shear response cannot close original P8")
    return True


def packets():
    return {
        "complete_constrained_Hamiltonian": hamiltonian.data(),
        "full_metric_vertices": vertices.data(),
        "actual_state_and_response": response.data(),
        "modewise_finite_amplitude": control.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, p in packets().items()
        for key, value in p["checks"].items()
    }
    rows.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_scoped_row_added": len(matching()) - len(previous.matching()) - 1,
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
            name + "_" + key: bool(v)
            for name, p in packets().items()
            for key, v in p.get("gates", {}).items()
        },
        "same_full_parent_on_unimodular_clock_shear": True,
        "temporal_constraint_not_extra_oscillator": True,
        "all_three_polarizations_and_mixing_retained": True,
        "exact_shear_exponentials_not_quadratic_only": True,
        "both_metric_vertices_and_contact_retained": True,
        "physical_Hamiltonian_current_sign_checked": True,
        "same_actual_all_order_initial_state": True,
        "prescribed_metric_family_not_quantum_initial_factorization": True,
        "Hadamard_propagation_not_a_full_quantum_metric_state": True,
        "exact_modewise_finite_amplitude_remainder": True,
        "all_momentum_formula_not_uniform_UV_majorant": True,
        "benchmark_band_not_a_physical_cutoff": True,
        "homogeneous_external_shear_not_all_spatial_response": True,
        "instantaneous_frequency_not_dynamic_gap_theorem": True,
        "continuum_contacts_and_renormalization_still_required": True,
        "full_feedback_inverse_and_remainder_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "Hamiltonian": "The exact constrained three-mode Proca Hamiltonian on homogeneous unimodular shear has positive mass, magnetic and electric exponential terms plus the unchanged temporal constraint term.",
        "response": "The exact covariance tangent and retarded commutator agree and retain the generally nonzero second-metric contact. Generic shear mixes the transverse and longitudinal coordinate modes.",
        "remainder": "The actual background energy propagator is below4 and gives an explicit all-momentum modewise finite-amplitude covariance Taylor remainder. A stated small-shear mode benchmark has relative interaction-picture remainder below1e-6.",
        "boundary": "The remainder majorant grows in momentum and does not justify a continuum renormalized stress-response limit, full tensor self-energy inverse or quantum background.",
        "remaining": "Uniform UV-subtracted response and metric contacts, spatial/mixed/nonlinear feedback, full interacting matching and physical heavy/cutoff bounds, and original V/G/B.",
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
        (0, modes.KAPPA * 2, modes.MASS, 1),
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
        raise ValueError(
            "Unsupported homogeneous shear response scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "constraint_and_three_modes_not_scalar_surrogate": True,
        "generic_shear_mode_mixing_retained": True,
        "contact_omission_is_a_negative_control": True,
        "actual_state_not_reset_to_W8": True,
        "finite_amplitude_tail_not_only_formal_coefficient": True,
        "benchmark_not_physical_cutoff": True,
        "no_unjustified_continuum_limit_or_full_feedback": True,
        "original_P8_not_closed": True,
    }
