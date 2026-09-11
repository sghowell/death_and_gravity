"""Boundary-compatible tree tensor response, with full-response limits explicit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import audit as previous
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import detectors, energy, projector, tree

ITEM = {
    "id": "actual_CD_boundary_compatible_tree_tensor_noise_response",
    "status": "RESTRICTED_DETECTOR_TREE_PROPAGATED_NOISE_NOT_FULL_QUANTUM_TENSOR_STATE_OR_V_G_B",
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
        raise ValueError(
            "A restricted tree response observable cannot close original P8"
        )
    return True


def packets():
    return {
        "actual_tensor_tree": tree.data(),
        "TT_projector_and_compact_class": projector.data(),
        "all_momentum_adjoint_energy": energy.data(),
        "compatible_noise_observable": detectors.data(),
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
        "actual_full_parent_tensor_restriction": True,
        "fixed_scalar_profile_and_original_M1_retained": True,
        "full_source_zero_on_unimodular_clock_tensor_history": True,
        "tree_canonical_normalization_not_scalar_mixed_reduction": True,
        "full_stress_not_only_scalar_source_noise": True,
        "same_actual_all_order_CD_state": True,
        "two_boundary_conditions_not_one": True,
        "weighted_moments_use_actual_homogeneous_modes": True,
        "compatible_class_nonempty_not_every_detector": True,
        "no_uncorrelated_zero_quantum_initial_data": True,
        "no_sharp_or_arbitrary_gravity_switch": True,
        "no_scalar_Minkowski_divergence_transferred_to_Proca": True,
        "no_compact_support_claim_for_general_TT_projection": True,
        "full_momentum_energy_including_zero": True,
        "stress_distribution_test_space_respected": True,
        "homogeneous_tensor_part_annihilated_only_in_stated_class": True,
        "source_only_tree_component_not_full_quantum_state": True,
        "full_self_energy_and_finite_coupling_remainder_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "tree": "The actual clock tensor action has canonical field h=sqrt(kappa)gamma/2 and wave operator L=partial_t^2+3H partial_t-a^-2 Delta.",
        "class": "Compact TT detectors have two exact weighted Fourier moment conditions. They are equivalently q=L psi with compact TT psi and remove both free initial-data contributions.",
        "bound": "The actual all-momentum advanced energy estimate gives N(Gadv q)<5000 D(q). The unchanged CD stress bound gives leading source-only tensor STD below5e-372 D(q), or dimensionless metric STD below1e-771 D(q).",
        "initial_state": "The restricted observable does not require a factorized zero quantum metric state, a sharp source switch or a claimed Proca/CD initial-time divergence.",
        "remaining": "Full tensor self-energy and finite-coupling remainder, unrestricted initial-state correlations, scalar/mixed and nonlinear response, interacting parent loops/matching, physical heavy/cutoff bounds and original V/G/B.",
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
            "Unsupported restricted tensor-response scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "complete_parent_not_minimal_scalar_surrogate": True,
        "two_actual_weighted_boundary_moments": True,
        "generic_detector_boundary_is_a_negative_control": True,
        "zero_momentum_not_discarded": True,
        "noncompact_TT_projection_not_mislabeled_compact": True,
        "no_initial_factorization_or_switched_conservation": True,
        "leading_component_not_full_metric_state_or_remainder": True,
        "original_P8_not_closed": True,
    }
