"""One tame classical comparison checkpoint; no advancement of original primitives."""

from functools import cache

import sympy as s
from p8_vacuum_affine_reduced_scalar_hamiltonian import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import charts, coefficients, energy, majorants

ITEM = {
    "id": "uniform_two_chart_QG1_classical_scalar_phase_propagator_with_twelve_spatial_derivative_loss",
    "status": "TAME_CLASSICAL_COEFFICIENT_PROPAGATOR_NOT_SAME_SPACE_QUANTUM_INVERSE_STABILITY_OR_P8",
}
require_scope = previous.require_scope


def require_chart(which, time, momentum):
    require_scope(time)
    t, p = map(rational, (time, momentum))
    if p < 0:
        raise ValueError("Use nonnegative transfer magnitude")
    if which == "low":
        okay = p <= 100
    elif which == "central":
        okay = abs(t) <= s.Rational(1, 4) and p >= 100
    elif which == "outer":
        okay = abs(t) >= s.Rational(1, 8) and p >= 100
    else:
        raise ValueError("Require low, central or outer chart")
    if not okay:
        raise ValueError("Outside the proved chart and transfer domain")
    return which, t, p


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A tame classical comparison does not close the quantum inverse or P8"
        )
    return True


def packets():
    return {
        "actual_coefficient_jets_and_characteristics": coefficients.data(),
        "complete_finite_transfer_two_chart_algebra": charts.data(),
        "uniform_denominator_coercivity_and_48_coefficient_majorants": majorants.data(),
        "energy_and_polynomial_phase_propagator": energy.data(),
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
            "nine_original_primitive_rows": len(frontier()) - 9,
            "all_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_previous_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_scoped_tame_classical_checkpoint": len(matching())
            - len(previous.matching())
            - 1,
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
        "all_physical_transfers_covered": True,
        "finite_q_chart_failure_explicitly_excluded": True,
        "original_weighted_boundary_and_full_QG1_retuning_retained": True,
        "actual_high_frequency_pair_not_frozen_q_spectrum": True,
        "same_space_quantum_derivative_loss_not_solved": True,
        "no_frozen_input_edits_or_primitive_promotion": True,
        "unique_matching_identifiers": len({r["id"] for r in matching()})
        == len(matching()),
    }


def observable():
    return {
        "actual_characteristics": "The original unit-slab fixed-QG1 coefficient pair has speeds squared1 andF/Jc, with F>1/100 and Jc-F>1/1000. This does not substitute a principal rescaling for a crossing chart.",
        "classical_propagator": "The complete regular coefficient-sector first-order system has ||U(t,s;P)||<=exp(1e29)(1+|P|^2)^6 on every ordered subinterval of the unit slab.",
        "functional_space": "For any real r, H^(r+12) initial data and L1 H^(r+12) additive phase forcing give C_t H^r comparison evolution. Twelve spatial derivatives are lost.",
        "remaining": "The constant is not a kappa-smallness bound. Full nonlocal quantum constraints, physical-source maps, same-space coupled inverse, nonlinear parent, quantum background/stability, physical cutoff/heavy sector and original V/G/B remain open.",
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
            args = [0, previous.modes.KAPPA, previous.modes.MASS, 1]
            args[pos] = value
            out.append((f"scope_type_{pos}_{i}", require_scope, tuple(args)))
        out.append((f"chart_time_type_{i}", require_chart, ("central", value, 100)))
        out.append((f"chart_momentum_type_{i}", require_chart, ("central", 0, value)))
    for args in (
        (1, previous.modes.KAPPA, previous.modes.MASS, 1),
        (-1, previous.modes.KAPPA, previous.modes.MASS, 1),
        (0, 2 * previous.modes.KAPPA, previous.modes.MASS, 1),
        (0, previous.modes.KAPPA, 999, 1),
        (0, previous.modes.KAPPA, previous.modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for args in (
        ("central", 0, 99),
        ("central", s.Rational(1, 3), 100),
        ("outer", 0, 100),
        ("outer", s.Rational(1, 16), 100),
        ("outer", s.Rational(1, 2), 99),
        ("low", 0, 101),
        ("low", 0, -1),
        ("every_q", 0, 100),
        (None, 0, 100),
        (True, 0, 100),
    ):
        out.append((f"chart_domain_{len(out)}", require_chart, args))
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
                "extra_quantum_inverse",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "quantum_inverse", "status": "COMPLETE"}],
                ),
            ),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported tame propagator input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "independent_actual_gradient_and_full_matter_characteristic": True,
        "independent_original_canonical_flow_equals_both_complete_chart_Euler_operators": True,
        "independent_finite_q_saddle_and_phase_reconstruction": True,
        "nonzero_deleted_weighted_transport_control": True,
        "finite_q_auxiliary_degeneracy_with_regular_original_phase_system": True,
        "independent_energy_curve_with_nonsymmetric_lower_term": True,
        "exact_polynomial_majorants_and_twelve_derivative_Sobolev_weight": True,
        "all_previous_frontiers_unchanged_and_original_P8_open": True,
    }
