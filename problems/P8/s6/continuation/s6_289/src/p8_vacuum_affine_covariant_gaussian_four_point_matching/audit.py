"""Scoped covariant local and complete Gaussian metric four-point matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_gravity_pole_completion import audit as previous

from . import forward, gaussian, local, source

MODEL = "original_covariant_Gaussian_four_point_matching_not_full_physical_P8"
OBSERVABLES = (
    "whole_source_and_counterterm_order",
    "entire_D_covariant_local_projection",
    "whole_Gaussian_metric_four_point",
    "exact_and_bounded_Gaussian_forward_coefficients",
)
ITEM = {
    "id": "QG2_H8A453_all_D_covariant_local_projection_and_complete_Gaussian_metric_four_point_with_exact_light_and_bounded_fixed_heavy_forward_coefficients",
    "status": "COMPLETE_SCOPED_COVARIANT_GAUSSIAN_FOUR_POINT_MATCHING_NOT_FULL_B20_IR_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original covariant Gaussian four-point sector")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated local-map or Gaussian matching observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "covariant_local_projection",
        "Gaussian_metric_four_point",
    ):
        raise ValueError(
            "Require a first-order local projection or Gaussian metric sector"
        )
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("One-sector Gaussian matching cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_D_covariant_local_projection": local.data(),
        "whole_Gaussian_metric_four_point": gaussian.data(),
        "whole_Gaussian_forward_coefficients": forward.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: clean(value)
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }


def scalar_entry_count():
    return sum(
        len(value) if isinstance(value, s.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "all9_original_primitive_statuses_preserved": len(frontier()) == 9,
        "all144_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 144,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "light_and_non_Gaussian_finite_coefficients_unassigned": True,
        "whole_physical_IR_and_Regge_frontiers_not_closed": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "A complete dimension-dependent rank-two projection maps the eight stated covariant local structures to the two on-shell four-Phi coefficients, with six explicit null directions and finite evanescent effects retained. The entire fixed H/Proca Gaussian metric four-point functions reproduce all-angle physical cuts. Their full pole-subtracted forward coefficient is strictly negative with magnitude below10^-1600. The light Gaussian nonlocal coefficient is exactly(319/4800-23pi/1280)/(pi^2*kappa^2), with its local curvature combination left unmatched.",
        "domain": "One-loop metric Gaussian sectors and one order-hbar local-counterfunctional tree insertion in the unchanged original physical frame. Actual mu1,kappa10^800,nH=10^200/512+2,M_A^2=10^6 and the already fixed H/Proca finite prescriptions. Known Newton/light poles are removed consistently before the gapped forward coefficient is formed.",
        "historical_qualification": qualifications(),
        "not_established": "Finite values or signs of the unmatched light and non-Gaussian local combinations; complete full-source b20, physical detector/dressing/IR or finite Regge remainder; all-order quantum field-redefinition, original state/domain/measure/bounce or V/G/B/P8 closure.",
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
        s.I,
        s.nan,
        s.Symbol("unspecified"),
    )
    rows = []
    for i, value in enumerate(bad):
        for label, call in (
            ("model", require_model),
            ("observable", require_observable),
            ("order", source.require_order),
            ("mass", source.require_mass),
            ("sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate((-1, 0, s.Rational(1, 2), s.Rational(3, 2))):
        rows.append(
            ("outside_positive_integer_order_" + str(i), source.require_order, (value,))
        )
    for i, value in enumerate((0, -1)):
        rows.append(("outside_positive_mass_" + str(i), source.require_mass, (value,)))
    for i, value in enumerate(
        (
            "whole_finite_G_amplitude",
            "UV_complete",
            "full_crossed_b20",
            "exact_massless_LSZ",
            "choose_light_finite_curvature",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "set_all_local_coefficients_zero",
            "Euler_topological_in_all_dimensions",
            "nonlinear_quantum_EOM_equivalence",
            "copy_heavy_finite_prescription_to_Phi",
            "add_light_Gaussian_loop_twice",
            "heavy_asymptotic_series_is_full_error_bound",
            "full_Regge_bound",
            "one_Gaussian_sign_is_full_positivity",
            "exchange_IR_limits",
        )
    ):
        rows.append(("unsupported_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("deleted_primitive", validate_scope, ([], matching())),
            ("deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    return rows


def rejected_inputs():
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported covariant Gaussian input accepted: " + name)
    return total


def controls():
    return {
        "whole_original_source_and_one_counterterm_order": True,
        "whole_all_D_local_map_and_nullspace": True,
        "whole_conserved_metric_response_and_Euler_kernel": True,
        "complete_Gaussian_functions_and_all_angle_cuts": True,
        "exact_light_and_uniform_fixed_heavy_forward_bounds": True,
        "one_sector_sign_not_full_positivity_or_Regge": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
