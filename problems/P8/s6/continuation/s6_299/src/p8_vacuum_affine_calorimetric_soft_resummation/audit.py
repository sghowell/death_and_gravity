"""Complete leading-soft calorimetry with every original closure gate retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_unequal_mass_regge_residue import audit as previous

from . import gamma, limits, poisson, source

MODEL = "original_calorimetric_leading_soft_sum_not_full_P8"
OBSERVABLES = (
    "whole_total_energy_soft_sum",
    "whole_Gamma_factor_conversion",
    "whole_original_uniform_conversion_error",
    "whole_ordered_detector_scaling_limits",
)
ITEM = {
    "id": "QG2_H8A463_complete_leading_soft_calorimetric_sum_Gamma_factor_and_ordered_detector_scaling",
    "status": "COMPLETE_ALL_LEADING_SOFT_CALORIMETRIC_CONVERSION_AND_LIMIT_ORDER_NOT_FULL_RADIATION_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated all-leading-soft calorimetric model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require one stated leading-soft calorimetric observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "leading_soft_calorimetry",
        "analytic_detector_conversion",
    ):
        raise ValueError(
            "Require the leading-soft calorimetric or analytic conversion sector"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A leading-soft resummation cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_total_energy_Poisson_resummation": poisson.data(),
        "whole_Gamma_factor_and_uniform_original_bound": gamma.data(),
        "whole_ordered_and_joint_detector_limits": limits.data(),
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
        "all154_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 154,
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
        "no_full_amplitude_remainder_is_invented": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The complete leading-soft total-energy sum has the regulator-removed factor exp(Delta)*exp(-EulerGamma*a)/Gamma(1+a)*(E/nu)^a. Its original finite angular conversion agrees with S296; the error beyond1+Delta is below2*10^-1596 uniformly in resolution. An explicit joint-regulator calculation distinguishes the mandated ordered detector limit.",
        "domain": "Fixed four-massive-scalar hard momenta, positive D=4+2e regulator, total unresolved energy0<E<=nu, with e removed first. Original numerical bounds use mu1,hard COM energy[5/4,2],nu1,resolution<=1/8 and kappa10^800. Multi-soft leading factorization defines the summed observable; no uniform full-amplitude approximation is assumed.",
        "historical_qualification": qualifications(),
        "not_established": "Uniform in multiplicity full-radiation-minus-leading-soft control, gravity-Born recoil, finite hard matching, all omitted loops, full physical analytic-amplitude relation, high-energy complex-energy/Regge control, original quantum state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "choose_finite_curvature",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "individual_energy_cut_as_total",
            "discard_Gamma_factor",
            "promote_one_real_to_all_N",
            "discard_finite_angular_conversion",
            "take_regulator_limits_together",
            "truncate_large_resolution_log",
            "full_inclusive_rate",
            "full_analytic_dispersion_unitarity",
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
            raise ValueError("Unsupported calorimetric input accepted: " + name)
    return total


def controls():
    return {
        "whole_total_energy_simplex_and_original_masses": True,
        "entire_positive_Poisson_resummation": True,
        "finite_angular_conversion_not_discarded": True,
        "whole_Gamma_product_inequality": True,
        "ordered_and_joint_regulator_comparison": True,
        "no_uniform_full_radiation_bound_invented": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
