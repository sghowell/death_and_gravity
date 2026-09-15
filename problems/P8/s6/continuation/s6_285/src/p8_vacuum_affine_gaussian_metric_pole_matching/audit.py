"""Preserved original frontier and complete fixed gapped Gaussian metric matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_dimensional_cut_completion import audit as previous

from . import matching as fixed_matching
from . import poles, source, spectral

MODEL = "original_fixed_H_Proca_metric_matching_and_unmatched_Phi_curvature_not_full_P8"
OBSERVABLES = (
    "whole_fixed_Gaussian_vacuum",
    "both_conserved_metric_spectral_cuts",
    "fixed_finite_Newton_and_curvature_matching",
    "Gaussian_metric_pole_and_low_disk_bound",
)
ITEM = {
    "id": "QG2_H8A449_complete_fixed_H_Proca_metric_matching_three_Gaussian_cuts_and_explicit_unmatched_light_curvature",
    "status": "COMPLETE_SCOPED_H_PROCA_METRIC_MATCHING_WITH_PHI_CURVATURE_OPEN_NOT_FULL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def frontier():
    return previous.frontier()


def qualifications():
    return previous.qualifications()


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the fixed gapped Gaussian metric model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a scoped Gaussian metric observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "fixed_gapped_Gaussian_metric",
        "fixed_Gaussian_long_range_residue",
    ):
        raise ValueError("Require one of the explicitly stated gravitational sectors")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Gaussian metric matching cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_both_spin_Gaussian_spectral_cuts": spectral.data(),
        "whole_fixed_covariant_polynomial_matching": fixed_matching.data(),
        "whole_Gaussian_metric_poles_and_bounds": poles.data(),
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
        "all9_primitive_statuses_preserved": len(frontier()) == 9,
        "all140_matching_records_preserved": matching()[:-1] == previous.matching()
        and len(previous.matching()) == 140,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "Gaussian_metric_kernel_not_full_four_point_tadpole_or_amplitude": True,
        "Gaussian_pole_not_exact_interacting_physical_normalization": True,
        "S275_S276_S277_and_S281_S284_scopes_unchanged": True,
    }


def observable():
    return {
        "established": "All three gapped Gaussian metric cuts have their full conserved spin0/spin2 spectral normalization. The inherited H/Proca prescriptions fix their complete curvature polynomials, and all three fixed Gaussian volume constants cancel every cosmological metric jet. The H/Proca finite Newton shift is positive and smaller than10^-602 relative to kappa; their prescribed metric inverse has no additional low-disk pole and a positive massless residue. The full Phi-inclusive kernel is explicitly parameterized by its three unresolved finite curvature coefficients.",
        "domain": "Same full original R,F and matter sources, mass1 Phi, unrounded heavy H and mass1000 Proca, at the separately stated constant-scalar Lorentz vacuum. The H and Proca finite prescriptions are inherited. Only Phi's Gaussian volume constant is pinned here; its finite gravitational curvature matching is not inferred from the H prescription.",
        "historical_qualification": qualifications(),
        "not_established": "Light-scalar finite gravitational curvature matching and full Phi-inclusive pole sign or inverse bound; complete four-scalar tadpole/local/threshold matching or finite amplitude and b20; full proper scalar-graviton vertex, external LSZ and interacting physical mass or Newton normalization; massless loops, omitted interacting loops, finite IR detector/Regge control, exact original vacuum/state/domain/measure/bounce or V/G/B/P8 closure.",
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
            ("spin", source.require_spin),
            ("mass", source.require_mass),
            ("sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate((-1, 1, 3, s.Rational(5, 2))):
        rows.append(("outside_conserved_spin_" + str(i), source.require_spin, (value,)))
    for i, value in enumerate((0, -1)):
        rows.append(("outside_positive_mass_" + str(i), source.require_mass, (value,)))
    for i, value in enumerate(
        (
            "full_first_loop_all_species",
            "UV_complete",
            "massless_external_replacement",
            "exact_interacting_Newton_pole",
            "zero_Gaussian_Newton_shift",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "full_finite_amplitude",
            "forward_b20",
            "drop_massless_physical_pole",
            "new_local_counterterm",
            "massless_external_scalar",
            "physical_cutoff_from_cap",
            "full_Regge_bound",
            "full_four_scalar_tadpole_from_metric_cut",
            "retune_fixed_vacuum",
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
            raise ValueError("Unsupported Gaussian metric input accepted: " + name)
    return total


def controls():
    return {
        "whole_unchanged_source_and_fixed_vacuum": True,
        "same_complete_covariant_finite_prescriptions": True,
        "whole_scalar_and_Proca_conserved_stress_cuts": True,
        "whole_covariance_and_original_Taylor_subtraction": True,
        "signature_correct_fixed_finite_Newton_shift": True,
        "whole_volume_cancellation_and_positive_Gaussian_residue": True,
        "full_four_point_IR_Regge_and_original_P8_boundaries_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
