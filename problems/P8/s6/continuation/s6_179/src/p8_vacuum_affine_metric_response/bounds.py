"""Source-pinned physical continuum bound and new uncancelled chart norm."""

import json
from functools import cache
from pathlib import Path

import sympy as s
from p8_proca_nonlocal_response import verify as old_response
from p8_vacuum_affine_proca_gaussian import verify as old_state
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA

from . import chart

RESPONSE_SHA = "fea66d1775ce50b7ef941af5b0987269f86ab6f1e6a5a3a23fc53ca63ba3bb01"
STATE_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"


def read_pinned(path, sha):
    path = Path(path)
    if old_response.sha(path) != sha:
        raise ValueError("A source-pinned physical conditional response report changed")
    return json.loads(path.read_text())


@cache
def data():
    old = read_pinned(old_response.REPORT, RESPONSE_SHA)
    state = read_pinned(old_state.REPORT, STATE_SHA)
    # Only the explicitly named PHYSICAL VECTOR block, not the old
    # Gaussian-plus-fixed-profile or background-cancelled clock fields.
    d = old["new_complete_physical_vector_response_bounds"]
    components = {
        key: s.Rational(value)
        for key, value in d[
            "complete_metric_response_C10_to_C0_component_bounds"
        ].items()
    }
    state0 = state["specified_vector_stress_and_reference_clock_response"]["stress"][
        "bounds_by_time_derivative_order"
    ]["0"]
    rho = s.Rational(state0["absolute_energy_derivative_upper"]) / KAPPA
    P = s.Rational(state0["absolute_pressure_derivative_upper"]) / KAPPA
    M = chart.envelopes()["C10_input_upper"]
    R = max(components.values())
    # |w1|<=1/2, |w1+w2|=r(1-r)<=1/4, |v|<=3/2.
    contact_N = s.Rational(9, 2) * rho + s.Rational(15, 2) * P
    contact_v = s.Rational(33, 2) * P
    contacts = max(contact_N, contact_v)
    current = 3 * M * R + contacts
    parts = {
        component: {name: s.Rational(value) for name, value in row.items()}
        for component, row in d["nonlocal_components"].items()
    }
    locals_ = {
        key: s.Rational(value)
        for key, value in d["finite_local_component_bounds"].items()
    }
    checks = {
        key + "_exact_reassembled_physical_response": components[key]
        - sum(parts[key].values())
        - locals_[key]
        for key in components
    }
    checks["physical_response_joint_reconstruction"] = R - s.Rational(
        d["joint_complete_metric_C10_to_C0_upper_bound"]
    )
    checks["all_current_terms_retained"] = current - 3 * M * R - contacts
    return {
        "source_pinned_old_response_sha": RESPONSE_SHA,
        "source_pinned_current_state_sha": STATE_SHA,
        "actual_conditional_physical_stress_C10_to_C0_over_kappa": components,
        "continuum_nonlocal_response_parts": parts,
        "matched_physical_local_response_parts": locals_,
        "uncancelled_reference_stress_C0_over_kappa": {"energy": rho, "pressure": P},
        "clock_input_C10_constant": M,
        "clock_current_C0_contact_rows_over_kappa": {"N": contact_N, "vhat": contact_v},
        "clock_current_C10_to_C0_over_kappa": current,
        "physical_uniform_display": "5e-795 per unit joint C10 physical lapse/log-scale norm",
        "clock_current_uniform_display": "1e-784 per unit joint C10 clock-chart lapse/log-scale norm; includes the nonzero background contacts",
        "normalization": "The fixed kappa=1e800 and fixed reference volume a0^3 are used. No division by H or the zero bounce density, scalar counter-profile, field-dependent state reset or ultraviolet momentum cutoff is used.",
        "checks": checks,
        "gates": {
            "physical_response_positive": all(x > 0 for x in components.values()),
            "physical_response_below_five_e_minus_795": R < 5 * s.Rational(1, 10**795),
            "both_nonlocal_parts_retained": all(
                x > 0 for row in parts.values() for x in row.values()
            ),
            "physical_local_parts_retained": all(x > 0 for x in locals_.values()),
            "uncancelled_state_bounds_positive": rho > 0 and P > 0,
            "nonstationary_contact_rows_positive": contact_N > 0 and contact_v > 0,
            "complete_clock_current_below_one_e_minus_784": current
            < s.Rational(1, 10**784),
        },
    }
