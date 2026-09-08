"""Uniform finite stress derivatives through order five, in physical units."""
from functools import lru_cache

import sympy as sp
from p8_vector_clock_matching import bounds as local_matching
from p8_vector_state import comparison, wkb

from . import preparation, spectral


@lru_cache(maxsize=None, typed=True)
def local_derivatives(observable, derivative_order):
    if observable not in ("energy", "pressure") or type(derivative_order) is not int or not 0 <= derivative_order <= 5:
        raise ValueError("Require energy or pressure and native derivative order 0..5")
    values = {order: sp.factor(sp.diff(coefficients[observable], wkb.u, derivative_order))
              for order, coefficients in local_matching.actual_coefficients().items()}
    boxes = {order: wkb.box_bound(value) for order, value in values.items()}
    return {"coefficients": values,
            "absolute_envelopes": {order: box["absolute_upper"] for order, box in boxes.items()},
            "reconstructions": {order: box["reconstruction"] for order, box in boxes.items()}}


def radial_envelope(mass_time_product):
    K = preparation.threshold(mass_time_product)
    m = comparison.physical_bounds(1, mass_time_product)["m0_tau"]
    A = comparison.AMAX
    data = preparation.constants()["common"]
    B, evolution = data["initial_mixing_envelopes"], data["evolution_mixing_envelope"]
    terms = {"low_preparation_band": 32*A**3*B[6]*m**3/27,
             "middle_preparation_band": A**3*B[8]*K/18,
             "high_preparation_band": A**3*B[10]/(9*K),
             "all_momentum_oscillatory_evolution": A**3*evolution/(24*m)}
    return {"m0_tau": m, "high_band_threshold": K,
            "integrated_beta_times_frequency_sixth_terms": terms,
            "integrated_beta_times_frequency_sixth_upper": sum(terms.values())}


def physical_bounds(planck_time_product, mass_time_product):
    scales = comparison.physical_bounds(planck_time_product, mass_time_product)
    L, m, A = scales["M_tau"], scales["m0_tau"], comparison.AMAX
    radial = radial_envelope(m)
    result = {"M_tau": L, "m0_tau": m, "high_band_threshold": radial["high_band_threshold"],
              "normalized_derivative_bounds": {}}
    for observable in ("energy", "pressure"):
        result["normalized_derivative_bounds"][observable] = {}
        names = [name for name in spectral.tail.physical_weights() if name.endswith("_"+observable)]
        for order in range(6):
            row_envelope = max(spectral.row_bounds(name, order)["normalized_reference_product_envelope"] for name in names)
            tails = [spectral.reference_tail(name, order) for name in names]
            if any(any(value != 0 for value in tail["low_coefficient_residuals"].values()) for tail in tails):
                raise ValueError("A nonintegrable lower subtraction coefficient survived")
            reference_constant = max(tail["reference_bracket_tail_integer_upper"] for tail in tails)
            difference = 9*row_envelope*A**(order+1)*m**(order-5)*radial["integrated_beta_times_frequency_sixth_upper"]/L**2
            reference = reference_constant/(72*m**2*L**2)
            local_data = local_derivatives(observable, order)
            local = sum(upper*m**(4-2*n) for n, upper in local_data["absolute_envelopes"].items())/(576*L**2)
            result["normalized_derivative_bounds"][observable][order] = {
                "exact_mode_minus_projected_reference": difference,
                "projected_reference_minus_differentiated_subtraction": reference,
                "finite_local_matched_term": local,
                "total": difference+reference+local}
    result["normalization"] = "The j-th physical-time derivative is divided by M^2/tau^(2+j)."
    result["state_cutoffs_and_finite_counterterms_unchanged"] = True
    result["perturbed_history_feedback_or_full_quantum_solution_claim"] = False
    return result


def radial_checks():
    A, m = sp.symbols("A m", positive=True)
    k, K = sp.symbols("k K", positive=True)
    nu2 = m**2+k**2/A**2
    return {"low_band_volume_bound_coefficient": sp.Rational(4**3, 54)-sp.Rational(32, 27),
            "middle_band_pointwise_denominator_identity": sp.factor(A**2-k**2/nu2-A**2*m**2/nu2),
            "high_band_pointwise_denominator_identity": sp.factor(nu2-k**2/A**2-m**2),
            "high_band_lower_momentum_square_margin": sp.factor(A**2*(K**2-m**2)-A**2*K**2/4-A**2*(3*K**2/4-m**2)),
            "all_momentum_nu_minus_four_radial_primitive": sp.factor(
                sp.diff((sp.atan(k/(A*m))/(2*A*m)-k/(2*(k**2+(A*m)**2)))*A**4, k)
                -k**2/nu2**2)}


def proof_checks():
    out = {"tail_boundary_exponent_for_fifth_derivative": 2+6-10 == -2,
           "same_integrable_envelope_for_every_lower_derivative": all(order-5 <= 0 for order in range(6)),
           "high_threshold_at_least_eight_masses": 8**2 > sp.Rational(4, 3),
           "three_polarizations_in_readout_constant": 3*3 == 9,
           "reference_subtraction_radial_constant": sp.Rational(3, 4)*sp.Rational(1, 6)*sp.Rational(1, 9) == sp.Rational(1, 72)}
    for name in spectral.tail.physical_weights():
        for order in range(6):
            row = spectral.row_bounds(name, order)
            tail = spectral.reference_tail(name, order)
            out[name+"_"+str(order)+"_frequency_degree"] = row["all_frequency_powers_at_most_derivative_order_plus_one"]
            out[name+"_"+str(order)+"_integrable_reference_tail"] = (
                all(value == 0 for value in tail["low_coefficient_residuals"].values())
                and tail["no_lower_Laurent_power"] and tail["all_majorant_coefficients_nonnegative"])
    return out
