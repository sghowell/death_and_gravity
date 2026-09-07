"""Continuous anchored C3 comparison to the actual A18 radiation history."""

import sympy as sp
from p8a_maxwell_cosmology import geometry as prior_geometry

from . import dynamics


def response_coefficients():
    r = prior_geometry.RATIO
    zero = 64*r  # 32*r/(1-8*lambda)<=64*r on lambda<=1/16.
    scheme = (sp.Integer(64), sp.Integer(6144), sp.Integer(589824))
    lipschitz = (sp.Integer(8), sp.Integer(96), sp.Integer(1536))
    return [zero, *(a+b*zero for a, b in zip(scheme, lipschitz, strict=True))]


def history(delta):
    """Uniform all-x errors, not sampled derivative values or asymptotic O()."""
    lam = dynamics.dimensionless_coupling(delta)
    coefficients = response_coefficients()
    errors = [value*lam for value in coefficients]
    margins = [cap-error for cap, error in zip(prior_geometry.ERROR_CAPS, errors, strict=True)]
    if min(*margins) <= 0:
        raise ValueError("the actual thermal history is not inside the pinned A18 C3 tube")
    return {"delta": sp.sympify(delta), "lambda": lam,
            "history_interval_x": [-prior_geometry.RATIO, sp.Integer(0)],
            "reference_power": sp.Rational(1, 2),
            "backward_solution_y_bounds": [sp.Rational(50, 27), sp.Integer(2)],
            "C3_response_coefficients_per_lambda": coefficients,
            "C3_error_upper": errors, "strict_A18_tube_margins": margins,
            "actual_observer_H0_times_tau": sp.Integer(2),
            "anchored_at_observer": True,
            "a0_one_is_a_chosen_spatial_normalization": True,
            "actual_state_and_SEE_history_realized": True,
            "A18_future_bounds_proved_beyond_this_branch_endpoint": False}


def identities():
    y, lam, z = sp.symbols("y lambda z", nonnegative=True)
    d = 1-2*lam*y*y
    zz = lam*y*y
    jets = dynamics.jet_functions(y, lam)
    differences = [2*lam*y**4/d,
                   8*lam*y**5*(6*zz*zz-8*zz+3)/d**3,
                   16*lam*y**6*(6*zz*zz-5*zz+2)*(14*zz*zz-22*zz+9)/d**5]
    result = {f"exact_scheme_jet_difference_{j+1}": sp.factor(actual-base-correction)
              for j, (actual, base, correction) in enumerate(zip(jets[1:], (2*y*y, 8*y**3, 48*y**4), differences, strict=True))}
    result.update({"reciprocal_clock_difference": sp.simplify(-dynamics.velocity(y, lam)/y**2+2+2*lam*y*y/d),
                   "response_coefficient_arithmetic": sum((a-b)**2 for a, b in zip(response_coefficients(),
                       (sp.Rational(16, 25), sp.Rational(1728, 25), sp.Rational(155136, 25), sp.Rational(14770176, 25)), strict=True)),
                   "first_numerator_upper": sp.expand(3-(6*z*z-8*z+3)-2*z*(4-3*z)),
                   "second_numerator_upper": sp.expand(2-(6*z*z-5*z+2)-z*(5-6*z)),
                   "third_numerator_upper": sp.expand(9-(14*z*z-22*z+9)-2*z*(11-7*z))})
    return result


def calibration():
    z = sp.Rational(1, 4)
    return {"actual_history_at_delta_upper": history(dynamics.DELTA_MAX),
            "generic_small_coupling_upper": sp.Rational(1, 16),
            "strict_named_coupling_margin": sp.Rational(1, 16)-dynamics.dimensionless_coupling(dynamics.DELTA_MAX),
            "scheme_only_jet_error_coefficients": [sp.Integer(64), sp.Integer(6144), sp.Integer(589824)],
            "radiation_jet_Lipschitz_coefficients": [sp.Integer(8), sp.Integer(96), sp.Integer(1536)],
            "numerator_lower_values_at_one_quarter": [6*z*z-8*z+3, 6*z*z-5*z+2, 14*z*z-22*z+9],
            "strict_backward_y_above_one": sp.Rational(50, 27)-1,
            "continuous_interval_not_a_sample_scan": True}
