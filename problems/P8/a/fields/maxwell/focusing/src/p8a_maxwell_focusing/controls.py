"""A smooth future-complete geometric countercontrol, not a Maxwell SEE model."""

import sympy as sp


def switch(value):
    value = sp.sympify(value)
    return 10*value**3-15*value**4+6*value**5


def identities():
    x, w = sp.symbols("x w", real=True)
    p = switch(x)
    residuals = {f"switch_endpoint_{endpoint}_jet_{order}":
                 sp.diff(p, x, order).subs(x, endpoint)-(endpoint if order == 0 else 0)
                 for endpoint in (0, 1) for order in range(3)}
    residuals.update({"positive_switch_slope": sp.expand(sp.diff(p, x)-30*x*x*(1-x)**2),
                     "second_jet_square_reduction": sp.expand(sp.diff(p, x, 2)**2-225*(2*x-1)**2*(1-(2*x-1)**2)**2),
                     "second_jet_maximum_certificate": sp.expand(sp.Rational(4, 27)-w*(1-w)**2
                         -(w-sp.Rational(1, 3))**2*(sp.Rational(4, 3)-w)),
                     "third_jet_quadratic": sp.expand(sp.diff(p, x, 3)-(60-360*x+360*x*x))})
    return residuals


def calibration():
    amplitude, width = sp.Rational(31, 20), sp.Rational(1, 1000)
    # Extend 1-p5 by1 before0 and0 after1. It is C2 with bounded weak
    # third derivative. Positive normalized convolution preserves these caps.
    derivative_caps = [amplitude, amplitude*sp.Rational(15, 8), amplitude*6, amplitude*60]
    theorem_caps = [sp.Integer(value) for value in (2, 4, 16, 96)]
    past_magnitude = amplitude*(1-10*width**3)
    delta = sp.Rational(1, 10**8)
    # In its constant past, beta I=0. The actual conformal Maxwell vacuum
    # would need this negative additional EED to satisfy SEE with Lambda0.
    required_source = 3*amplitude**2-delta*sp.Rational(31, 60)*amplitude**4
    return {"amplitude": amplitude, "mollifier_width_over_tau": width,
            "global_Hubble_jet_caps": derivative_caps,
            "cap_margins": [cap-value for cap, value in zip(theorem_caps, derivative_caps, strict=True)],
            "past_contraction_magnitude_lower": past_magnitude,
            "history_margin": past_magnitude-sp.Rational(3, 2),
            "second_switch_derivative_square_margin": 36-sp.Rational(100, 3),
            "future_static_by_time_over_tau": 1+width,
            "future_log_scale_factor_loss_upper": amplitude*(1+width),
            "reference_vacuum_required_sigma_at_Lambda_zero": required_source,
            "required_sigma_exceeds_one_by": required_source-1,
            "future_timelike_and_null_complete_geometry": True,
            "initial_pointwise_SEC_holds": False,
            "actual_allowed_Maxwell_SEE_solution": False,
            "geometric_assumptions_alone_imply_the_conclusion": False}
