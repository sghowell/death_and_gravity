"""Uniform finite-parameter C0 physical response bound with causal memory.

All norms are specified in notes/proof.md. Nothing here is a local derivative
expansion, delta-convergence claim, rolling gap or cutoff.
"""

import sympy as sp

from .exact import R_MAX, nonnegative, parameters, radius


def calibration(r=R_MAX):
    r = radius(r)
    ratio = sp.Rational(101, 100)
    eta = sp.Rational(322, 3)*ratio*r**2
    feedback = 320*r**4*(10+60*r**2)
    l_over_r2 = 5*(1+64*r**2)/(1-feedback)
    l_error_over_r4 = 320*(1+60*r**2+360*r**4)
    f0 = 1+50*r**2+300*r**4
    q_error_over_r4 = 292*sp.Rational(101, 100)+sp.Rational(101, 75)*r**2*(10+60*r**2)*323
    exponent = sp.Rational(1, 2)+eta
    exp_upper = 1+exponent+exponent**2/(2*(1-exponent/3))
    return {"r": r, "rho_over_r": ratio, "heavy_perturbation_exponent": eta,
            "full_heavy_exponent": exponent, "full_heavy_exp_upper": exp_upper,
            "GH_over_r2": sp.Rational(101, 75),
            "u_GHprime_over_r2": sp.Rational(14, 3)*ratio,
            "Bstar_GH_over_r2": sp.Rational(188, 3)*ratio,
            "GH_minus_G0_over_r4": sp.Rational(16*161, 9)*ratio**2,
            "GL_over_r2": 4/(1-18*r**2), "feedback": feedback,
            "l_over_r2": l_over_r2, "l_error_over_r4": l_error_over_r4,
            "F0_over_S": f0, "Q_error_over_r4": q_error_over_r4,
            "g_error_over_r4": sp.Integer(915),
            "stiffness_lower_times_r2": sp.Rational(8000, 801)-161*r**2,
            "stiffness_endpoint_upper_times_r2": 10+161*r**2,
            "RMS_over_stiffness_proxy_lower": sp.Rational(6, 13)}


def physical_response_error(delta, source_sup, r=R_MAX, momentum_squared=1):
    _, r, _ = parameters(delta, r, momentum_squared)
    return 1000*r**4*nonnegative(source_sup, "source L-infinity norm")


def checks():
    c = calibration()
    margins = {
        "rho_factor": c["rho_over_r"]**2-(1+sp.Rational(1, 800)),
        "heavy_exponent": sp.Rational(2, 3)-c["full_heavy_exponent"],
        "heavy_propagator": 2-c["full_heavy_exp_upper"],
        "adjoint_heavy": 64-c["Bstar_GH_over_r2"],
        "kernel_difference": 292-c["GH_minus_G0_over_r4"],
        "light_Green": 5-c["GL_over_r2"],
        "feedback": sp.Rational(1, 10000)-c["feedback"],
        "full_light": 6-c["l_over_r2"],
        "light_error": 323-c["l_error_over_r4"],
        "approximate_forcing": sp.Rational(101, 100)-c["F0_over_S"],
        "heavy_error": 296-c["Q_error_over_r4"],
        "physical_error": 1000-c["g_error_over_r4"],
        "positive_stiffness": c["stiffness_lower_times_r2"]-9,
        "endpoint_stiffness": sp.Rational(169, 16)-c["stiffness_endpoint_upper_times_r2"],
    }
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A uniform causal/source bound failed")
    return margins
