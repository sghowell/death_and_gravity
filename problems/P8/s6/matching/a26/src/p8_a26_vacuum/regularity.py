"""Literal-coordinate regularity, with the full reconstructed-source cancellation."""

from functools import cache

import sympy as sp


def tilt_data():
    """Chosen audit value, not a numerical b attributed to the paper."""
    alpha = sp.Rational(509, 250)
    return {"b": sp.Rational(9, 500), "alpha": alpha, "n_s_from_source_equation_21": sp.Rational(24, 25),
            "numerical_b_printed_by_source": False,
            "Q_leading_absolute_power": -sp.Rational(1, 4),
            "Q_regular_derivatives": 2, "Q_second_holder_exponent": alpha-2,
            "Q_third_leading_coefficient": -alpha*(alpha-1)*(alpha-2)/4,
            "W2_regular_derivatives_if_Bprime_nonzero": 1,
            "W2_leading_sign_power_coefficient_in_Bprime_over_tau_B": alpha/4,
            "physical_instability_established": False}


def entropy_background_at_zero():
    return {"B": sp.Rational(29, 1000)+sp.Rational(7, 50)/sp.cosh(10),
            "tau_Bprime": sp.Rational(7, 500)*sp.tanh(10)/sp.cosh(10)}


def reconstructed_source(q, qp, b, bp, bpp, h, x, y, chi):
    """All Q/W2 terms of F-Xf1 and the interacting-chi source, not just W2."""
    w = -(q*bpp+(3*h*q+qp)*bp)/b
    return {"W2": w, "full_Q_part": q*(x*bp**2-y)+w*(b**2-chi**2)}


@cache
def checks():
    q, qp, b, bp, bpp, h, x, y, chi, xi, z, box = sp.symbols(
        "Q Qp B Bp Bpp H X Y chi xi gradphi_gradxi Boxphi", real=True)
    yy = sp.Symbol("gradxi_squared", real=True)
    source = reconstructed_source(q, qp, b, bp, bpp, h, x, y, chi)
    w, full = source["W2"], source["full_Q_part"]
    shifted = sp.expand(full.subs({chi: b+xi, y: bp**2*x+2*bp*z+yy}))
    expected_shift = -q*yy-2*q*bp*z-w*(2*b*xi+xi**2)
    # Integrate -2QB' partial(phi).partial(xi) once. No field equation was
    # used; divergence(QB' partial phi)=(Q'B'+QB'')X+QB' Box(phi).
    integrated = shifted.subs(z, 0)+2*xi*((qp*bp+q*bpp)*x+q*bp*box)
    expected = (-q*yy-w*xi**2
                +2*xi*((qp*bp+q*bpp)*(x+1)+q*bp*(box+3*h)))
    alpha = tilt_data()["alpha"]
    return {"all_Q_terms_cancel_on_reconstructed_trajectory": sp.simplify(full.subs({chi: b, y: bp**2*x})),
            "exact_entropy_shift": sp.expand(shifted-expected_shift),
            "one_IBP_full_linear_source": sp.simplify(integrated-expected),
            "linear_source_vanishes_on_clock_background": sp.expand(
                ((qp*bp+q*bpp)*(x+1)+q*bp*(box+3*h)).subs({x: -1, box: -3*h})),
            "constant_B_W2_control": sp.simplify(w.subs({bp: 0, bpp: 0})),
            "chosen_fractional_exponent": alpha-2*(1+sp.Rational(9, 500)),
            "chosen_b_spectral_dictionary": 1-2*sp.Rational(9, 500)*10/9-sp.Rational(24, 25),
            "finite_regular_Q2": alpha-2-sp.Rational(9, 250)}
