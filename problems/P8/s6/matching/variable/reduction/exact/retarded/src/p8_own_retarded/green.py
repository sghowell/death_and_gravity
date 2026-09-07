"""Finite-window retarded own-f Green bounds, with exact rational margins.

The prescribed physical metric is an off-shell input, not a separately
conserved matter source or an on-shell g solution.  The interval is in
x=T/(tau*sqrt(delta)); no fixed physical-frequency band is asserted.
"""

from fractions import Fraction

import sympy as sp

K_MIN, K_MAX = Fraction(19, 5), Fraction(21, 5)
S_MIN, S_MAX = Fraction(42), Fraction(65)
LEFT, RIGHT = -Fraction(1, 4), Fraction(0)
PULSE_RIGHT = -Fraction(1, 8)
DELTA_MAX = Fraction(1, 625)
X = sp.Symbol("x", real=True)
ETA, TAU, DELTA = sp.symbols("eta tau delta", positive=True)


def exact(value):
    if type(value) is int or isinstance(value, Fraction):
        return Fraction(value)
    raise TypeError("Green bounds require exact int/Fraction inputs")


def volterra_bounds(k_min=K_MIN, k_max=K_MAX, s_max=S_MAX, width=RIGHT-LEFT):
    """First-zero positive-flux constants; validity is an explicit margin.

    For k_min<=k<=k_max, 0<=s<=s_max, h<=width, the bootstrap gives
    k*R_x>=1-s_max*h²/(2*k_min). Its positivity is needed before using
    k_max in the lower flux integral. A longer interval can fail this
    sufficient criterion without implying an actual negative Green kernel.
    """
    k_min, k_max, s_max, width = map(exact, (k_min, k_max, s_max, width))
    if k_min <= 0 or k_max < k_min or s_max < 0 or width < 0:
        raise ValueError("Positive ordered kinetic bounds and nonnegative s/width are required")
    flux = 1-s_max*width**2/(2*k_min)
    derivative = flux/k_max
    sharp = (1-s_max*width**2/(6*k_min))/k_max
    direct = 1/k_max-s_max*width**2/(6*k_min**2)
    return {"k_min": k_min, "k_max": k_max, "s_max": s_max, "width": width,
            "flux_lower": flux, "R_x_lower": derivative,
            "R_over_distance_lower": sharp,
            "direct_Volterra_R_over_distance_lower": direct,
            "positive_flux_bootstrap": flux > 0}


def green_constants():
    values = volterra_bounds()
    values.update({"s_min": S_MIN, "left": LEFT, "right": RIGHT,
                   "R_target": Fraction(3, 16), "R_x_target": Fraction(1, 10),
                   "R_margin": values["R_over_distance_lower"]-Fraction(3, 16),
                   "R_x_margin": values["R_x_lower"]-Fraction(1, 10),
                   "direct_R_margin": values["direct_Volterra_R_over_distance_lower"]-Fraction(3, 16)})
    return values


def response_constants():
    """Amplitude-normalized bounds for the stated positive smooth pulses."""
    g = green_constants()
    width, distance = PULSE_RIGHT-LEFT, -PULSE_RIGHT
    area_min = Fraction(1, 16)
    moment_max = (LEFT**2-PULSE_RIGHT**2)/2
    return {"pulse_left": LEFT, "pulse_right": PULSE_RIGHT,
            "pulse_width": width, "endpoint_distance_min": distance,
            "pulse_area_min_over_eta": area_min,
            "Q_final_lower_over_eta": g["R_target"]*distance*S_MIN*area_min,
            "Q_x_final_lower_over_eta": g["R_x_target"]*S_MIN*area_min,
            "Q_final_sharper_lower_over_eta": g["R_over_distance_lower"]*distance*S_MIN*area_min,
            "Q_x_final_sharper_lower_over_eta": g["R_x_lower"]*S_MIN*area_min,
            "Q_whole_window_upper_over_eta": S_MAX*moment_max/K_MIN,
            "Q_x_whole_window_upper_over_eta": S_MAX*width/K_MIN,
            "physical_Q_time_derivative_lower": sp.Rational(21, 80)*ETA/(TAU*sp.sqrt(DELTA)),
            "physical_window_duration": TAU*sp.sqrt(DELTA)/4,
            "delta_max": DELTA_MAX}


def smooth_step(argument):
    """C-infinity step: zero through 0, one from 1, flat at both ends."""
    argument = sp.sympify(argument)
    middle = sp.exp(-1/argument)/(sp.exp(-1/argument)+sp.exp(-1/(1-argument)))
    return sp.Piecewise((0, argument <= 0), (1, argument >= 1), (middle, True))


def pulse():
    """An explicit nonempty smooth pulse class, with an exact plateau area.

    No transcendental integration is used: q=eta on a plateau of length
    5/64, strictly greater than the required 1/16. The support closure is
    strictly inside (-1/4,-1/8), so both the initial and final germs vanish.
    """
    left, right, ramp = -Fraction(31, 128), -Fraction(17, 128), Fraction(1, 64)
    rational = lambda value: sp.Rational(value.numerator, value.denominator)
    profile = ETA*smooth_step((X-rational(left))/rational(ramp))*smooth_step((rational(right)-X)/rational(ramp))
    return {"x": X, "eta": ETA, "support_left": left, "support_right": right,
            "ramp_width": ramp, "plateau_left": left+ramp, "plateau_right": right-ramp,
            "plateau_width": right-left-2*ramp,
            "area_floor_over_eta": right-left-2*ramp,
            "area_margin_over_eta": right-left-2*ramp-Fraction(1, 16),
            "profile": profile}


def flat_derivative_polynomials(order):
    """D_t^n exp(-1/t)=exp(-1/t)*P_n(1/t), on t>0.

    P_(n+1)=z²(P_n-P_n') proves finite-order polynomial growth. Exponential
    decay dominates every such polynomial, establishing all flat endpoint
    derivatives. Actual integer validation precedes every computation.
    """
    if isinstance(order, bool) or not isinstance(order, (int, sp.Integer)):
        raise TypeError("A derivative order must be an actual integer")
    if order < 0:
        raise ValueError("A derivative order must be nonnegative")
    z = sp.Symbol("z", real=True)
    values = [sp.Poly(1, z)]
    for _ in range(int(order)):
        expression = values[-1].as_expr()
        values.append(sp.Poly(sp.expand(z**2*(expression-sp.diff(expression, z))), z))
    return {"z": z, "polynomials": tuple(values)}


def identities():
    """Exact rational equalities used in both Green lower-bound routes."""
    g, r, p = green_constants(), response_constants(), pulse()
    return {"first_zero_flux": g["flux_lower"]-Fraction(283, 608),
            "derivative_floor": g["R_x_lower"]-Fraction(1415, 12768),
            "positive_flux_integral_floor": g["R_over_distance_lower"]-Fraction(7495, 38304),
            "direct_Volterra_floor": g["direct_Volterra_R_over_distance_lower"]-Fraction(46385, 242592),
            "sharp_kernel_margin": g["R_margin"]-Fraction(313, 38304),
            "derivative_margin": g["R_x_margin"]-Fraction(691, 63840),
            "direct_kernel_margin": g["direct_R_margin"]-Fraction(899, 242592),
            "endpoint_Q_lower": r["Q_final_lower_over_eta"]-Fraction(63, 1024),
            "endpoint_Qx_lower": r["Q_x_final_lower_over_eta"]-Fraction(21, 80),
            "whole_window_Q_upper": r["Q_whole_window_upper_over_eta"]-Fraction(975, 2432),
            "whole_window_Qx_upper": r["Q_x_whole_window_upper_over_eta"]-Fraction(325, 152),
            "explicit_plateau_width": p["plateau_width"]-Fraction(5, 64),
            "explicit_pulse_area_margin": p["area_margin_over_eta"]-Fraction(1, 64),
            "inherited_v_radius": DELTA_MAX*(RIGHT-LEFT)**2-Fraction(1, 10000)}


def checks():
    g, r, p = green_constants(), response_constants(), pulse()
    values = {name: value == 0 for name, value in identities().items()}
    values.update({"first_zero_flux_positive": g["flux_lower"] > 0,
                   "sharp_R_above_target": g["R_margin"] > 0,
                   "R_x_above_target": g["R_x_margin"] > 0,
                   "direct_R_above_target": g["direct_R_margin"] > 0,
                   "pulse_support_inside": LEFT < p["support_left"] < p["support_right"] < PULSE_RIGHT,
                   "pulse_area_strict_margin": p["area_margin_over_eta"] > 0,
                   "positive_history_response": r["Q_final_lower_over_eta"] > 0,
                   "positive_history_velocity": r["Q_x_final_lower_over_eta"] > 0})
    return values


def calibration():
    return {"claim": "P8-S6.34.OWN_RETARDED", "green": green_constants(),
            "response": response_constants(), "pulse": pulse(),
            "normalization": "(k Q_x)_x+s(Q-q)=0; k(a)R_x(a,a)=1; Q=Q_x=0 at x=-1/4",
            "physical_clock": "x=T/(tau*sqrt(delta)), 0<delta<=1/625",
            "scope": "linear own-f response to prescribed off-shell physical metric; no physical matter-source, full-parent, fixed low-frequency EFT or Green-selected action claim"}
