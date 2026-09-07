"""A fixed physical-time pulse that provably loads the own-f tensor.

The source is a prescribed metric tensor q, not canonical matter stress.
The estimates use the actual S6.33/S6.35 kinetic, pump, and source maps.
No coefficient limit at the singular center is needed for loading.
"""

from fractions import Fraction
from functools import cache

import sympy as sp

L = Fraction(1, 100)
A = L / 2
DELTA_MAX = Fraction(1, 10**6)
REMAINDER_CAP = Fraction(44)
PUMP_CAP = Fraction(22)
COFRAME_MIN = Fraction(19, 10)
COFRAME_MAX = Fraction(21, 10)
U = sp.Symbol("u", real=True)
ETA = sp.Symbol("eta", positive=True)


def rational(value):
    if isinstance(value, bool) or value is sp.true or value is sp.false:
        raise TypeError("exact rational data exclude booleans")
    if type(value) is int or isinstance(value, Fraction):
        return Fraction(value)
    if isinstance(value, sp.Rational):
        return Fraction(int(value.p), int(value.q))
    raise TypeError("exact int/Fraction/SymPy Rational data are required")


def source_bounds():
    """Bounds derived from V=16/(delta+8u^2)+R and mass=V+pump.

    Valid on [-L,-L/2], 0<delta<=1e-6, for the source-pinned branch.
    k=b^3/N>1 makes sqrt(k)*mass >= 1/L^2 a strict coefficient bound.
    Multiplication by q>=0 is weak at points where q vanishes.
    """
    pole_min = 16 / (DELTA_MAX + 8 * L**2)
    return {
        "left": -L, "load_endpoint": -A,
        "delta_max": DELTA_MAX,
        "kinetic_lower": COFRAME_MIN**3 / COFRAME_MAX,
        "kinetic_upper": COFRAME_MAX**3 / COFRAME_MIN,
        "potential_lower": pole_min - REMAINDER_CAP,
        "potential_upper": 8 / L**2 + REMAINDER_CAP,
        "mass_lower": pole_min - REMAINDER_CAP - PUMP_CAP,
        "canonical_source_lower": 1 / L**2,
        "punctured_D_lower": 2 - 2 / (1 + A**2)**4,
        "remainder_cap": REMAINDER_CAP,
        "pump_cap": PUMP_CAP,
    }


def volterra_bounds(potential_max, lag):
    """Sufficient first-zero bound for G''+V G=0, G(a)=0, G'(a)=1.

    The assumption is 0<=V<=potential_max.  A failed positive-derivative
    margin only invalidates this sufficient argument, not the actual ODE.
    """
    potential_max, lag = rational(potential_max), rational(lag)
    if potential_max < 0 or lag < 0:
        raise ValueError("potential ceiling and lag must be nonnegative")
    derivative = 1 - potential_max * lag**2 / 2
    ratio = 1 - potential_max * lag**2 / 6
    return {
        "potential_max": potential_max, "lag": lag,
        "derivative_lower": derivative,
        "G_over_lag_lower": ratio,
        "positive_derivative_bootstrap": derivative > 0,
    }


def smooth_step(argument):
    argument = sp.sympify(argument)
    middle = sp.exp(-1 / argument) / (sp.exp(-1 / argument) + sp.exp(-1 / (1 - argument)))
    return sp.Piecewise((0, argument <= 0), (1, argument >= 1), (middle, True))


@cache
def pulse():
    """One delta-independent C-infinity pulse, with support strictly inside.

    Its plateau alone has area 5*eta*L/64, greater than eta*L/16.
    Hence existence does not rely on a numerical bump integral.
    """
    left, right, ramp = -95 * L / 128, -81 * L / 128, L / 64
    convert = lambda value: sp.Rational(value.numerator, value.denominator)
    profile = ETA * smooth_step((U - convert(left)) / convert(ramp)) * smooth_step(
        (convert(right) - U) / convert(ramp)
    )
    return {
        "u": U, "eta": ETA, "q": profile,
        "support_left": left, "support_right": right,
        "ramp": ramp,
        "plateau_left": left + ramp,
        "plateau_right": right - ramp,
        "plateau_width": right - left - 2 * ramp,
        "required_area_over_eta": L / 16,
        "source_to_load_lag_min": -A - right,
        "source_to_load_lag_max": -A - left,
    }


def flat_derivative_polynomial(order):
    if isinstance(order, bool) or not isinstance(order, (int, sp.Integer)):
        raise TypeError("derivative order must be an actual integer")
    order = int(order)
    if order < 0:
        raise ValueError("derivative order must be nonnegative")
    return _flat_derivative_polynomial(order)


@cache
def _flat_derivative_polynomial(order):
    z = sp.Symbol("z", real=True)
    polynomial = sp.Integer(1)
    for _ in range(order):
        polynomial = sp.expand(z**2 * (polynomial - sp.diff(polynomial, z)))
    return polynomial


def calibration():
    bounds = source_bounds()
    green = volterra_bounds(bounds["potential_upper"], L / 4)
    y_floor = Fraction(9, 10) * (L / 8) * bounds["canonical_source_lower"] * (L / 16)
    yp_floor = Fraction(2, 3) * bounds["canonical_source_lower"] * (L / 16)
    return {
        **bounds,
        "L": L, "a": A,
        "green_lag_max": L / 4,
        "green_derivative_lower": green["derivative_lower"],
        "green_ratio_lower": green["G_over_lag_lower"],
        "green_derivative_gate": Fraction(2, 3),
        "green_ratio_gate": Fraction(9, 10),
        "Y_load_lower_over_eta": y_floor,
        "Y_u_load_lower_over_eta": yp_floor,
        "Y_load_sharper_lower_over_eta": green["G_over_lag_lower"] / 128,
        "Y_u_load_sharper_lower_over_eta": green["derivative_lower"] / (16 * L),
        "canonical_source": "sqrt(k)*(A_profile/D)*q; k=b^3/N",
        "initial_data": "Q(-L)=Q_u(-L)=0, equivalently Y(-L)=Y_u(-L)=0",
        "loading_limit": "ordinary coefficient convergence on the fixed punctured interval only",
        "physical_clock": "u=T/tau; Q_u=tau*Q_T; tau is fixed as delta tends to zero",
        "matter_source_response": False,
    }


@cache
def identities():
    k = sp.Function("k", positive=True)(U)
    y, q = sp.Function("Y")(U), sp.Function("q")(U)
    mass = sp.Function("mass")(U)
    root = sp.sqrt(k)
    old = sp.diff(k * sp.diff(y / root, U), U) + k * mass * (y / root - q)
    new = sp.diff(y, U, 2) + (mass - sp.diff(root, U, 2) / root) * y - root * mass * q
    h, ceiling = sp.symbols("h Vmax", nonnegative=True)
    s = sp.symbols("s", nonnegative=True)
    v, delta = sp.symbols("v delta", nonnegative=True)
    D = 2 + delta - 2 / (1 + v)**4
    return {
        "literal_canonical_source": sp.simplify(old / root - new),
        "positive_flux_integral": sp.integrate(ceiling * s, (s, 0, h)) - ceiling * h**2 / 2,
        "positive_G_integral": sp.integrate(ceiling * s**2 / 2, (s, 0, h)) - ceiling * h**3 / 6,
        "G_bound_derivative": sp.diff(h - ceiling * h**3 / 6, h) - (1 - ceiling * h**2 / 2),
        "punctured_D_monotonicity": sp.diff(D, v) - 8 / (1 + v)**5,
        "pulse_is_delta_independent": sp.diff(pulse()["q"], delta),
        "zero_initial_whitening_value": root * 0,
        "zero_initial_whitening_derivative": sp.diff(root, U) * 0 + root * 0,
    }


def checks():
    data, p = calibration(), pulse()
    return {
        "coframes_positive": COFRAME_MIN > 0 and COFRAME_MAX > COFRAME_MIN,
        "kinetic_exceeds_one": data["kinetic_lower"] > 1,
        "kinetic_below_five": data["kinetic_upper"] < 5,
        "positive_potential": data["potential_lower"] > 0,
        "mass_exceeds_source_gate": data["mass_lower"] > data["canonical_source_lower"],
        "punctured_denominator_positive": data["punctured_D_lower"] > 0,
        "potential_upper_value": data["potential_upper"] == 80044,
        "positive_derivative_bootstrap": data["green_derivative_lower"] > 0,
        "green_derivative_strict": data["green_derivative_lower"] > data["green_derivative_gate"],
        "green_ratio_strict": data["green_ratio_lower"] > data["green_ratio_gate"],
        "green_derivative_value": data["green_derivative_lower"] == Fraction(59989, 80000),
        "green_ratio_value": data["green_ratio_lower"] == Fraction(219989, 240000),
        "support_left_strict": p["support_left"] > -3 * L / 4,
        "support_right_strict": p["support_right"] < -5 * L / 8,
        "source_finishes_before_loading": p["support_right"] < -A,
        "plateau_area_strict": p["plateau_width"] > p["required_area_over_eta"],
        "plateau_value": p["plateau_width"] == 5 * L / 64,
        "lag_above_L_over_eight": p["source_to_load_lag_min"] > L / 8,
        "lag_below_L_over_four": p["source_to_load_lag_max"] < L / 4,
        "positive_Y_loading": data["Y_load_lower_over_eta"] == Fraction(9, 1280),
        "positive_Y_u_loading": data["Y_u_load_lower_over_eta"] == Fraction(1, 24) / L,
        "strict_loading_survives_limit": data["Y_load_sharper_lower_over_eta"] > data["Y_load_lower_over_eta"],
    }
