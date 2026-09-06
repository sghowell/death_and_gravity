"""Derived finite actual-stress and SEE-residual bounds, with domain guards.

The numerical theorem keeps the SAME prepared A.7 metric/state, named
lambda=2sqrt(2)A eta_star^2 and gamma=0. It is not a nearby-solution theorem.
"""

from functools import cache

import sympy as sp
from p8a_remainder import geometry as prior_geometry
from p8a_remainder import mode_bounds as prior_modes
from p8a_remainder import stress as prior_stress
from p8a_remainder.mode_bounds import exact_nonnegative

from . import kernel, plateau

NAMES = ("density", "pressure", "EED")
ROUNDED_FIRST = dict(zip(NAMES, (14000, 236000, 360000), strict=True))
ROUNDED_SECOND = dict(zip(NAMES, (8*10**14, 9*10**14, 17*10**14), strict=True))
FROZEN_NUMERATORS = dict(zip(NAMES, (3, 5, 9), strict=True))


@cache
def calibration():
    data = prior_geometry.calibration()
    b = list(map(sp.Rational, data["potential_bounds"]))
    cap = sp.Rational(data["delta_bar"])
    delta = sp.Symbol("delta", nonnegative=True)
    modes = prior_modes.polynomial_bounds([delta*v for v in b], sp.Integer(3), sp.Integer(1),
                                          sp.Integer(2), sp.Integer(2))["integrated_remainder_derivatives"]
    m = [sp.cancel(value/delta**2).subs(delta, cap) for value in modes]
    remainder = prior_stress.calibration_bounds(m, b, cap)
    k = kernel.calibration()
    s = [k[j]/2+b[j]/10 for j in range(3)]
    derivative_piece = prior_stress.calibration_bounds(s, b, cap)
    euler = plateau.euler_error_constants()
    first = {name: derivative_piece[name]+euler[name] for name in NAMES}
    curvature_squared = {"density": 3*b[0]**2/128, "pressure": 7*b[0]**2/768, "EED": 13*b[0]**2/512}
    second = {name: remainder[name]+curvature_squared[name] for name in NAMES}
    if any(first[name] >= ROUNDED_FIRST[name] or second[name] >= ROUNDED_SECOND[name] for name in NAMES):
        raise ValueError("The claimed rounded actual-stress coefficients are not proved")
    return {"delta_bar": cap, "potential_bounds": b, "mode_remainder_constants": m,
            "kernel_norm_constants": k, "effective_S_Born_norm_constants": s,
            "Euler_first_constants": euler, "derivative_first_constants": derivative_piece,
            "R_squared_history_second_constants": curvature_squared,
            "A7_stress_remainder_second_constants": remainder,
            "first_constants": first, "second_constants": second,
            "rounded_first_constants": ROUNDED_FIRST.copy(), "rounded_second_constants": ROUNDED_SECOND.copy()}


def amplitude(epsilon, d, *, normalization=1, eta_star=1):
    """d is the resolved kappa*hbar/(46080*pi^2), as in the immutable A.5."""
    epsilon = exact_nonnegative(epsilon)
    d = exact_nonnegative(d, positive=True)
    normalization = exact_nonnegative(normalization, positive=True)
    eta_star = exact_nonnegative(eta_star, positive=True)
    return prior_geometry.check_amplitude(16*epsilon*d/(normalization**2*eta_star**4))


def coefficients(*, rounded=True):
    if not isinstance(rounded, bool):
        raise TypeError("The exact/rounded coefficient selector must be a boolean")
    data = calibration()
    prefix = "rounded_" if rounded else ""
    return data[prefix+"first_constants"], data[prefix+"second_constants"]


def stress_error(delta, *, normalization=1, eta_star=1, hbar=1, rounded=True):
    """Actual minus positive A.3 reference at the SAME proper-time label y.

    Only the already named lambda/gamma and derived preparation are covered.
    These bounds need not bound the actual stress relative to itself.
    """
    delta = prior_geometry.check_amplitude(delta)
    normalization = exact_nonnegative(normalization, positive=True)
    eta_star = exact_nonnegative(eta_star, positive=True)
    hbar = exact_nonnegative(hbar)
    first, second = coefficients(rounded=rounded)
    pref = hbar/(sp.pi**2*normalization**4*eta_star**8)
    return {name: pref*(delta*first[name]+delta**2*second[name]) for name in NAMES}


def see_residual(delta, *, normalization=1, eta_star=1, rounded=True):
    """Actual physical Einstein+radiation residual, not only a frozen defect.

    Requires epsilon*kappa*hbar/pi^2=2880*delta*A^2*eta_star^4, zero Lambda,
    no extra explicit curvature sources, and the target 2<=y<=3. The actual
    prepared metric is not asserted to solve SEE or be close to a solution.
    """
    delta = prior_geometry.check_amplitude(delta)
    normalization = exact_nonnegative(normalization, positive=True)
    eta_star = exact_nonnegative(eta_star, positive=True)
    first, second = coefficients(rounded=rounded)
    pref = delta**2/(normalization**2*eta_star**4)
    return {name: pref*(sp.Rational(FROZEN_NUMERATORS[name], 1922)
                       + 2880*(first[name]+delta*second[name])) for name in NAMES}


def reference_accuracy(delta, *, rounded=True):
    """Uniform error/reference ratios; references are positive pinned A.3 data."""
    delta = prior_geometry.check_amplitude(delta)
    first, second = coefficients(rounded=rounded)
    denominators = {"density": 960, "pressure": 576, "EED": 320}
    return {name: denominators[name]*3**8*(delta*first[name]+delta**2*second[name]) for name in NAMES}
