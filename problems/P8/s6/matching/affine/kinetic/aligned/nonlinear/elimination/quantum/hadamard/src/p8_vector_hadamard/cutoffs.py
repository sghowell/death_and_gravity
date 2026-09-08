"""Explicit locally finite all-order cutoffs and initial Gaussian matching."""
from functools import cache, lru_cache

import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_vector_state import wkb

from . import series

INITIAL_MIXING_CONSTANT = sp.Integer(10_000_000)


@lru_cache(maxsize=None, typed=True)
def threshold(kind, order, mass_time_product):
    m = exact(mass_time_product, "m0_tau")
    if type(order) is not int or order < 3 or not bool(m >= wkb.MASS_TIME_MIN):
        raise ValueError("Require order>=3 and exact m0*tau>=1000")
    if order == 3:
        series.coefficient(kind, 0)
        return 2*m
    data = series.coefficient_bounds(kind, order)
    Cn = max(data["coefficient_upper"], data["frequency_coefficient_slope_upper"])
    return max(2*threshold(kind, order-1, m), 2**order*(1+Cn))


def turn_on(s):
    if (not isinstance(s, sp.Expr) or s.is_real is not True
            or s.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan)):
        raise ValueError("Require an exact finite cutoff argument")
    bump_left, bump_right = sp.exp(-1/(s-1)), sp.exp(-1/(2-s))
    return sp.Piecewise((0, s <= 1), (1, s >= 2), (bump_left/(bump_left+bump_right), True))


def bump_derivative_polynomial(order):
    if type(order) is not int or order < 0:
        raise ValueError("Require a nonnegative native integer derivative order")
    y = sp.Symbol("inverse_positive_argument", positive=True)
    poly = sp.Integer(1)
    for _ in range(order):
        poly = sp.expand(y**2*(poly-sp.diff(poly, y)))
    return y, poly


@cache
def initial_data_identities():
    W, V = sp.symbols("old_frequency new_frequency", positive=True)
    J, K = sp.symbols("old_half_log_rate new_half_log_rate", real=True)
    f, v = 1/sp.sqrt(2*W), 1/sp.sqrt(2*V)
    fd, vd = (-J-sp.I*W)*f, (-K-sp.I*V)*v
    A = (W+V-sp.I*(K-J))/(2*sp.sqrt(W*V))
    B = (W-V+sp.I*(K-J))/(2*sp.sqrt(W*V))
    return {"Cauchy_value_matches_initial_Bogoliubov_coefficients": sp.simplify(A*f+B*f-v),
            "Cauchy_derivative_matches_initial_Bogoliubov_coefficients": sp.simplify(A*fd+B*sp.conjugate(fd)-vd),
            "initial_state_change_preserves_CCR": sp.simplify(A*sp.conjugate(A)-B*sp.conjugate(B)-1)}


@cache
def proof_checks():
    m = wkb.MASS_TIME_MIN
    out = {}
    for kind in ("transverse", "longitudinal"):
        d = series.leading_correction_bounds(kind)
        F = d["full_frequency_correction_over_inverse_frequency_fifth_upper"]
        G = d["full_frequency_slope_correction_over_inverse_frequency_fifth_upper"]
        out[kind+"_positive_all_order_Cauchy_frequency"] = bool(F/m**6 < sp.Rational(1, 4))
        out[kind+"_initial_state_mixing_constant"] = bool(d["initial_beta_over_inverse_frequency_sixth_upper"] <= INITIAL_MIXING_CONSTANT)
        out[kind+"_complex_initial_state_data_envelope"] = bool(4*F+8*(G+4*F)/m < INITIAL_MIXING_CONSTANT)
    out["initial_mixing_below_one"] = bool(INITIAL_MIXING_CONSTANT/m**6 < 1)
    out["higher_coefficient_geometric_tail_below_one"] = sp.Rational(1, 16)/(1-sp.Rational(1, 2)) < 1
    return {name: bool(value) for name, value in out.items()}


def checks():
    out = dict(initial_data_identities())
    s = sp.Symbol("positive_argument", positive=True)
    for n in range(6):
        y, polynomial = bump_derivative_polynomial(n)
        out["bump_derivative_polynomial_"+str(n)] = sp.simplify(
            sp.diff(sp.exp(-1/s), s, n)-sp.exp(-1/s)*polynomial.subs(y, 1/s))
    return out


def structural_checks():
    t, new = sp.symbols("t new_coefficient")
    sigma = sp.Symbol("relative_root", nonzero=True)
    out = {"complex_diagonal_amplitude_has_quadratic_root_remainder":
           sp.factor((1+sigma**2)/(2*sigma)-1-(sigma-1)**2/(2*sigma))}
    for n in range(1, 7):
        previous = sp.symbols("P1:"+str(n))
        S = 1+sum(value*t**j for j, value in enumerate(previous, 1))
        change = sp.expand(-((S+new*t**n)**2-S**2)/t)
        out["new_WKB_coefficient_first_residual_order_"+str(n)] = change.coeff(t, n-1)+2*new
    return out
