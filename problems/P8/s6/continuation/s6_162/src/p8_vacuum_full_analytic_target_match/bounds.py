"""Exact classical coefficient majorants and an even Cauchy remainder."""

from fractions import Fraction
from functools import cache

import sympy as s

from . import jets


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, s.Rational)):
        raise TypeError("Require a finite exact rational")
    return s.Rational(value)


def weighted_majorant(polynomial, caps=(1, 4, 64, 64, 256)):
    if not isinstance(caps, tuple) or len(caps) != 5:
        raise TypeError("Require five ordered exact jet-invariant caps")
    caps = tuple(map(exact, caps))
    if min(caps) < 0:
        raise ValueError("Invariant caps must be nonnegative")
    p = s.Poly(polynomial, *jets.VARIABLES)
    total = s.Integer(0)
    for powers, coefficient in p.terms():
        c = exact(coefficient)
        total += abs(c) * s.prod(v**n for v, n in zip(caps, powers))
    return total


def even_tail(circle_upper, radius, last_degree=8):
    M, rho = map(exact, (circle_upper, radius))
    if M < 0 or rho <= 1:
        raise ValueError(
            "Require a nonnegative circle bound and radius greater than one"
        )
    if type(last_degree) is not int or last_degree < 0 or last_degree % 2:
        raise TypeError("Require a nonnegative native even last degree")
    return M / (rho ** (last_degree + 2) * (1 - rho**-2))


def integrated_error(sixth, eighth, tail, jet_count=21):
    values = tuple(map(exact, (sixth, eighth, tail)))
    if min(values) < 0:
        raise ValueError("Density remainder allowances must be nonnegative")
    if type(jet_count) is not int or jet_count < 1:
        raise TypeError("Require a positive native jet-component count")
    return jet_count * sum(values)


def compose(old_error, target_error):
    old_error, target_error = map(exact, (old_error, target_error))
    if min(old_error, target_error) < 0:
        raise ValueError("Both common-class action allowances must be nonnegative")
    return old_error + target_error


@cache
def data():
    rho, M, W, e6, e8, e10 = s.symbols("rho M W E6 E8 E10", positive=True)
    j = s.Symbol("j", integer=True, nonnegative=True)
    geometric = M / (rho**10 * (1 - rho**-2))
    return {
        "ordered_jet_component_count": 21,
        "invariant_caps": dict(
            zip(("Phi", "Y", "L3", "L4", "L5"), (1, 4, 64, 64, 256))
        ),
        "even_Cauchy_tail_after_degree_eight": geometric,
        "pointwise_remainder": e6 * W**6 + e8 * W**8 + e10 * W**10,
        "integrated_coefficient": 21 * (e6 + e8 + e10),
        "analytic_to_integrated_argument": "For W=max of the 21 canonical jet moduli, rescale each jet by W. At W=0 the density remainder vanishes. For 0<W<=1 the homogeneous sixth/eighth pieces and the even analytic tail are bounded by E6 W^6+E8 W^8+E10 W^10 <=(E6+E8+E10)W^2. Parseval on the Euclidean Fourier unit ball gives integral W^2 <= sum of the 21 jet L2 norms squared <=21 ||Psi||_2^2. No infinite-volume constant is integrated.",
        "checks": {
            "ordered_scalar_first_second_components": 1 + 4 + 16 - 21,
            "Lorentz_gradient_absolute_sum": 4 - 4,
            "L3_Box_times_Z_absolute_sum": 4 * 16 - 64,
            "L4_three_contracted_indices_absolute_sum": 4**3 - 64,
            "L5_squared_Z_absolute_sum": 16**2 - 256,
            "even_geometric_tail": s.simplify(geometric * (1 - rho**-2) - M / rho**10),
            "first_unretained_even_degree": 8 + 2 - 10,
            "homogeneous_rescaling_tail_power": s.simplify(
                (W / rho) ** (10 + 2 * j) - W**10 * rho**-10 * (W * W / rho**2) ** j
            ),
        },
        "scope": "These are classical field-amplitude and common-function-class action bounds, not a physical momentum cutoff, quantum omitted-order estimate or bound on metric variations.",
    }
