"""All-order local WKB recurrence and bounded Cauchy-frequency correction."""
from functools import cache, lru_cache

import sympy as sp
from p8_vector_clock_matching import continuation
from p8_vector_state import wkb

u, z = wkb.u, wkb.z


def D0(expression):
    return sp.factor(sp.diff(expression, u)+wkb.background()["z_prime"]*sp.diff(expression, z))


@lru_cache(maxsize=None, typed=True)
def coefficient(kind, order):
    if type(order) is not int or order < 0:
        raise ValueError("Require a nonnegative native integer WKB coefficient order")
    U = continuation.coefficients(kind)["U"].subs(continuation.local.dimension, 3)
    if order == 0:
        return sp.Integer(1)
    lam = wkb.background()["lambda"]
    previous = [coefficient(kind, n) for n in range(order)]
    inverse = [sp.Integer(1)]
    for n in range(1, order):
        inverse.append(sp.factor(-sum(previous[j]*inverse[n-j] for j in range(1, n+1))))
    DS = [sp.Integer(0)]+[D0(previous[n])-2*n*lam*previous[n] for n in range(1, order)]
    rate = [lam]+[sp.factor(sum(DS[j]*inverse[n-j] for j in range(1, n+1)))
                  for n in range(1, order)]
    product = -sum(previous[j]*previous[order-j] for j in range(1, order))
    derivative = -D0(rate[-1])/2+(order-1)*lam*rate[-1]
    square = sum(rate[j]*rate[order-1-j] for j in range(order))/4
    residual = product+derivative+square-(U if order == 1 else 0)
    return sp.factor(residual/2)


@lru_cache(maxsize=None, typed=True)
def coefficient_bounds(kind, order):
    P = coefficient(kind, order)
    slope = sp.factor(D0(P)+(1-2*order)*wkb.background()["lambda"]*P)
    values = wkb.box_bound(P), wkb.box_bound(slope)
    return {"coefficient_upper": values[0]["absolute_upper"],
            "frequency_coefficient_slope_upper": values[1]["absolute_upper"],
            "coefficient_reconstruction": values[0]["reconstruction"],
            "slope_reconstruction": values[1]["reconstruction"]}


@cache
def leading_correction_bounds(kind):
    d = coefficient_bounds(kind, 3)
    # Every higher term is assigned at most 2^-n in the same
    # omega^-5 norm, for value and slope; their sum is below one.
    F, G = d["coefficient_upper"]+1, d["frequency_coefficient_slope_upper"]+1
    m = wkb.MASS_TIME_MIN
    return {"sixth_order_coefficient_upper": d["coefficient_upper"],
            "sixth_order_frequency_slope_upper": d["frequency_coefficient_slope_upper"],
            "full_frequency_correction_over_inverse_frequency_fifth_upper": F,
            "full_frequency_slope_correction_over_inverse_frequency_fifth_upper": G,
            "initial_beta_over_inverse_frequency_sixth_upper": 2*F+4*(G+4*F)/m}


@cache
def low_order_checks():
    out = {}
    for kind in ("transverse", "longitudinal"):
        prior = continuation.coefficients(kind)
        for n, name in ((1, "P2"), (2, "P4")):
            out[kind+"_recurrence_replays_frozen_"+name] = sp.factor(coefficient(kind, n)-prior[name].subs(continuation.local.dimension, 3))
        p, r = coefficient(kind, 1), coefficient(kind, 2)
        b2 = D0(p)-2*wkb.background()["lambda"]*p
        B = prior["B"].subs(continuation.local.dimension, 3)
        out[kind+"_sixth_order_cancels_leading_frozen_residual"] = sp.factor(
            2*coefficient(kind, 3)+p**3+4*p*r-B-sp.Rational(3, 4)*b2**2)
        out.update({kind+"_"+name: value for name, value in coefficient_bounds(kind, 3).items()
                    if name.endswith("reconstruction")})
    return out
