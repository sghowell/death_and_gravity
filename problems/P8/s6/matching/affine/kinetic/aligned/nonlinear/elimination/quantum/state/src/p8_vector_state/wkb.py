"""Fourth-order WKB algebra and continuous rational coefficient bounds."""
from functools import cache

import sympy as sp
from p8_affine_aligned import modes

u = modes.u
z = sp.Symbol("physical_momentum_fraction", nonnegative=True)
t = sp.Symbol("inverse_frequency_squared", positive=True)
MASS_TIME_MIN = sp.Integer(1000)


@cache
def background():
    H = modes.canonical()["H"]
    return {"H": H, "Hd": sp.diff(H, u), "z_prime": -2*H*z*(1-z), "lambda": -H*z}


def D(expression):
    b = background()
    return sp.factor(sp.diff(expression, u)+b["z_prime"]*sp.diff(expression, z)
                     -2*b["lambda"]*t*sp.diff(expression, t))


@cache
def frequency(kind):
    b = background()
    H, Hd, lam = b["H"], b["Hd"], b["lambda"]
    if kind == "transverse":
        U = Hd/2+H**2/4
    elif kind == "longitudinal":
        U = (sp.Rational(1, 2)+z)*Hd+(sp.Rational(1, 4)-z+3*z**2)*H**2
    else:
        raise ValueError("Require transverse or longitudinal mode")
    P2 = sp.factor(-U/2-D(lam)/4+lam**2/8)
    P4 = sp.factor(-P2**2/2-D(D(P2))/4+sp.Rational(5, 4)*lam*D(P2)
                   +(D(lam)/2-sp.Rational(3, 2)*lam**2)*P2)
    B2, B4 = sp.factor(D(P2)-2*lam*P2), sp.factor(D(P4)-4*lam*P4)
    C4 = sp.factor(D(B4)-4*lam*B4)
    A, B = sp.factor(P2**2+2*P4), sp.factor(-C4/2+lam*B4/2)
    S = 1+P2*t+P4*t**2
    rate = lam+D(S)/S
    actual_residual = (1-S**2)/t-U-D(rate)/2+rate**2/4
    expected_residual = t**2*(-A*(P2+P4*t)/S-2*P2*P4-P4**2*t+B/S
                                   +sp.Rational(3, 4)*(B2+t*B4)**2/S**2)
    return {"U": U, "P2": P2, "P4": P4, "B2": B2, "B4": B4, "A": A, "B": B,
            "S": S, "reference_log_rate": rate, "residual": expected_residual,
            "exact_fourth_order_WKB_residual": sp.factor(actual_residual-expected_residual),
            "WKB_second_order_cancellation": sp.factor(-U-D(lam)/2+lam**2/4-2*P2),
            "WKB_fourth_order_cancellation": sp.factor(-D(B2)/2+sp.Rational(3, 2)*lam*B2-A)}


def box_bound(expression):
    """Exact L1 polynomial bound on |u|<=1/2 and 0<=z<=1."""
    numerator, denominator = sp.fraction(sp.cancel(expression))
    poly_den = sp.Poly(denominator, u)
    n, factor = poly_den.degree()//2, poly_den.LC()
    if factor.is_positive is not True or sp.expand(denominator-factor*(1+u**2)**n) != 0:
        raise ValueError("Require a positive quadratic-power denominator")
    poly = sp.Poly(numerator/factor, u, z)
    if any(not coefficient.is_Rational for _, coefficient in poly.terms()):
        raise ValueError("Require exact fixed rational coefficients")
    bound = sum(abs(coefficient)/sp.Integer(2)**power[0] for power, coefficient in poly.terms())
    return {"denominator_power": n, "absolute_upper": bound,
            "reconstruction": sp.factor(expression-poly.as_expr()/(1+u**2)**n)}


@cache
def bounds(kind):
    d = frequency(kind)
    b = {name: box_bound(d[name])["absolute_upper"] for name in ("P2", "P4", "B2", "B4", "A", "B")}
    tmax = 1/MASS_TIME_MIN**2
    deviation = b["P2"]*tmax+b["P4"]*tmax**2
    log_rate = box_bound(background()["lambda"])["absolute_upper"]+(b["B2"]*tmax+b["B4"]*tmax**2)/(1-deviation)
    constant = (2*b["A"]*(b["P2"]+b["P4"]*tmax)+2*b["P2"]*b["P4"]+b["P4"]**2*tmax
                +2*b["B"]+3*(b["B2"]+b["B4"]*tmax)**2)
    return {"coefficient_envelopes": b, "S_deviation_upper": deviation,
            "reference_log_rate_absolute_upper": log_rate,
            "residual_over_inverse_frequency_fourth_upper": constant,
            "residual_constant_integer_upper": sp.ceiling(constant)}


@cache
def checks():
    out = {}
    for kind in ("transverse", "longitudinal"):
        d = frequency(kind)
        out.update({kind+"_"+name: d[name] for name in
                    ("exact_fourth_order_WKB_residual", "WKB_second_order_cancellation", "WKB_fourth_order_cancellation")})
        out.update({kind+"_coefficient_box_reconstruction_"+name: box_bound(d[name])["reconstruction"]
                    for name in ("P2", "P4", "B2", "B4", "A", "B")})
    return out


@cache
def physical_chain_rule_checks():
    bg = modes.canonical()
    H = bg["H"]
    q, mass2 = sp.symbols("physical_q mass_squared", positive=True)
    actual_z, actual_t = q/(q+mass2), 1/(q+mass2)
    out = {"fixed_comoving_momentum_fraction_derivative": sp.factor(sp.diff(actual_z, q)*(-2*H*q)
                                        -background()["z_prime"].subs(z, actual_z)),
           "fixed_comoving_inverse_frequency_derivative": sp.factor(sp.diff(actual_t, q)*(-2*H*q)
                                        +2*background()["lambda"].subs(z, actual_z)*actual_t)}
    for kind in ("transverse", "longitudinal"):
        out[kind+"_original_canonical_frequency_correction"] = sp.factor(frequency(kind)["U"].subs(z, bg["r"]/(1+bg["r"]))
                                                            -bg[kind]["correction"])
    return out
