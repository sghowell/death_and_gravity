"""Uniform complex-dimension bounds for the subtracted exact mode family."""
from functools import cache

import sympy as sp
from p8_vector_dimensional import local
from p8_vector_state import wkb

u, z = wkb.u, wkb.z
eta = sp.Symbol("dimension_minus_three")
RESIDUAL_CONSTANT = sp.Integer(10_000_000)
TAIL_CONSTANT = sp.Integer(1_000_000)


def derivative(expression):
    H = wkb.background()["H"]
    return sp.factor(sp.diff(expression, u)-2*H*z*(1-z)*sp.diff(expression, z))


@cache
def coefficients(kind):
    H = wkb.background()["H"]
    mapping = {local.z: z}
    mapping.update({symbol: sp.diff(H, u, j) for j, symbol in
                    enumerate((local.H, local.Hd, local.Hdd, local.Hddd, local.Hdddd))})
    data = local.reference(kind)
    U, p, r = [sp.factor(data[name].subs(mapping)) for name in ("U", "P2", "P4")]
    lam = -H*z
    b2, b4 = derivative(p)-2*lam*p, derivative(r)-4*lam*r
    A, B = p**2+2*r, -derivative(b4)/2+sp.Rational(5, 2)*lam*b4
    return {"U": U, "P2": p, "P4": r, "B2": sp.factor(b2), "B4": sp.factor(b4),
            "A": sp.factor(A), "B": sp.factor(B),
            "second_order_cancellation": sp.factor(-U-derivative(lam)/2+lam**2/4-2*p),
            "fourth_order_cancellation": sp.factor(-derivative(b2)/2+sp.Rational(3, 2)*lam*b2-A)}


def box_bound(expression):
    if (not isinstance(expression, sp.Expr) or expression.has(sp.Float)
            or expression.is_finite is False
            or expression.has(sp.oo, -sp.oo, sp.zoo, sp.nan)):
        raise ValueError("Require an exact finite rational expression")
    value = sp.cancel(expression.subs(local.dimension, 3+eta))
    numerator, denominator = sp.fraction(value)
    poly_den = sp.Poly(denominator, u)
    power, factor = poly_den.degree()//2, poly_den.LC()
    if factor.is_positive is not True or sp.expand(denominator-factor*(1+u**2)**power) != 0:
        raise ValueError("Require a positive quadratic-power denominator")
    polynomial = sp.Poly(numerator/factor, u, z, eta)
    if any(not coefficient.is_Rational for _, coefficient in polynomial.terms()):
        raise ValueError("Require exact rational polynomial coefficients")
    bound = sum(abs(coefficient)/sp.Integer(2)**powers[0]/sp.Integer(4)**powers[2]
                for powers, coefficient in polynomial.terms())
    return {"absolute_upper": bound, "denominator_power": power,
            "reconstruction": sp.factor(value-polynomial.as_expr()/(1+u**2)**power)}


@cache
def abstract_residual():
    p, r, lam, lamd, b2, b4, B = sp.symbols("p r lambda lambda_dot b2 b4 B")
    t = sp.Symbol("inverse_frequency_squared")
    S, A = 1+p*t+r*t**2, p**2+2*r
    DS = b2*t+b4*t**2
    DDS = (lam*b2-2*A)*t+(lam*b4-2*B)*t**2
    U = -2*p-lamd/2+lam**2/4
    lograte, lograte_d = lam+DS/S, lamd+DDS/S-DS**2/S**2
    actual = (1-S**2)/t-U-lograte_d/2+lograte**2/4
    expected = t**2*(-A*(p+r*t)/S-2*p*r-r**2*t+B/S+sp.Rational(3, 4)*(b2+t*b4)**2/S**2)
    return {"abstract_exact_WKB_residual_given_checked_cancellations": sp.factor(actual-expected)}


@cache
def bilinear_mode_checks():
    W, Wd, phase, rho = sp.symbols("W W_dot phase rho")
    plus, minus = sp.exp(-sp.I*phase)/sp.sqrt(2*W), sp.exp(sp.I*phase)/sp.sqrt(2*W)
    dp, dm = (-Wd/(2*W)-sp.I*W)*plus, (-Wd/(2*W)+sp.I*W)*minus
    matrix = sp.Matrix([[-sp.I*rho*plus*minus, -sp.I*rho*minus**2],
                        [sp.I*rho*plus**2, sp.I*rho*plus*minus]])
    return {"analytic_pair_Wronskian_without_complex_dimension_conjugation": sp.simplify(plus*dm-minus*dp-sp.I),
            "analytic_pair_variation_constraint": (sp.Matrix([[plus, minus]])*matrix).applyfunc(sp.simplify),
            "analytic_pair_exact_forced_equation": (sp.Matrix([[dp, dm]])*matrix+rho*sp.Matrix([[plus, minus]])).applyfunc(sp.simplify)}


@cache
def bounds(kind):
    data = coefficients(kind)
    b = {name: box_bound(data[name])["absolute_upper"] for name in ("P2", "P4", "B2", "B4", "A", "B")}
    P, Q, B2, B4 = [b[name] for name in ("P2", "P4", "B2", "B4")]
    m, tmax = wkb.MASS_TIME_MIN, 1/wkb.MASS_TIME_MIN**2
    deviation = P*tmax+Q*tmax**2
    phase = P/m+Q/m**3
    lograte = 2+(B2*tmax+B4*tmax**2)/(1-deviation)
    residual = 2*b["A"]*(P+Q*tmax)+2*P*Q+Q**2*tmax+2*b["B"]+3*(B2+B4*tmax)**2
    R1 = 2*(2*P*Q+P**3+Q*(Q+P**2)*tmax)
    R2 = 2*(P**2+Q+P*Q*tmax)
    R3 = 4*(B4+B2*(2*P+(P**2+2*Q)*tmax+2*P*Q*tmax**2+Q**2*tmax**3))
    tail = sp.Rational(3, 2)*R1+sp.Rational(5, 4)*(sp.Rational(81, 16)*R2+sp.Rational(9, 4)*R3+2*(B2+B4*tmax)**2)
    return {"coefficient_envelopes": b, "S_deviation_upper": deviation, "phase_imaginary_part_upper": phase,
            "reference_log_rate_upper": lograte, "residual_constant_upper": residual,
            "residual_integer_upper": sp.ceiling(residual), "reference_tail_upper": tail,
            "reference_tail_integer_upper": sp.ceiling(tail)}


@cache
def proof_checks():
    out = {}
    for kind in ("transverse", "longitudinal"):
        d = bounds(kind)
        out[kind+"_no_reference_zero_on_complex_dimension_disk"] = bool(d["S_deviation_upper"] < sp.Rational(1, 2))
        out[kind+"_phase_bound_gives_mode_multiplier_below_two"] = bool(d["phase_imaginary_part_upper"] < sp.Rational(1, 2))
        out[kind+"_uniform_log_rate_bound"] = bool(d["reference_log_rate_upper"] < 4)
        out[kind+"_uniform_complex_residual_constant"] = bool(d["residual_integer_upper"] <= RESIDUAL_CONSTANT)
        out[kind+"_uniform_complex_reference_tail_constant"] = bool(d["reference_tail_integer_upper"] <= TAIL_CONSTANT)
    out.update({"dimension_disk_excludes_zero": bool(3-sp.Rational(1, 4) > 0),
                "complex_pressure_kinetic_weight_bound": bool(sp.Rational(13, 11) < sp.Rational(5, 4)),
                "complex_pressure_potential_weight_bound": bool(sp.Rational(5, 11) < sp.Rational(3, 2)),
                "uniform_complex_variation_exponent_below_quarter": bool(5*RESIDUAL_CONSTANT/wkb.MASS_TIME_MIN**5 < sp.Rational(1, 4)),
                "reference_physical_derivative_below_two_frequencies": bool(sp.Rational(3, 2)+sp.Rational(21, 4)/wkb.MASS_TIME_MIN < 2),
                "complex_reference_product_energy_constant": (sp.Rational(5, 4)*16+sp.Rational(3, 2)*4)/2 == 13,
                "complex_total_evolution_envelope_constant": bool(sp.Rational(13, 4)*52*5 < 1000),
                "complex_total_reference_tail_constant": bool(sp.Rational(13, 4)/4 < 1),
                "high_momentum_evolution_envelope_integrable": bool(sp.Rational(13, 4)-1-4 < -1),
                "high_momentum_reference_tail_integrable": bool(sp.Rational(13, 4)-1-5 < -1),
                "low_momentum_envelope_integrable": bool(sp.Rational(11, 4)-1 > -1)})
    return out


@cache
def checks():
    out = {**abstract_residual(), **bilinear_mode_checks()}
    for kind in ("transverse", "longitudinal"):
        data = coefficients(kind)
        for key in ("second_order_cancellation", "fourth_order_cancellation"):
            out[kind+"_complex_dimension_"+key] = data[key]
        for key in ("P2", "P4", "B2", "B4", "A", "B"):
            out[kind+"_complex_dimension_box_"+key] = box_bound(data[key])["reconstruction"]
    return out
