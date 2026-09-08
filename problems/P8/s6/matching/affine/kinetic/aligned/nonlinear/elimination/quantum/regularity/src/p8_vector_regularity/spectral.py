"""Higher-time-derivative readouts and exact projected subtraction algebra."""
from functools import lru_cache

import sympy as sp
from p8_vector_hadamard import series
from p8_vector_state import wkb
from p8_vector_subtraction import tail

u, z = wkb.u, wkb.z
omega2 = sp.Symbol("frequency_squared", positive=True)


def clean(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.factor))


def D0(value):
    return sp.diff(value, u)+wkb.background()["z_prime"]*sp.diff(value, z)


def D(value):
    return D0(value)+2*wkb.background()["lambda"]*omega2*sp.diff(value, omega2)


@lru_cache(maxsize=None, typed=True)
def rows(name, derivative_order):
    if name not in tail.physical_weights() or type(derivative_order) is not int or not 0 <= derivative_order <= 5:
        raise ValueError("Require an actual physical readout and native derivative order 0..5")
    weights = tail.physical_weights()[name]
    H, lam = wkb.background()["H"], wkb.background()["lambda"]
    d = weights["c1"]-lam/2
    matrix = sp.Matrix([[2*d, 2, 0], [-omega2, 0, 1], [0, -2*omega2, -2*d]])
    row = sp.Matrix([[weights["B"]*omega2, 0, weights["A"]]])
    for _ in range(derivative_order):
        row = clean(row.applyfunc(D)-3*H*row+row*matrix)
    return clean(row)


@lru_cache(maxsize=None, typed=True)
def row_bounds(name, derivative_order):
    row = rows(name, derivative_order)
    coefficients, majorants = [], []
    reconstructions = []
    degree_controls = []
    for column, component in enumerate(row):
        poly = sp.Poly(component, omega2)
        exact, upper = {}, {}
        for (power,), value in poly.terms():
            bound = wkb.box_bound(value)
            exact[-power] = value
            upper[-power] = bound["absolute_upper"]
            reconstructions.append(bound["reconstruction"])
            if value != 0:
                degree_controls.append(2*power+(-1, 0, 1)[column] <= derivative_order+1)
        coefficients.append(exact)
        majorants.append(upper)
    readout = sum(sum(poly.values())*(1, 2, 3)[column] for column, poly in enumerate(majorants))
    return {"coefficients_in_t": coefficients, "coefficient_majorants": majorants,
            "normalized_reference_product_envelope": readout,
            "all_frequency_powers_at_most_derivative_order_plus_one": all(degree_controls),
            "box_reconstructions": reconstructions}


@lru_cache(maxsize=None, typed=True)
def adiabatic_rows(name, derivative_order):
    rows(name, derivative_order)
    terms = tail.reference_terms()[name]
    lam, H = wkb.background()["lambda"], wkb.background()["H"]
    out = {}
    for order, key in enumerate(("zero", "second", "fourth")):
        value = terms[key]
        for _ in range(derivative_order):
            value = sp.factor(D0(value)+((1-2*order)*lam-3*H)*value)
        out[order] = value
    return out


def add(*polynomials):
    out = {}
    for polynomial in polynomials:
        for degree, value in polynomial.items():
            out[degree] = out.get(degree, sp.Integer(0))+value
    return out


def multiply(left, right):
    out = {}
    for j, a in left.items():
        for k, b in right.items():
            out[j+k] = out.get(j+k, sp.Integer(0))+a*b
    return out


def scale(polynomial, factor, shift=0):
    return {degree+shift: factor*value for degree, value in polynomial.items()}


def numerator(S, R, row, c, adiabatic, absolute=False):
    X, Y, Z = row
    S2, S3 = multiply(S, S), multiply(multiply(S, S), S)
    S4 = multiply(S3, S)
    first = add(X, scale(Y, c if absolute else -c), scale(Z, c**2))
    second = add(Y if absolute else scale(Y, -1), scale(Z, 2*c))
    terms = add(scale(multiply(first, S2), 4, 1),
                scale(multiply(multiply(second, R), S), 2, 1),
                scale(multiply(Z, multiply(R, R)), 1, 1),
                scale(multiply(Z, S4), 4),
                scale(multiply(adiabatic, S3), 4 if absolute else -4))
    # For derivative order <=5, the most negative power is t^-2.
    return scale(terms, 1, 2)


@lru_cache(maxsize=None, typed=True)
def reference_tail(name, derivative_order, wkb_order=4):
    if type(wkb_order) is not int or not 2 <= wkb_order <= 4:
        raise ValueError("Require native WKB coefficient order 2..4")
    data = row_bounds(name, derivative_order)
    weights = tail.physical_weights()[name]
    kind = weights["mode"]
    lam = wkb.background()["lambda"]
    S = {0: sp.Integer(1)}
    R = {}
    Sb, Rb = {0: sp.Integer(1)}, {}
    for order in range(1, wkb_order+1):
        P = series.coefficient(kind, order)
        B = sp.factor(D0(P)-2*order*lam*P)
        S[order], R[order] = P, B
        Sb[order], Rb[order] = wkb.box_bound(P)["absolute_upper"], wkb.box_bound(B)["absolute_upper"]
    adiabatic = adiabatic_rows(name, derivative_order)
    polynomial = numerator(S, R, data["coefficients_in_t"], weights["c1"], adiabatic)
    low = {order: sp.factor(polynomial.get(order, 0)) for order in range(5)}
    upper = numerator(Sb, Rb, data["coefficient_majorants"], sp.Integer(2),
                      {order: wkb.box_bound(value)["absolute_upper"] for order, value in adiabatic.items()},
                      absolute=True)
    mass = wkb.MASS_TIME_MIN
    bound = 2*sum(value/mass**(2*(order-5)) for order, value in upper.items() if order >= 5)
    return {"low_coefficient_residuals": low,
            "integrable_tail_certified": all(value == 0 for value in low.values()),
            "reference_bracket_tail_over_t_cubed_upper": bound,
            "reference_bracket_tail_integer_upper": sp.ceiling(bound),
            "no_lower_Laurent_power": min(polynomial) >= 0,
            "all_majorant_coefficients_nonnegative": all(value >= 0 for value in upper.values())}


def algebra_checks():
    vr, vi, pr, pi = sp.symbols("vr vi pr pi", real=True)
    d, omega2_generic = sp.symbols("d omega_squared", real=True)
    bilinears = sp.Matrix([vr**2+vi**2, vr*pr+vi*pi, pr**2+pi**2])
    exact_derivative = bilinears.jacobian((vr, vi, pr, pi))*sp.Matrix([
        pr+d*vr, pi+d*vi, -omega2_generic*vr-d*pr, -omega2_generic*vi-d*pi])
    generator = sp.Matrix([[2*d, 2, 0], [-omega2_generic, 0, 1], [0, -2*omega2_generic, -2*d]])
    S, R, c1, t, X, Y, Z, ad = sp.symbols("S R c1 t X Y Z ad", nonzero=True)
    c = c1+R/(2*S)
    projected = (t*(X-c*Y+c**2*Z)+S**2*Z)/S-ad
    polynomial = numerator({0: S}, {0: R}, [{0: X}, {0: Y}, {0: Z}], c1, {0: ad})
    actual_numerator = sum(value*t**degree for degree, value in polynomial.items())
    return {"real_canonical_bilinear_generator": (exact_derivative-generator*bilinears).applyfunc(sp.expand),
            "projected_readout_common_denominator_identity": sp.factor(4*t**2*S**3*projected-actual_numerator)}


def low_reference_failure():
    data = reference_tail("longitudinal_energy", 2, 2)
    return sp.factor(data["low_coefficient_residuals"][4].subs({u: 0, z: 1}))
