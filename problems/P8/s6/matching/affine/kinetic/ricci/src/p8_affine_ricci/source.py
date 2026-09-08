"""Full rolling stationary connection and its Ricci-difference source."""
from functools import cache
from itertools import product

import sympy as sp
from p8_affine import connection as old

from . import geometry

u = sp.Symbol("u", real=True)
H, h, hp, n, nd = sp.symbols("Hubble h h_prime n n_dot", real=True)
DN = sp.symbols("dn1:4", real=True)


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.factor))


@cache
def background():
    a = (1+u**2)**2
    return {"a": a, "h": (1+u**2)**3, "H": sp.factor(sp.diff(a, u)/a)}


@cache
def stationary_variation():
    """Differentiate all 64 generic components, including every coefficient jet."""
    p, x = sp.symbols("p x", real=True)
    px = -1/(8*h*p)
    cx = 2*(1-4*p)*px/x-2*(p-2*p**2)/x**2
    cxx = sp.factor((sp.diff(cx, x)+px*sp.diff(cx, p)).subs({p: sp.Rational(1, 2), x: -1}))
    clock = {old.P: sp.Rational(1, 2), old.S: 1, old.PP: 0,
             old.PX: -1/(4*h), old.CP: 0, old.CX: -1/(2*h), old.F3: 0}
    coefficient_jets = {old.P: -n/(2*h), old.S: -n, old.PP: hp*n/(2*h**2),
                        old.PX: -n/(4*h**2), old.CP: hp*n/h**2,
                        old.CX: 2*n*cxx, old.F3: hp*n/(2*h**2)}
    hessian_background = {old.H[i, j]: (-H if i == j and i > 0 else 0)
                          for i in range(4) for j in range(i, 4)}
    hessian_jets = {old.H[0, 0]: -nd}
    hessian_jets.update({old.H[0, i]: -DN[i-1] for i in range(1, 4)})
    hessian_jets.update({old.H[i, j]: sp.Symbol(f"delta_H{i}{j}", real=True)
                         for i in range(1, 4) for j in range(i, 4)})
    stationary = old.eliminate()["solution"]
    actual = []
    for value in stationary:
        coefficient = sum(sp.diff(value, key).subs(clock).subs(hessian_background)*jet
                          for key, jet in coefficient_jets.items())
        coefficient += sum(sp.diff(value, key).subs(clock).subs(hessian_background)*jet
                           for key, jet in hessian_jets.items())
        actual.append(sp.factor(coefficient.subs(hp, 3*H*h/2)))
    expected = {index: sp.S.Zero for index in old.INDICES}
    expected[0, 0, 0] = -15*H*n/(8*h)
    for i in range(1, 4):
        expected[0, i, 0] = expected[i, 0, 0] = -3*DN[i-1]/(2*h)
        expected[0, i, i] = expected[i, 0, i] = (H*n-nd)/(2*h)
        expected[i, i, 0] = 5*H*n/(8*h)
    for i, j, k in product(range(1, 4), repeat=3):
        expected[i, j, k] = (int(j == k)*DN[i-1]-int(i == k)*DN[j-1])/(2*h)
    return {"actual": sp.ImmutableMatrix(actual),
            "expected": sp.ImmutableMatrix([expected[index] for index in old.INDICES]),
            "background": _matrix(stationary.subs(clock).subs(hessian_background)),
            "coefficient_jets": coefficient_jets, "hessian_jets": hessian_jets,
            "c_XX_clock": cxx}


def ricci_difference(kappa, christoffel, metric, inverse_metric, coordinates):
    """Literal linear curvature in the displayed diagonal FLRW coordinates."""
    @cache
    def derivative(a, b, c, d):
        return sp.diff(kappa[a, b, c], coordinates[d])+sum(
            christoffel[a, e, d]*kappa[e, b, c]
            -christoffel[e, b, d]*kappa[a, e, c]
            -christoffel[e, c, d]*kappa[a, b, e] for e in range(4))

    def curvature(a, b, c, d):
        return derivative(a, b, d, c)-derivative(a, b, c, d)

    def difference(mu, nu):
        return sum(curvature(a, nu, a, mu)
                   -metric[nu]*inverse_metric[a]*curvature(nu, a, a, mu) for a in range(4))

    return sp.ImmutableMatrix([sp.factor((difference(mu, nu)-difference(nu, mu))/2)
                               for mu, nu in geometry.PAIRS])


@cache
def rolling():
    coordinates = (u, *sp.symbols("x y z", real=True))
    lapse = sp.Function("n")(*coordinates)
    bg = background()
    a, h_actual, hubble = bg["a"], bg["h"], bg["H"]
    metric = (-sp.S.One, a**2, a**2, a**2)
    inverse = tuple(1/value for value in metric)
    christoffel = {index: sp.S.Zero for index in old.INDICES}
    for i in range(1, 4):
        christoffel[0, i, i] = a**2*hubble
        christoffel[i, 0, i] = christoffel[i, i, 0] = hubble
    # The background stationary distortion is zero, so first-order frame
    # changes multiply zero. Only the background coframe factors remain.
    mapping = {n: lapse, nd: sp.diff(lapse, u), H: hubble, h: h_actual}
    mapping.update({DN[i-1]: sp.diff(lapse, coordinates[i])/a for i in range(1, 4)})
    ortho = stationary_variation()["actual"]
    coordinate = {index: sp.factor(ortho[position].subs(mapping)
                                   *a**(int(index[1] > 0)+int(index[2] > 0)-int(index[0] > 0)))
                  for position, index in enumerate(old.INDICES)}
    actual = ricci_difference(coordinate, christoffel, metric, inverse, coordinates)
    coefficient = -5*hubble/(2*h_actual)
    expected = sp.Matrix([coefficient*sp.diff(lapse, coordinates[i]) for i in range(1, 4)]+[0, 0, 0])
    return {"coordinates": coordinates, "lapse": lapse, "metric": metric,
            "inverse_metric": inverse, "Christoffel": christoffel,
            "coordinate_stationary": coordinate, "Cstar": actual,
            "expected": expected, "coefficient": sp.factor(coefficient)}


@cache
def checks():
    data, roll, bg = stationary_variation(), rolling(), background()
    return {"all_64_background_components": data["background"],
            "all_64_full_first_variations": _matrix(data["actual"]-data["expected"]),
            "full_curved_Cstar": _matrix(roll["Cstar"]-roll["expected"]),
            "Cstar_coefficient": sp.factor(roll["coefficient"]+10*u/(1+u**2)**4),
            "actual_h_derivative": sp.factor(sp.diff(bg["h"], u)-3*bg["H"]*bg["h"]/2),
            "full_coefficient_c_XX": sp.factor(data["c_XX_clock"]+(4*h-1)/(4*h**2))}
