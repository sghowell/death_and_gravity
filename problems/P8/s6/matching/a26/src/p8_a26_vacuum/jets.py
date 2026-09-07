"""Full symmetric-Hessian Euler and independent-of-boundary bulk controls."""

import math
from functools import cache
from itertools import combinations_with_replacement

import sympy as sp

from .coefficients import exact


@cache
def flat_jet():
    grad = sp.Matrix(sp.symbols("p q s", real=True))
    hs = sp.symbols("H00 H11 H22 H01 H02 H12", real=True)
    h = sp.Matrix([[hs[0], hs[3], hs[4]], [hs[3], hs[1], hs[5]],
                   [hs[4], hs[5], hs[2]]])
    eta = sp.diag(-1, 1, 1)
    raised = eta*grad
    x = (grad.T*eta*grad)[0]
    l3 = (raised.T*h*raised)[0]*sp.trace(eta*h)
    l4 = (raised.T*h*eta*h*raised)[0]
    lag = (l3-l4)/x
    return grad, hs, h, eta, x, lag


@cache
def diagonal_euler():
    """Differentiate ALL H entries first; only then evaluate diagonal H.

    An independent symmetric off-diagonal entry already differentiates both
    slots of the matrix. Adding another factor of two would be incorrect.
    The profiles are quadratic locally, so third and fourth jets vanish.
    """
    grad, hs, _, _, x, lag = flat_jet()
    ds = sp.symbols("d0 d1 d2", real=True)
    at = dict(zip(hs, (*ds, 0, 0, 0), strict=True))
    result = -sum(ds[i]*sp.diff(lag, grad[i], 2).subs(at) for i in range(3))
    result += sum(ds[i]**2*sp.diff(lag, hs[i], grad[i], 2).subs(at) for i in range(3))
    result += sum(ds[i]*ds[j]*sp.diff(lag, hs[k], grad[i], grad[j]).subs(at)
                  for i, j, k in ((0, 1, 3), (0, 2, 4), (1, 2, 5)))
    return grad, ds, x, sp.factor(result)


def evaluate_euler(gradient, diagonal):
    if len(gradient) != 3 or len(diagonal) != 3:
        raise ValueError("Three covariant gradient and diagonal-Hessian entries are required")
    gradient, diagonal = tuple(map(exact, gradient)), tuple(map(exact, diagonal))
    x = -gradient[0]**2+gradient[1]**2+gradient[2]**2
    if x.is_negative is not True:
        raise ValueError("This replay uses the common strictly timelike directional cone")
    grad, ds, _, result = diagonal_euler()
    return sp.simplify(result.subs(dict(zip(grad, gradient, strict=True))
                                   | dict(zip(ds, diagonal, strict=True))))


@cache
def ibp_checks():
    grad, hs, h, eta, x, lag = flat_jet()
    triples = list(combinations_with_replacement(range(3), 3))
    third = dict(zip(triples, sp.symbols(f"T0:{len(triples)}"), strict=True))
    pairs = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))

    def total_derivative(expr, i):
        return sum(h[j, i]*sp.diff(expr, grad[j]) for j in range(3))+sum(
            third[tuple(sorted((*pair, i)))]*sp.diff(expr, hh)
            for pair, hh in zip(pairs, hs, strict=True))

    current = eta*grad*sp.trace(eta*h)-eta*h*eta*grad
    div = sum(total_derivative(current[i], i) for i in range(3))
    d = sp.trace(eta*h)**2-sp.trace(eta*h*eta*h)
    dxj = sum(total_derivative(x, i)*current[i] for i in range(3))
    return {"full_third_jet_divergence": sp.expand(div-d),
            "gradient_of_X_current": sp.factor(dxj-2*x*lag)}


def _average(expr, variables):
    """Exact normalized torus average of a sin/cos polynomial."""
    out = 0
    for powers, coefficient in sp.Poly(sp.expand(expr), *variables).terms():
        if any(power % 2 for power in powers):
            continue
        moment = sp.S.One
        for c, s in ((powers[0]//2, powers[1]//2), (powers[2]//2, powers[3]//2)):
            moment *= sp.Rational(math.factorial(2*c)*math.factorial(2*s),
                                  4**(c+s)*math.factorial(c)*math.factorial(s)*math.factorial(c+s))
        out += coefficient*moment
    return sp.expand(out)


@cache
def fourier_bulk():
    """Expand the literal rational L on psi=t+l*T(t)*(cos x+cos y)."""
    lam, t, tp, tpp = sp.symbols("lambda T Tp Tpp", real=True)
    cx, sx, cy, sy = sp.symbols("cx sx cy sy", real=True)
    g = cx+cy
    grad = sp.Matrix([1+lam*tp*g, -lam*t*sx, -lam*t*sy])
    h = sp.Matrix([[lam*tpp*g, -lam*tp*sx, -lam*tp*sy],
                   [-lam*tp*sx, -lam*t*cx, 0], [-lam*tp*sy, 0, -lam*t*cy]])
    eta = sp.diag(-1, 1, 1)
    v = eta*grad
    numerator = sp.expand((v.T*h*v)[0]*sp.trace(eta*h)-(v.T*h*eta*h*v)[0])
    # Exact inverse coefficients through lambda^2 suffice: numerator starts at 2.
    inv0, inv1 = -1, 2*tp*g
    inv2 = -3*tp**2*g**2-t**2*(sx**2+sy**2)
    coeffs = {2: numerator.coeff(lam, 2)*inv0,
              3: numerator.coeff(lam, 3)*inv0+numerator.coeff(lam, 2)*inv1,
              4: numerator.coeff(lam, 4)*inv0+numerator.coeff(lam, 3)*inv1
                 +numerator.coeff(lam, 2)*inv2}
    average = {order: _average(value, (cx, sx, cy, sy)) for order, value in coeffs.items()}
    return {"T": t, "Tp": tp, "Tpp": tpp, "averages": average,
            "quartic_after_time_IBP": -sp.Rational(1, 2),
            "fourth_difference_multiplier": sp.Integer(24)}


@cache
def checks():
    _, ds, x, euler = diagonal_euler()
    bulk = fourier_bulk()
    t, tp, tpp = bulk["T"], bulk["Tp"], bulk["Tpp"]
    return {"generic_diagonal_euler": sp.factor(x*euler+2*sp.prod(ds)),
            "unit_diagonal": sp.factor(euler.subs(dict.fromkeys(ds, 1))+2/x),
            "additivity_defect": evaluate_euler((3, 0, 0), (2, 2, 2))
            -evaluate_euler((1, 0, 0), (1, 1, 1))-evaluate_euler((2, 0, 0), (1, 1, 1))
            +sp.Rational(13, 18),
            "rank_one_homogeneous_control": evaluate_euler((1, 0, 0), (1, 0, 0)),
            "bulk_quadratic_is_a_time_boundary": sp.factor(bulk["averages"][2]-t*tpp-tp**2),
            "quadratic_boundary_chain_rule": sp.diff(t*tp, t)*tp+sp.diff(t*tp, tp)*tpp-t*tpp-tp**2,
            "bulk_cubic_spatial_average": sp.factor(bulk["averages"][3]),
            "bulk_quartic_raw": sp.expand(bulk["averages"][4]-(t*tpp/2+tp**2)*t**2),
            "bulk_quartic_time_IBP": -sp.Rational(3, 2)+1-bulk["quartic_after_time_IBP"],
            **ibp_checks()}
