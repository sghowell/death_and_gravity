"""Independent coordinate-curvature and diagonal-ADM action derivation."""
from functools import cache

import sympy as sp

TIME = sp.Symbol("t", real=True)
R = sp.Function("r")(TIME)
Q = sp.Function("q")(TIME)
A, AP, RV, L, L1, L2 = sp.symbols("A A_prime r_value ell ell_prime ell_second", real=True)
Q1, Q2, Q3 = sp.symbols("q1 q2 q3", real=True)
A2, A3, A4, L3 = sp.symbols("A_second A_third A_fourth ell_third", real=True)
QJ = sp.symbols("qj0:7", real=True)


@cache
def coordinate_schouten():
    """Derive Ricci from all 4D coordinate Christoffels, not supplied jets."""
    g = sp.diag(1, -sp.exp(Q), -sp.exp(-Q), -1)
    metric = R**2*g
    inverse = metric.inv()

    def partial(value, index):
        return sp.diff(value, TIME) if index == 0 else sp.S.Zero

    gamma = [[[sp.simplify(sum(inverse[a, b]*(partial(metric[b, k], j)
                + partial(metric[b, j], k)-partial(metric[j, k], b))
                for b in range(4))/2) for k in range(4)] for j in range(4)] for a in range(4)]
    ricci = sp.zeros(4)
    for i in range(4):
        for j in range(4):
            ricci[i, j] = sp.simplify(sum(
                partial(gamma[a][i][j], a)-partial(gamma[a][i][a], j)
                + sum(gamma[a][a][b]*gamma[b][i][j]-gamma[a][j][b]*gamma[b][i][a]
                      for b in range(4)) for a in range(4)))
    scalar = sp.simplify(sp.trace(inverse*ricci))
    correction = g.inv()*(ricci-metric*scalar/6)
    substitutions = {sp.diff(R, TIME, 2): R*(L1+L**2), sp.diff(R, TIME): R*L,
                     sp.diff(Q, TIME): Q1, sp.diff(Q, TIME, 2): Q2}
    mixed = correction.applyfunc(lambda value: sp.simplify(value.subs(substitutions)))
    if any(mixed[i, j] != 0 for i in range(4) for j in range(4) if i != j):
        raise ValueError("The literal diagonal probe developed a mixed off-diagonal component")
    return {"metric": metric, "ricci": ricci, "scalar": scalar,
            "mixed_correction_over_A": mixed,
            "h": tuple(A*mixed[i, i] for i in range(4))}


def _time_derivative(expression):
    if sp.sympify(expression).has(A4, L3, Q3, QJ[-1]):
        raise ValueError("The independent finite jet map needs another derivative")
    mapping = {A: AP, AP: A2, A2: A3, A3: A4,
               RV: RV*L, L: L1, L1: L2, L2: L3, Q1: Q2, Q2: Q3}
    mapping.update({QJ[index]: QJ[index+1] for index in range(6)})
    return sp.expand(sum(sp.diff(expression, key)*value for key, value in mapping.items()))


def _quadratic(expression):
    polynomial = sp.Poly(sp.expand(expression), Q1, Q2, Q3)
    return sp.Add(*(coefficient*Q1**powers[0]*Q2**powers[1]*Q3**powers[2]
                    for powers, coefficient in polynomial.terms() if sum(powers) == 2))


@cache
def diagonal_action():
    """Coefficient epsilon^2 of exact diagonal ADM EH plus potential cubic.

    f_i=f0_i*(1+epsilon*h_i). The series for volume/lapse is obtained
    from its logarithm; the logarithmic scale velocities are differentiated
    before extracting epsilon^2. This independently retains h'_0 where it
    enters the volume/lapse prefactor, without giving the lapse a new mode.
    All densities here are divided by the hidden squared Planck mass.
    """
    h = coordinate_schouten()["h"]
    hd = tuple(_time_derivative(value) for value in h)
    log_factor_1 = (sum(h[1:])-h[0])/2
    log_factor_2 = (-sum(value**2 for value in h[1:])+h[0]**2)/4
    prefactor = (sp.S.One, log_factor_1, log_factor_2+log_factor_1**2/2)
    velocities = [(base, hd[i]/2, -h[i]*hd[i]/2)
                  for i, base in enumerate((L+Q1/2, L-Q1/2, L), start=1)]
    pair_series = [sum(velocities[i][j]*velocities[k][degree-j]
                       for i in range(3) for k in range(i+1, 3) for j in range(degree+1))
                   for degree in range(3)]
    eh = -RV**2*sum(prefactor[j]*pair_series[2-j] for j in range(3))
    # The beta2=beta3=0 own-root cubic follows independently from the
    # literal diagonal sqrt and sqrt(det) series (Newton identities).
    epsilon = sp.Symbol("epsilon", real=True)
    sqrt_series = [1+epsilon*x/2-epsilon**2*x*x/8+epsilon**3*x**3/16 for x in h]
    determinant_series = sp.Poly(1, epsilon)
    for root in sqrt_series:
        product = determinant_series*sp.Poly(root, epsilon)
        determinant_series = sp.Poly(sum(product.nth(k)*epsilon**k for k in range(4)), epsilon)
    trace_cubic = sum(value**3 for value in h)/16
    det_cubic = determinant_series.nth(3)
    potential = -2*RV**2*(trace_cubic-det_cubic)/A
    raw = _quadratic(eh+potential)
    raw_coefficients = {f"a{i}{j}": sp.factor(sp.Poly(raw, Q1, Q2, Q3).coeff_monomial(left*right))
                        for i, left in enumerate((Q1, Q2, Q3), start=1)
                        for j, right in enumerate((Q1, Q2, Q3), start=1) if j >= i}
    return {"h": h, "hd": hd, "prefactor": prefactor,
            "EH_quadratic": sp.factor(_quadratic(eh)),
            "potential_quadratic": sp.factor(_quadratic(potential)),
            "raw_quadratic": sp.factor(raw), "raw_coefficients": raw_coefficients}


@cache
def normal_and_euler():
    raw = diagonal_action()["raw_coefficients"]
    B3 = raw["a33"]
    B2 = sp.factor(raw["a22"]-raw["a13"]-_time_derivative(raw["a23"])/2)
    B1 = sp.factor(raw["a11"]+_time_derivative(_time_derivative(raw["a13"]))/2
                   -_time_derivative(raw["a12"])/2)
    normal = B1*QJ[1]**2+B2*QJ[2]**2+B3*QJ[3]**2
    euler = sp.S.Zero
    for order in (1, 2, 3):
        term = sp.diff(normal, QJ[order])
        for _ in range(order):
            term = _time_derivative(term)
        euler += (-1)**order*term
    return {"B1": B1, "B2": B2, "B3": B3,
            "euler": sp.expand(euler),
            "euler_coefficients": {n: sp.factor(sp.expand(euler).coeff(QJ[n])) for n in range(7)}}


@cache
def independent_center():
    """Power-series root of a different rational profile representation.

    This uses finite binomial/quotient series in v=u^2, not the agent's
    chain-rule jet formulas. The degree-two v series suffices for the
    derivatives actually present in the center Euler coefficients.
    """
    c = sp.Symbol("c", positive=True)
    v = sp.Symbol("v", real=True)
    d = 1+v
    W = c*(1-7*v)*d**4+12*v
    r3 = 8*c*(1-v)/(d**8*W)
    x = sp.series(r3/8-1, v, 0, 3).removeO()
    r = sp.series(2*(1+x/3-x*x/9), v, 0, 3).removeO().expand()
    abar = sp.series(r*c*d**10*(c*d**4-2)/(32*(1-v)), v, 0, 3).removeO().expand()
    ell_over_u = sp.series(2*sp.diff(r, v)/r, v, 0, 2).removeO().expand()
    substitutions = {
        RV: 2, A: abar.coeff(v, 0), AP: 0, A2: 2*abar.coeff(v, 1),
        A3: 0, A4: 24*abar.coeff(v, 2),
        L: 0, L1: ell_over_u.coeff(v, 0), L2: 0, L3: 6*ell_over_u.coeff(v, 1)}
    values = normal_and_euler()
    return {"c": c, "v": v, "r_cubed": r3, "r_series": r,
            "A_bar_series": abar, "substitution": substitutions,
            "B": {n: sp.factor(values[f"B{n}"].subs(substitutions)) for n in (1, 2, 3)},
            "E": {n: sp.factor(value.subs(substitutions))
                  for n, value in values["euler_coefficients"].items()}}
