"""Literal four-dimensional jet check of off-shell quartic field redefinition."""

from functools import cache
from itertools import product

import sympy as sp

DIM = 4
ORDER = 6
SIGNS = (1, -1, -1, -1)
INDICES = tuple(a for a in product(range(ORDER + 1), repeat=DIM) if sum(a) <= ORDER)
JETS = {a: sp.Symbol("jet_" + "".join(map(str, a)), real=True) for a in INDICES}
BY_SYMBOL = {v: k for k, v in JETS.items()}
PHI = JETS[(0, 0, 0, 0)]


def derivative(expr, axis):
    if type(axis) is not int or not 0 <= axis < DIM:
        raise ValueError("A native four-dimensional coordinate index is required")
    return _derivative(sp.sympify(expr), axis)


@cache
def _derivative(expr, axis):
    terms = []
    for s in expr.free_symbols.intersection(BY_SYMBOL):
        a = list(BY_SYMBOL[s])
        a[axis] += 1
        if tuple(a) not in JETS:
            raise ValueError("The exact sixth-derivative jet domain was exceeded")
        terms.append(sp.diff(expr, s) * JETS[tuple(a)])
    return sp.expand(sum(terms, sp.Integer(0)))


def box(expr):
    return sp.expand(
        sum(SIGNS[i] * derivative(derivative(expr, i), i) for i in range(DIM))
    )


def K(expr):
    return sp.expand(box(expr) + 2 * expr)


def gradient(expr):
    return sp.Matrix([SIGNS[i] * derivative(expr, i) for i in range(DIM)])


def divergence(v):
    return sp.expand(sum(derivative(v[i], i) for i in range(DIM)))


@cache
def data():
    phi = PHI
    lam, gamma, c = sp.symbols("quartic_lambda quartic_gamma redundant_c", real=True)
    v = sp.Matrix([derivative(phi, i) for i in range(DIM)])
    H = sp.Matrix([[derivative(v[i], j) for j in range(DIM)] for i in range(DIM)])
    X = sp.expand(sum(SIGNS[i] * v[i] ** 2 for i in range(DIM)))
    raised = sp.diag(*SIGNS) * v
    Z = sp.expand((raised.T * H * raised)[0])
    L3 = sp.expand(box(phi) * Z)
    L4 = sp.expand((raised.T * H * sp.diag(*SIGNS) * H * raised)[0])
    E = box(phi) + phi
    P = sp.expand(phi * E)
    J = phi**2
    kj, k2j = K(J), K(K(J))
    k3j = K(k2j)
    target = lam * X**2 - gamma * phi**4 / 3 - 2 * gamma * (L3 - L4)
    low = c * J * J / 12 - c * J * kj / 8
    truncated = low + lam * J * k2j / 4 - gamma * J * k3j / 8
    R = (
        -c * phi**3 / 6
        + lam * (2 * phi * X + phi**2 * E)
        + gamma * (2 * Z + phi * X + phi**3 / 3)
        - gamma * phi * K(X)
        - gamma * phi * K(P) / 2
    )
    current = (
        -c * phi**3 * gradient(phi) / 12
        + lam * (J * gradient(kj) - kj * gradient(J)) / 4
        - gamma
        * (
            J * gradient(k2j)
            - k2j * gradient(J)
            + 4 * (X * gradient(P) - P * gradient(X))
        )
        / 8
        - gamma * X * gradient(X) / 2
        - gamma * (phi * X + phi**3 / 3) * gradient(phi)
    )
    split3 = (
        2 * gamma * L4
        - gamma * X**2
        - (
            gamma * E * (2 * Z + phi * X + phi**3 / 3)
            - 2 * gamma * (L3 - L4)
            - gamma * phi**4 / 3
        )
        - divergence(-gamma * (phi * X + phi**3 / 3) * gradient(phi))
    )
    return {
        "dimension": DIM,
        "maximum_jet_derivative": ORDER,
        "independent_jet_count": len(JETS),
        "signature": SIGNS,
        "phi": phi,
        "free_equation": E,
        "X": X,
        "Z": Z,
        "L3": L3,
        "L4": L4,
        "target_quartic": target,
        "literal_centered_resolvent_truncation": truncated,
        "cubic_field_redefinition": sp.expand(R),
        "full_boundary_current": current.applyfunc(sp.expand),
        "field_map": "old Phi = new Psi + R(Psi); tangent identity, quartic action equivalence only",
        "checks": {
            "actual_product_K_on_phi_squared": sp.expand(kj - 2 * (X + P)),
            "gradient_X_squared_is_four_L4": sp.expand(
                (gradient(X).T * sp.diag(*SIGNS) * gradient(X))[0] - 4 * L4
            ),
            "low_two_resolvent_terms_are_EOM_plus_boundary": sp.expand(
                low + c * phi**3 * E / 6 - divergence(-c * phi**3 * gradient(phi) / 12)
            ),
            "quadratic_centered_operator_green_identity": sp.expand(
                J * k2j - kj * kj - divergence(J * gradient(kj) - kj * gradient(J))
            ),
            "cubic_centered_operator_green_identity": sp.expand(
                J * k3j - kj * k2j - divergence(J * gradient(k2j) - k2j * gradient(J))
            ),
            "cross_green_identity": sp.expand(
                (X + P) * K(X + P)
                - X * K(X)
                - 2 * P * K(X)
                - P * K(P)
                - divergence(X * gradient(P) - P * gradient(X))
            ),
            "massive_gamma_EOM_and_boundary": sp.expand(split3),
            "full_pointwise_quartic_action_identity": sp.expand(
                truncated - target - E * R - divergence(current)
            ),
            "free_action_cubic_variation": sp.expand(
                (gradient(phi).T * sp.diag(*SIGNS) * gradient(R))[0]
                - phi * R
                + E * R
                - divergence(R * gradient(phi))
            ),
        },
    }
