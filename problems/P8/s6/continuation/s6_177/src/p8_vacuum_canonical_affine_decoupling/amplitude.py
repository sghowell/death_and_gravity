"""Literal full four-scalar tree amplitude of the fixed nongravitational limit."""

from functools import cache
from itertools import permutations

import sympy as s

from .family import GAMMA, LAMBDA


def pair_matrix(S, T, U, m):
    a, b, c = (S - 2 * m) / 2, (T - 2 * m) / 2, (U - 2 * m) / 2
    return s.Matrix([[m, a, b, c], [a, m, c, b], [b, c, m, a], [c, b, a, m]])


def labelled_vertices(dots, lam, gam):
    """Sum all 24 labelled assignments; i^6=-1 retained in L3 and L4."""
    if dots.shape != (4, 4) or dots != dots.T:
        raise ValueError("Require a symmetric four-leg dot-product matrix")
    perms = list(permutations(range(4)))
    lower = lam * sum(dots[i, j] * dots[k, l] for i, j, k, l in perms)
    L3 = -sum(dots[i, i] * dots[j, k] * dots[k, l] for i, j, k, l in perms)
    L4 = -sum(dots[i, j] * dots[j, k] * dots[k, l] for i, j, k, l in perms)
    gal = gam * sum(
        dots[i, j] * (dots[k, l] ** 2 - dots[k, k] * dots[l, l]) for i, j, k, l in perms
    )
    return lower, 2 * gam * (L4 - L3), gal, -8 * gam


@cache
def data():
    S, T, U, v = s.symbols("s t u crossing_v", real=True)
    m = s.Symbol("mass_squared", positive=True)
    lam, gam = s.symbols("lambda gamma", real=True)
    dots = pair_matrix(S, T, U, m)
    lower, original, gal, potential = labelled_vertices(dots, lam, gam)
    shell = {U: 4 * m - S - T}
    expected = 2 * lam * sum((x - 2 * m) ** 2 for x in (S, T, U))
    A = (lower + original + potential).subs(m, 1)
    shell_one = {U: 4 - S - T}
    formula = expected.subs(m, 1) + 3 * gam * S * T * U - 8 * gam
    forward = s.expand(formula.subs({S: 2 + v, T: 0, U: 2 - v}, simultaneous=True))
    checks = {
        "literal_derivative_contact_pairing": s.expand(lower - expected),
        "original_DHOST_quartic_to_boundary_equivalent_Galileon": s.factor(
            (original - gal).subs(shell)
        ),
        "massive_on_shell_Galileon_stu": s.factor(
            (gal - 3 * gam * S * T * U).subs(shell)
        ),
        "complete_mass_one_tree_amplitude": s.factor((A - formula).subs(shell_one)),
        "crossing_s_t": s.expand(formula - formula.xreplace({S: T, T: S})),
        "crossing_s_u": s.expand(formula - formula.xreplace({S: U, U: S})),
        "potential_four_leg_factorial": potential + s.factorial(4) * gam / 3,
        "forward_full_tree": s.expand(forward - (4 * lam * v * v + 8 * lam - 8 * gam)),
        "forward_b2": s.diff(forward, v, 2) / 2 - 4 * lam,
        "nonzero_gamma_not_erased": s.diff(formula, gam) - (3 * S * T * U - 8),
        "actual_gamma_is_positive": s.cancel(GAMMA - s.Rational(1024, 10**800)),
        "actual_lambda_is_fixed": s.cancel(LAMBDA - s.Rational(1, 10**600)),
    }
    return {
        "canonical_mass_squared": s.S.One,
        "four_field_action": "lambda Y^2-gamma Phi^4/3+2gamma(L4-L3)",
        "literal_labelled_DHOST_vertex": original,
        "complete_on_shell_tree_amplitude": formula,
        "forward_amplitude": forward,
        "actual_forward_b2": 4 * LAMBDA,
        "actual_fixed_gamma": GAMMA,
        "tree_ownership": "No scalar cubic vertex exists. Higher scalar vertices start at six fields; gravity-dependent Ia and the source contact start at eight, and the vector-source vertex has four scalar legs plus one vector. Thus none changes the tree four-scalar amplitude of the decoupled vacuum. This is not all-order matching.",
        "boundary_identity": "div[Y(grad(Phi) Box(Phi)-Hess(Phi).grad(Phi))]=2(L3-L4)+Y[(Box(Phi))^2-Hess(Phi)^2] on flat space; compact/Schwartz variations justify its integrated use only for the constant quartic coefficient.",
        "scope": "The complete higher-field independent functions remain in the limiting action. A positive four-point tree b2 and this exact amplitude do not bound quantum remainders, prove a full contour relation or give a UV completion.",
        "checks": checks,
    }
