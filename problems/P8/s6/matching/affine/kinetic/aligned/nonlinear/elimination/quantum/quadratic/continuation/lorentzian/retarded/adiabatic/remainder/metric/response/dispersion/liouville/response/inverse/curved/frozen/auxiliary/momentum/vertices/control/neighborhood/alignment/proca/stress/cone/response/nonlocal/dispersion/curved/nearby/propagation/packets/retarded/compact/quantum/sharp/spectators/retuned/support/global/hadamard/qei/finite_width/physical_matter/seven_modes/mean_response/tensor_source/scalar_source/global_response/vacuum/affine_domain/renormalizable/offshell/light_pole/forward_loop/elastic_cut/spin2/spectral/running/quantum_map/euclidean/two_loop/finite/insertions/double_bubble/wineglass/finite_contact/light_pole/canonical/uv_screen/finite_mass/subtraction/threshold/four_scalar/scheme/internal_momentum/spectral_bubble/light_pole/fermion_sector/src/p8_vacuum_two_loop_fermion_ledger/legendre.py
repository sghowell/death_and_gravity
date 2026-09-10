"""Explicit Euclidean Gaussian and Legendre cancellation of reducible terms."""

from functools import cache

import sympy as sp


def gaussian_moment(power, a):
    if power % 2:
        return sp.Integer(0)
    return sp.factorial2(power - 1) / a ** (power // 2) if power else sp.Integer(1)


@cache
def data():
    a = sp.symbols("positive_tree_Hessian", positive=True)
    b, c, f1, f2, eta = sp.symbols(
        "tree_cubic tree_quartic fermion_first_derivative fermion_second_derivative fluctuation"
    )
    odd = b * eta**3 / 6 + f1 * eta
    even = c * eta**4 / 24 + f2 * eta**2 / 2
    polynomial = sp.Poly(sp.expand(-even + odd**2 / 2), eta)
    W2 = sum(coef * gaussian_moment(power[0], a) for power, coef in polynomial.terms())
    shift = (f1 + b / (2 * a)) / a
    Gamma2 = sp.factor(-W2 + a * shift**2 / 2)
    expected = c / (8 * a * a) - b * b / (12 * a**3) + f2 / (2 * a)
    D = sp.Matrix(
        [[sp.Rational(2), sp.Rational(1, 3)], [sp.Rational(1, 3), sp.Rational(3)]]
    )
    F = sp.Matrix(sp.symbols("first_derivative1 first_derivative2"))
    tad = sp.Matrix(sp.symbols("bosonic_tadpole1 bosonic_tadpole2"))
    connected = (F.T * D * F)[0] / 2 + (F.T * D * tad)[0] + (tad.T * D * tad)[0] / 2
    correction = ((F + tad).T * D * (F + tad))[0] / 2
    checks = {
        "Gaussian_second_moment": gaussian_moment(2, a) - 1 / a,
        "Gaussian_fourth_moment": gaussian_moment(4, a) - 3 / a**2,
        "Gaussian_sixth_moment": gaussian_moment(6, a) - 15 / a**3,
        "connected_saddle_coefficient": sp.factor(
            W2
            + c / (8 * a * a)
            + f2 / (2 * a)
            - 5 * b * b / (24 * a**3)
            - b * f1 / (2 * a * a)
            - f1 * f1 / (2 * a)
        ),
        "Legendre_transform_removes_reducible_terms": sp.factor(Gamma2 - expected),
        "no_two_fermion_dumbbell_survives": sp.diff(Gamma2, f1, 2),
        "no_fermion_bosonic_tadpole_dumbbell_survives": sp.diff(Gamma2, f1, b),
        "mixed_effective_action_coefficient": sp.diff(Gamma2, f2) - 1 / (2 * a),
        "two_component_source_bilinear_cancellation": sp.expand(connected - correction),
        "tree_background_shift_formula": sp.factor(a * shift - f1 - b / (2 * a)),
    }
    return {
        "connected_saddle_two_loop_coefficient": W2,
        "classical_to_mean_background_shift": shift,
        "Legendre_effective_action_two_loop_coefficient": Gamma2,
        "mixed_fermion_part": f2 / (2 * a),
        "finite_index_source_bilinear": connected,
        "finite_index_Legendre_cancellation": correction,
        "arbitrary_index_identity": "After integrating the quadratic fermions, F=-Tr log D_F enters the bosonic action at order h. At order h^2 its genuinely mixed contribution to the bosonic 1PI action is one half Tr[(S_B^(0)'')^-1 F'']. Reducible F' and bosonic-tadpole bilinears cancel under the Legendre transform.",
        "scope": "A formal loop expansion with regulator retained. This is an ownership identity, not evaluation or an error bound for its continuum diagrams.",
        "checks": checks,
    }
