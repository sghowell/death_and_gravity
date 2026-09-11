"""Independent sphere, curvature and local Euler-variation fixtures."""

from functools import cache
from itertools import product
from math import factorial

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_adiabatic_action import (
    action,
    angular,
    audit,
    curvature,
    matching,
)


@cache
def matrices(n):
    raw = s.Matrix(n, n, lambda i, j: (i + j) % 3 - 1)
    F = raw - s.eye(n) * s.trace(raw) / n
    rawG = s.Matrix(n, n, lambda i, j: (i * j + i + j) % 5 - 2)
    G = rawG + s.eye(n) * (s.trace(F * F) - s.trace(rawG)) / n
    return F, G


def trace_values(expr, F, G):
    replacements = {angular.d: F.rows}
    for symbol in expr.free_symbols:
        if str(symbol).startswith("tr_"):
            word = str(symbol)[3:]
            value = s.eye(F.rows)
            for letter in word:
                value = value * ({"F": F, "G": G}[letter])
            replacements[symbol] = s.trace(value)
    return s.factor(expr.subs(replacements))


def literal_sphere_average(polynomial, variables):
    n = len(variables)
    out = 0
    for powers, coefficient in s.Poly(s.expand(polynomial), *variables).terms():
        if any(p % 2 for p in powers):
            continue
        half = sum(powers) // 2
        moment = s.prod(s.factorial2(p - 1) for p in powers) / s.prod(
            n + 2 * j for j in range(half)
        )
        out += coefficient * moment
    return s.factor(out)


WORDS = (
    (("F",),),
    (("G",),),
    (("F", "F"),),
    (("F",), ("F",)),
    (("F",), ("G",)),
    (("F", "F"), ("F",)),
    (("F",), ("F",), ("F",)),
    (("F",), ("F",), ("G",)),
    (("F",), ("F",), ("F",), ("F",)),
    (("F", "G"), ("G", "F")),
    (("G",), ("G",)),
)


@pytest.mark.parametrize("dimension", (3, 4, 5))
@pytest.mark.parametrize("words", WORDS)
def test_independent_expanded_polynomial_sphere_contractions(dimension, words):
    F, G = matrices(dimension)
    variables = s.symbols("n0:" + str(dimension), real=True)
    n = s.Matrix(variables)
    polynomial = s.Integer(1)
    for word in words:
        M = s.eye(dimension)
        for letter in word:
            M = M * ({"F": F, "G": G}[letter])
        polynomial *= (n.T * M * n)[0]
    actual = angular.projection_average(angular.projection_key(words))
    assert trace_values(actual, F, G) == literal_sphere_average(polynomial, variables)


def direct_curvature(F, G, H, Hd):
    n = F.rows
    D = n + 1
    metric = s.diag(-1, *([1] * n))
    first = s.zeros(D)
    first[1:, 1:] = 2 * H * s.eye(n) + F
    second = s.zeros(D)
    second[1:, 1:] = (2 * Hd + 4 * H * H) * s.eye(n) + 4 * H * F + G
    inverse = metric
    inverse_first = -inverse * first * inverse
    gamma = {}
    gamma_dot = {}
    for a, b, c in product(range(D), repeat=3):
        value = derivative = 0
        for j in range(D):
            t0 = (
                (first[j, c] if b == 0 else 0)
                + (first[j, b] if c == 0 else 0)
                - (first[b, c] if j == 0 else 0)
            )
            t1 = (
                (second[j, c] if b == 0 else 0)
                + (second[j, b] if c == 0 else 0)
                - (second[b, c] if j == 0 else 0)
            )
            value += inverse[a, j] * t0 / 2
            derivative += (inverse_first[a, j] * t0 + inverse[a, j] * t1) / 2
        gamma[a, b, c] = value
        gamma_dot[a, b, c] = derivative
    Riem = {}
    for a, b, c, d in product(range(D), repeat=4):
        value = (gamma_dot[a, d, b] if c == 0 else 0) - (
            gamma_dot[a, c, b] if d == 0 else 0
        )
        value += sum(
            gamma[a, c, e] * gamma[e, d, b] - gamma[a, d, e] * gamma[e, c, b]
            for e in range(D)
        )
        Riem[a, b, c, d] = s.expand(value)
    Ric = s.Matrix(D, D, lambda b, d: sum(Riem[a, b, a, d] for a in range(D)))
    scalar = sum(inverse[a, b] * Ric[a, b] for a, b in product(range(D), repeat=2))
    Ric2 = sum(
        inverse[a, a] * inverse[b, b] * Ric[a, b] ** 2
        for a, b in product(range(D), repeat=2)
    )
    Riem2 = sum(
        inverse[a, a]
        * inverse[b, b]
        * inverse[c, c]
        * inverse[d, d]
        * (metric[a, a] * Riem[a, b, c, d]) ** 2
        for a, b, c, d in product(range(D), repeat=4)
    )
    return {"R": scalar, "Ricci2": Ric2, "Riemann2": Riem2}


@pytest.mark.parametrize("dimension", (3, 4, 5))
@pytest.mark.parametrize("name", ("R", "Ricci2", "Riemann2"))
def test_independent_full_nondiagonal_Christoffel_Riemann(dimension, name):
    F, G = matrices(dimension)
    H, Hd = s.Rational(2, 7), -s.Rational(3, 11)
    direct = direct_curvature(F, G, H, Hd)
    target = trace_values(matching.curvatures()[name], F, G).subs(
        {angular.h: H, angular.h1: Hd}
    )
    assert s.factor(direct[name] - target) == 0
    assert F * G != G * F


def mp_matrix(M):
    return mp.matrix([[mp.mpf(str(s.N(x, 85))) for x in row] for row in M.tolist()])


def mp_trace(M):
    return sum(M[j, j] for j in range(M.rows))


def local_density(K, Kd, Kdd, w, wd, wdd, order):
    A = Kd * K**-1
    Ad = Kdd * K**-1 - A * A
    p = wd / w
    pd = wdd / w - p * p
    S = p * mp.eye(K.rows) / 2 - A / 2
    T = (pd - p * p) * mp.eye(K.rows) / 2 - Ad / 2 + p * A / 2
    if order == 0:
        return -K.rows * w / 2
    if order == 2:
        return mp_trace(S * S) / (4 * w)
    if order == 4:
        return mp_trace(T * T + S**4) / (16 * w**3)
    raise ValueError("Unsupported local order")


@cache
def background(which):
    F, G = matrices(3)
    K0 = s.eye(3) if which == 0 else s.diag(1, 4, 9)
    jets = [K0, F / 31, G / 37, (F + G) / 41, (2 * F - G) / 43]
    freq = [s.Integer(1000), s.Integer(2), s.Integer(-3), s.Integer(5), s.Integer(-7)]
    B = s.eye(3) if which == 0 else s.diag(1, 2, 3)
    normalized = {"K" + str(j): B.inv() * jets[j] * B.inv() for j in range(1, 5)}
    substitution = dict(zip(action.w, freq))
    checks = action.variational_check()
    E = {
        order: checks["matrix_Euler"][order].evaluate(normalized).subs(substitution)
        for order in (0, 2, 4)
    }
    Ew = {
        order: checks["frequency_Euler"][order].subs(substitution).subs(angular.d, 3)
        for order in (0, 2, 4)
    }
    for order in (0, 2, 4):
        for sym in list(Ew[order].free_symbols):
            if str(sym).startswith("tr_"):
                encoded = str(sym)[3:]
                letters = ["K" + x for x in encoded.split("K") if x]
                value = s.eye(3)
                for letter in letters:
                    value = value * normalized[letter]
                Ew[order] = Ew[order].subs(sym, s.trace(value))
    return jets, freq, B, E, Ew


@pytest.mark.parametrize("which", (0, 1))
@pytest.mark.parametrize("order", (0, 2, 4))
@pytest.mark.parametrize("direction", ("metric_F", "metric_G", "frequency"))
def test_independent_literal_local_action_Euler_derivative(which, order, direction):
    with mp.workdps(80):
        jets, freq, B, E, Ew = background(which)
        Kjet = list(map(mp_matrix, jets))
        wjet = [mp.mpf(int(x)) for x in freq]
        F, G = matrices(3)
        D = F / 17 if direction == "metric_F" else G / 19

        def perturbed(t, b, slot):
            Ks = [
                sum(
                    (t**j * Kjet[j + n] / factorial(j) for j in range(5 - n)),
                    mp.zeros(3),
                )
                for n in range(3)
            ]
            ws = [
                sum(t**j * wjet[j + n] / factorial(j) for j in range(5 - n))
                for n in range(3)
            ]
            if direction == "frequency":
                ws[slot] += b
            else:
                Ks[slot] += b * mp_matrix(D)
            return local_density(*Ks, *ws, order)

        actual = mp.mpf(0)
        for slot in range(3):
            actual += (-1) ** slot * mp.diff(
                lambda t, b, slot=slot: perturbed(t, b, slot),
                (mp.mpf(0), mp.mpf(0)),
                (slot, 1),
            )
        expected = (
            Ew[order]
            if direction == "frequency"
            else s.trace(E[order] * B.inv() * D * B.inv())
        )
        expected = mp.mpf(str(s.N(expected, 75)))
        assert abs(actual - expected) < mp.mpf("1e-55") * max(1, abs(expected))


@pytest.mark.parametrize("order", (0, 2, 4))
def test_all_noncommuting_word_variation_identities(order):
    c = action.variational_check()
    assert not c["checks"]["matrix_variation_order_" + str(order)]
    assert c["checks"]["frequency_variation_order_" + str(order)] == 0


def test_complete_fourth_order_covariant_matching_before_dimension_limit():
    c = matching.calculation()
    assert c["second_remaining"] == 0 and c["fourth_remaining"] == 0


@pytest.mark.parametrize("key", ("vacuum", "einstein", "curvature"))
def test_independent_MSbar_laurent_finite_coefficients(key):
    e = curvature.epsilon
    L = curvature.ell
    c = curvature.dimensional_coefficients()
    expected = {
        "vacuum": 3 / e + s.Rational(5, 2) - 3 * L,
        "einstein": 1 / e + s.Rational(5, 3) - L,
        "curvature": 2 * c["fixed_vector"] / e
        - 2 * L * c["fixed_vector"]
        - 4 * c["fixed_scalar"],
    }
    assert s.expand(c[key] - expected[key]) == 0


def test_omitting_evanescent_polarizations_changes_finite_terms():
    e = curvature.epsilon
    weight = s.exp(e * (s.EulerGamma - curvature.ell))
    wrong = s.series(6 * weight * s.gamma(e - 2), e, 0, 1).removeO().expand()
    assert s.simplify(wrong - curvature.dimensional_coefficients()["vacuum"]) == 2


@pytest.mark.parametrize(
    "name,value",
    audit.residuals().items(),
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_exact_identity(name, value):
    assert s.cancel(value) == 0, name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda x: x if isinstance(x, str) else None
)
def test_unsupported_scope(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", tuple(audit.packets()))
def test_all_packet_gates(name):
    assert all(audit.packets()[name]["gates"].values())


def test_counts_and_original_frontier_open():
    assert len(audit.residuals()) == 172 and audit.scalar_entry_count() == 172
    assert len(audit.gates()) == 36 and len(audit.controls()) == 9
    assert audit.rejected_inputs() == 115 and len(audit.frontier()) == 9
    assert audit.matching()[-1] == audit.ITEM


@cache
def independent_root_jets(which):
    with mp.workdps(80):
        jets, _, B, _, _ = background(which)
        normalized = [mp.eye(3)] + [
            mp_matrix(B.inv() * jets[j] * B.inv()) for j in range(1, 5)
        ]

        def root(t):
            K = sum(
                (t**j * normalized[j] / factorial(j) for j in range(5)), mp.zeros(3)
            )
            vals, O = mp.eigsy((K + K.T) / 2)
            return O * mp.diag([mp.sqrt(x) for x in vals]) * O.T

        return {
            "B": [mp.diff(root, mp.mpf(0), j) for j in range(5)],
            "inverseB": [
                mp.diff(lambda t: root(t) ** -1, mp.mpf(0), j) for j in range(5)
            ],
        }


@pytest.mark.parametrize("which", (0, 1))
@pytest.mark.parametrize("order", range(5))
@pytest.mark.parametrize("kind", ("B", "inverseB"))
def test_independent_positive_root_and_inverse_raw_jets(which, order, kind):
    with mp.workdps(65):
        jets, _, B, _, _ = background(which)
        normalized = {"K" + str(j): B.inv() * jets[j] * B.inv() for j in range(1, 5)}
        exact = action.ordered_frame()[kind][order].evaluate(normalized)
        actual = independent_root_jets(which)[kind][order]
        assert mp.norm(actual - mp_matrix(exact)) < mp.mpf("1e-40")


@pytest.mark.parametrize("dimension", (3, 4, 5))
@pytest.mark.parametrize("momentum", (0, 1000, 10**20))
def test_complete_actual_dimension_dependent_Hamiltonian_constraint(
    dimension, momentum
):
    n = dimension
    v = s.Matrix(range(1, n + 1))
    O = s.eye(n) - 2 * v * v.T / (v.dot(v))
    E = O * s.diag(1, 4, s.Rational(1, 4), *([1] * (n - 3))) * O.T
    Q = E.inv()
    k = momentum * v
    a, m = s.Rational(25, 16), s.Integer(1000)
    K = a ** (2 - n) * E + k * k.T / (a**n * m**2)
    V = a ** (n - 2) * m**2 * Q + a ** (n - 4) * (
        (k.T * Q * k)[0] * Q - Q * k * k.T * Q
    )
    omega2 = m * m + (k.T * Q * k)[0] / a**2
    assert (K * V - omega2 * s.eye(n)).applyfunc(s.cancel) == s.zeros(n)
    assert E.det() == 1
    if momentum:
        assert ((a ** (2 - n) * E) * V - omega2 * s.eye(n)) != s.zeros(n)


def test_evanescent_Einstein_factor_cannot_be_frozen_before_subtraction():
    e = curvature.epsilon
    weight = s.exp(e * (s.EulerGamma - curvature.ell))
    wrong = s.series(-weight * s.gamma(e - 1), e, 0, 1).removeO().expand()
    assert s.simplify(
        curvature.dimensional_coefficients()["einstein"] - wrong
    ) == s.Rational(2, 3)


def test_evanescent_curvature_contact_cannot_be_discarded():
    e = curvature.epsilon
    c = curvature.dimensional_coefficients()
    weight = s.exp(e * (s.EulerGamma - curvature.ell))
    wrong = (
        s.series(2 * weight * s.gamma(e) * c["fixed_vector"], e, 0, 1)
        .removeO()
        .expand()
    )
    assert s.simplify(wrong - c["curvature"]) == 4 * c["fixed_scalar"]


def test_commuting_trace_replacement_changes_oriented_contractions():
    F, G = matrices(3)
    assert s.trace(F * G * F * G) != s.trace(F * F * G * G)
    expr = angular.projection_average((("F", "G"), ("G", "F")))
    changed = expr.subs(
        angular.trace_word(("F", "G", "F", "G")),
        angular.trace_word(("F", "F", "G", "G")),
    )
    assert trace_values(expr - changed, F, G) != 0


@pytest.mark.parametrize("rates,accelerations", (([], []), ([1], []), ([], [1])))
def test_invalid_diagonal_curvature_inputs_rejected(rates, accelerations):
    with pytest.raises(ValueError):
        curvature.diagonal_curvature_jets(rates, accelerations)


def test_matching_not_full_feedback_or_original_closure():
    assert "not established" in audit.observable()["boundary"]
    assert "unchanged" in audit.observable()["conclusion"]
    assert all(row["status"] != "COMPLETE" for row in audit.frontier())
