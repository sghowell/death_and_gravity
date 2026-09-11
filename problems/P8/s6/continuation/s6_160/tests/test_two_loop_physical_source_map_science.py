"""Independent nonlinear pullbacks, coupled Ward identities and controls."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_two_loop_physical_source_map import (
    audit,
    calibration,
    topology,
    transport,
    vertices,
    ward,
)


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_invalid_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("r", (".02", ".2", "1"))
@pytest.mark.parametrize("n", (0, 2, 4, 8))
def test_direct_nonlinear_Gaussian_pullback(r, n):
    with mp.workdps(42):
        r, K, h = mp.mpf(r), mp.mpf(3), mp.mpf(".7")

        def transformed(x):
            F = x + r * x**3
            return (1 + 3 * r * x * x) * F**n * mp.exp(-K * F * F / (2 * h))

        integral = 2 * mp.quad(transformed, [0, mp.mpf(".5"), 1, 2, mp.inf])
        moment = (
            mp.factorial(n)
            / (2 ** (n // 2) * mp.factorial(n // 2))
            * (h / K) ** (n // 2)
        )
        expected = mp.sqrt(2 * mp.pi * h / K) * moment
        assert abs(integral / expected - 1) < mp.mpf("1e-32")


@cache
def coupled_map_data():
    x, y, r = s.symbols("x y r")
    # A finite-difference-like cubic coupling between sites, not a
    # product of independent scalar maps.
    fx = x + r * (x**3 / 2 + (x - y) ** 3 / 5)
    fy = y + r * (y**3 / 3 + (y - x) ** 3 / 5)
    change = s.expand((2 * fx**2 + 3 * fy**2 - 2 * x * x - 3 * y * y) / 2)
    w1, w2 = change.coeff(r, 1), change.coeff(r, 2)
    weight = 1 - r * w1 + r * r * (w1 * w1 / 2 - w2) + r**3 * (w1 * w2 - w1**3 / 6)
    jac = s.det(s.Matrix([fx, fy]).jacobian([x, y]))

    def mean(poly):
        result = 0
        for (a, b), c in s.Poly(s.expand(poly), x, y).terms():
            if a % 2 == 0 and b % 2 == 0:
                result += (
                    c
                    * s.factorial2(a - 1)
                    * s.factorial2(b - 1)
                    / s.Integer(2) ** (a // 2)
                    / s.Integer(3) ** (b // 2)
                )
        return s.factor(result)

    rows = {}
    for a, b in ((0, 0), (2, 0), (0, 2), (1, 1), (2, 2), (4, 0)):
        expr = s.expand(jac * weight * fx**a * fy**b)
        rows[(a, b)] = [
            mean(expr.coeff(r, j)) - (mean(x**a * y**b) if j == 0 else 0)
            for j in range(4)
        ]
    missing = mean(s.expand(weight * fx**2).coeff(r, 1))
    return rows, missing


@pytest.mark.parametrize("a,b", ((0, 0), (2, 0), (0, 2), (1, 1), (2, 2), (4, 0)))
def test_coupled_two_variable_third_order_Ward_identity(a, b):
    rows, missing = coupled_map_data()
    assert rows[(a, b)] == [0, 0, 0, 0]
    assert missing != 0


@pytest.mark.parametrize("degree", (1, 3, 5))
def test_formal_inverse_and_physical_source_composition(degree):
    z, r = s.symbols("z r")
    inv = z - r * z**3 + 3 * r * r * z**5 - 12 * r**3 * z**7
    composed = s.series(inv + r * inv**3, r, 0, 4).removeO()
    assert s.expand(composed - z) == 0
    assert s.series(composed**degree - z**degree, r, 0, 4).removeO() == 0


def test_octic_control_from_independent_shifted_Gaussian():
    J, h, K, q = s.symbols("J h K q", positive=True)
    # E_J[x^8] for mean hJ/K and covariance h/K.
    mu, v = h * J / K, h / K
    eighth = (
        mu**8 + 28 * mu**6 * v + 210 * mu**4 * v * v + 420 * mu * mu * v**3 + 105 * v**4
    )
    defect = s.diff(q * eighth / (4 * h), J, 4).subs(J, 0)
    assert defect == 1260 * q * h**5 / K**6
    assert s.degree(defect, h) == 5  # E=4 => h^(E-1+L), L=2.


def test_Yukawa_control_from_normalized_cumulant():
    Y, r, v, m = s.symbols("Y r v m", positive=True)
    moments = {
        0: 1 - Y * v / m**2 + 6 * Y * r * v * v / m**2,
        2: v - 3 * Y * v * v / m**2 + 30 * Y * r * v**3 / m**2,
        4: 3 * v * v - 15 * Y * v**3 / m**2 + 210 * Y * r * v**4 / m**2,
    }
    cumulant = moments[4] / moments[0] - 3 * (moments[2] / moments[0]) ** 2
    coefficient = s.diff(cumulant, Y, r).subs({Y: 0, r: 0})
    assert s.factor(coefficient) == 48 * v**4 / m**2


@pytest.mark.parametrize("E", (0, 2, 4))
@pytest.mark.parametrize("L", (0, 1, 2))
def test_independent_integer_loop_support(E, L):
    # For pure even bosons include action map replacements and cubic
    # physical-source endpoints; verify the halfedge degree directly.
    bound = topology.max_map_degree(E, L)
    observed = []
    for base4 in range(4):
        for source_R in range(E + 1):
            for action_R in range(4):
                p = source_R + action_R
                effective = 1 - E // 2 + base4 + p
                if effective == L and p >= 0:
                    observed.append(p)
                    assert p <= bound
    if L - 1 + E // 2 >= 0:
        assert max(observed) == bound


def test_all_required_counterterm_and_Yukawa_support_present():
    rows = {
        (r["operator"], r["counterterm_loop_weight"], r["map_degree"])
        for r in topology.support_rows()
    }
    assert {
        ("Phi_quartic", 0, 2),
        ("Phi_quadratic", 1, 2),
        ("Phi_quartic", 1, 1),
        ("Phi_quadratic", 2, 1),
        ("Yukawa", 0, 1),
        ("Yukawa", 1, 1),
    } <= rows
    assert len(rows) == len(topology.support_rows()) == 42


@pytest.mark.parametrize("o1,o2", ((2, 3), (-2, 7), (10**30, -(10**60))))
def test_formal_LSZ_cancellation_does_not_assume_small_ordinary_overlap(o1, o2):
    h = s.Symbol("h")
    o = 1 + h * o1 + h * h * o2
    A = 2 + 3 * h + 5 * h * h
    residue = s.series(o**2, h, 0, 3).removeO()
    numerator = s.series(o**4 * A, h, 0, 3).removeO()
    assert s.series(numerator / residue**2 - A, h, 0, 3).removeO() == 0


def test_common_D_dimensional_coefficient_lift_preserves_relation():
    mu, e, lam, gamma = s.symbols("mu epsilon lambda gamma", positive=True)
    c = 4 * lam * lam / gamma
    assert (
        s.simplify(
            4 * (mu ** (2 * e) * lam) ** 2 / (mu ** (2 * e) * gamma) - mu ** (2 * e) * c
        )
        == 0
    )


def test_actual_literal_map_and_parent_enclosures():
    d = calibration.data()
    assert d["literal_R_jet_monomial_count"] == 61
    assert d["literal_R_field_degrees"] == [3]
    assert max(d["literal_R_total_derivative_degrees"]) == 4
    assert all(d["bounds"].values())
    assert len(vertices.data()["all_scalar_generated_degrees"]) == 6
    assert len(ward.data()["omission_controls"]) == 5
    assert "No numerical bound" in transport.data()["not_claimed"]


def test_scope_preserves_all_broad_obligations():
    m = audit.matching()
    assert len(m) == 16
    assert sum(r["status"] == "OPEN" for r in m) == 5
    assert m[-1] == audit.ITEM
    assert audit.controls()["original_P8_not_closed"] is True
