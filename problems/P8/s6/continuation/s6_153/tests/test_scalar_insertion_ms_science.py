"""Independent dimensional slopes, finite pole products and exact reference shifts."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_scalar_insertion_ms import (
    audit,
    bounds,
    calibration,
    forward,
    inner,
    outer,
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


def intervals(M):
    return [0, 1 / M, 1]


def alpha_D(e, M, ell):
    return (
        mp.exp(mp.euler * e + ell * e)
        * mp.gamma(1 + e)
        * mp.quad(
            lambda x: x * (1 - x) / (x * M + (1 - x) ** 2) ** (1 + e), intervals(M)
        )
    )


def Pi_D(z, e, M, ell):
    return (
        mp.exp(mp.euler * e + ell * e)
        * mp.gamma(e)
        * mp.quad(lambda x: (x * M + 1 - x - x * (1 - x) * z) ** (-e), intervals(M))
    )


@pytest.mark.parametrize("M", (2, 10, 40, 1000))
@pytest.mark.parametrize("e", (".17", ".4"))
def test_dimensional_slope_is_fixed_scale_bubble_derivative(M, e):
    with mp.workdps(36):
        M, e, ell = mp.mpf(M), mp.mpf(e), mp.mpf(2)
        derivative = mp.diff(lambda z: Pi_D(z, e, M, ell), 1)
        assert abs(derivative - alpha_D(e, M, ell)) < mp.mpf("1e-26")


@pytest.mark.parametrize("M", (2, 40, 1000))
@pytest.mark.parametrize("ell", (0, 2, 10))
def test_complex_finite_part_retains_alpha_epsilon_coefficient(M, ell):
    with mp.workdps(36):
        M = mp.mpf(M)
        b = lambda x: x * (1 - x) / (x * M + (1 - x) ** 2)
        A0 = mp.quad(b, intervals(M))
        A1 = mp.quad(
            lambda x: b(x) * (ell - mp.log(x * M + (1 - x) ** 2)), intervals(M)
        )
        wanted = mp.quad(
            lambda x: b(x) * (2 * ell - mp.log(x * M + (1 - x) ** 2)), intervals(M)
        )
        pole, finite = 0, 0
        for k in range(24):
            e = mp.mpf(".04") * mp.exp(2j * mp.pi * k / 24)
            # e times the complete alpha_D I0 product, in Q^2/g units.
            regular = (
                mp.exp(2 * mp.euler * e + 2 * ell * e)
                * mp.gamma(1 + e) ** 2
                * mp.quad(
                    lambda x, e=e: b(x) * (x * M + (1 - x) ** 2) ** (-e), intervals(M)
                )
            )
            pole += regular / 24
            finite += regular / e / 24
        assert abs(pole - A0) < mp.mpf("1e-25")
        assert abs(finite - wanted) < mp.mpf("1e-24")
        assert abs(wanted - (ell * A0 + A1)) < mp.mpf("1e-28")


@pytest.mark.parametrize("M", (2, 40, 1000))
@pytest.mark.parametrize("e", (".2", ".4"))
def test_inner_OS_asymptote_retains_alpha_D(M, e):
    with mp.workdps(35):
        M, e, ell = mp.mpf(M), mp.mpf(e), mp.mpf(2)
        alpha = alpha_D(e, M, ell)
        p1 = Pi_D(1, e, M, ell)
        for factor in (mp.mpf("1e4"), mp.mpf("1e8")):
            t = M * factor
            p = Pi_D(-t, e, M, ell)
            Q = (p - p1 + (t + 1) * alpha) / (t + 1)
            assert 0 < Q < alpha
            assert 0 < alpha - Q <= p1 / (t + 1)
        assert abs(Q - alpha) / alpha < mp.mpf("1e-4")


@pytest.mark.parametrize("M", (32, 40, 1000, 1000000))
def test_positive_parameter_enclosure_at_scale_above_M(M):
    with mp.workdps(35):
        M = mp.mpf(M)
        ell = mp.log(M) + 2
        b = lambda x: x * (1 - x) / (x * M + (1 - x) ** 2)
        A0 = mp.quad(b, intervals(M))
        F = mp.quad(
            lambda x: b(x) * (2 * ell - mp.log(x * M + (1 - x) ** 2)), intervals(M)
        )
        assert 0 < A0 < 1 / (2 * M)
        assert ell * A0 <= F <= 2 * ell * A0 <= ell / M


@pytest.mark.parametrize("first", (1, -2, s.Rational(3, 7)))
def test_early_alpha_truncation_loses_finite_reference(first):
    e, ell, a0 = s.symbols("e ell a0")
    alpha = a0 + e * first
    I0 = 1 / e + ell + e * (ell**2 / 2 + s.pi**2 / 12)
    defect = s.limit((alpha - a0) * I0, e, 0)
    assert defect == first != 0


@pytest.mark.parametrize("M", (10, 40, 100))
def test_C_squared_forward_coefficient_from_independent_contour(M):
    with mp.workdps(36):
        L, g = mp.mpf(2), mp.mpf(".5")
        D = M - 2
        value = 0
        for k in range(32):
            nu = mp.exp(2j * mp.pi * k / 32)
            amp = (
                (-L + g / (D - nu)) ** 2 + (-L + g / (D + nu)) ** 2 + (-L + g / M) ** 2
            )
            value += amp / nu**2 / 32
        btree = 2 * g / D**3
        expected = btree * (-2 * L + 3 * g / D)
        assert abs(value - expected) < mp.mpf("1e-24")
        assert -2 * L < expected / btree < 0


@pytest.mark.parametrize("F", (s.Rational(1, 5), s.Rational(-2, 7)))
def test_literal_bare_coordinate_shift(F):
    L, G, M, z, h = s.symbols("L G M z h")
    g = G**2
    shiftedG = G - h**2 * L * G * F
    shiftedM = M - h**2 * g * F
    shiftedL = L - h**2 * 3 * L**2 * F
    invariants = (z, 4 - z, 0)
    before = -L + g * sum(1 / (M - q) for q in invariants)
    after = -shiftedL + shiftedG**2 * sum(1 / (shiftedM - q) for q in invariants)
    second = s.diff(after - before, h, 2).subs(h, 0) / 2
    target = F * sum((-L + g / (M - q)) ** 2 for q in invariants)
    assert s.factor(second - target) == 0


def test_inner_OS_annihilates_all_affine_mass_and_source_terms():
    z, alpha, beta = s.symbols("z alpha beta")
    f = alpha + beta * z
    assert s.expand(f - f.subs(z, 1) - (z - 1) * s.diff(f, z)) == 0
    assert outer.data()["checks"]["constant_tadpole_inner_OS_annihilation"] == 0
    assert inner.data()["checks"]["same_actual_parent_normalized_weight"] == 0


def test_four_raw_families_and_no_global_matching_overclaim():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(d["four_raw_scalar_interaction_family_bounds"]) == 4
    assert d["four_raw_scalar_MS_interaction_families_upper"] == sum(
        d["four_raw_scalar_interaction_family_bounds"].values()
    )
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 11
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 4
    assert len(audit.gates()) == 23
    assert audit.controls()["rejected_inputs"] == 90
    assert forward.data()["checks"]["fundamental_G_reference_is_linear"] == 0
    assert bounds.enclosure(1, 1, 40, 1)["finite_F_alpha_upper"] > 0
