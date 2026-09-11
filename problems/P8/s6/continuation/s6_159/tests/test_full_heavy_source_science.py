"""Independent fixed-reference, regulated-source and mass-derivative checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_full_heavy_source import (
    audit,
    bounds,
    calibration,
    differentiation,
    ownership,
    source,
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


@pytest.mark.parametrize("L,g,G", ((2, 3, 5), (3, 2, 7), (1, 4, 2)))
def test_fixed_counterterm_vacuum_derivative_equals_covariance_source(L, g, G):
    a = s.Symbol("a", positive=True)
    Tfun = s.Function("T")
    Sfun = s.Function("S")
    T, I2, Sa, B0, r = s.symbols("T I2 Sa B0 r")
    dz = -r
    dm = -L * T / 2 + g * B0 - r
    vacuum = L * Tfun(a) ** 2 / 8 - g * Sfun(a) / 4 + (dm - a * dz) * Tfun(a) / 2
    derivative = G * s.diff(vacuum, a).subs(a, 1)
    derivative = derivative.subs(
        {
            Tfun(1): T,
            s.Subs(s.Derivative(Tfun(a), a), a, 1): -I2,
            s.Subs(s.Derivative(Sfun(a), a), a, 1): Sa,
        }
    )
    assert s.expand(derivative - G * (-g * Sa / 2 - g * B0 * I2 + r * T) / 2) == 0


def test_retuning_OS_counterterms_with_background_changes_the_source():
    a, L, g, G, T, I2, B0, r, b1, r1, Sa = s.symbols("a L g G T I2 B0 r b1 r1 Sa")
    Ta = T - I2 * (a - 1)
    Ba = B0 + b1 * (a - 1)
    ra = r + r1 * (a - 1)
    raw = L * Ta**2 / 8 - g * (s.Symbol("S0") + Sa * (a - 1)) / 4
    fixed = (-L * T / 2 + g * B0 - r + a * r) * Ta / 2
    moving = (-L * Ta / 2 + g * Ba - a * ra + a * ra) * Ta / 2
    defect = s.expand(G * s.diff(moving - fixed, a).subs(a, 1))
    assert defect != 0
    assert s.diff(defect, L) == G * T * I2 / 4
    assert s.diff(G * s.diff(raw + fixed, a).subs(a, 1), L) == 0


@pytest.mark.parametrize(
    "k0,k1", ((".2", ".3"), ("-.2", ".3"), (".2", "-.4"), ("0", ".3"))
)
def test_full_bare_first_source_map_on_complex_regulator_circle(k0, k1):
    with mp.workdps(42):
        G, mu, Q = mp.mpf(".3"), mp.mpf(5), 16 * mp.pi**2
        k0, k1 = mp.mpf(k0), mp.mpf(k1)
        finite, pole = 0, 0
        for j in range(32):
            e = mp.mpf(".02") * mp.exp(2j * mp.pi * j / 32)
            Tad = mp.exp(mp.euler * e) * mu ** (2 * e) * mp.gamma(e - 1) / Q

            def first_source_reexpanded(h, e=e, Tad=Tad):
                gs = G - h * (k0 + e * k1) * G
                return -gs * Tad / 2

            second = mp.diff(first_source_reexpanded, 0)
            finite += second / 32
            pole += e * second / 32
        expected = -G * (k0 * (2 * mp.log(mu) + 1) + k1) / (2 * Q)
        assert abs(finite - expected) < mp.mpf("1e-34")
        assert abs(pole + G * k0 / (2 * Q)) < mp.mpf("1e-36")
        early = -G * k0 * (2 * mp.log(mu) + 1) / (2 * Q)
        assert abs(finite - early) > mp.mpf("1e-5")


@pytest.mark.parametrize("k1", (s.Rational(1, 3), s.Rational(-2, 7)))
def test_positive_epsilon_source_pole_product_is_required(k1):
    e, G, Q, k0, ell = s.symbols("e G Q k0 ell")
    T = -1 / (Q * e) - (ell + 1) / Q
    full = s.expand((k0 + e * k1) * G * T / 2).coeff(e, 0)
    early = s.expand(k0 * G * T / 2).coeff(e, 0)
    assert s.expand(full - early) == -G * k1 / (2 * Q)


@pytest.mark.parametrize("B0,Bp,rval", ((2, 3, 5), (1, 4, 2)))
def test_literal_OS_covariance_insertion_retains_correlated_slope(B0, Bp, rval):
    u, g = s.symbols("u g", positive=True)
    B = B0 + Bp * (u + 1) + rval * (u + 1) ** 2
    Pi = -g * B
    PiR = Pi - Pi.subs(u, -1) - (u + 1) * s.diff(Pi, u).subs(u, -1)
    inserted = -PiR / (u + 1) ** 2
    r = -g * Bp
    assert s.cancel(inserted - (g * (B - B0) / (u + 1) ** 2 + r / (u + 1))) == 0
    without_r = g * (B - B0) / (u + 1) ** 2
    assert s.cancel(inserted - without_r) != 0


@cache
def gauss_nodes():
    with mp.workdps(42):
        nodes, weights = mp.gauss_quadrature(40, "legendre")
        return tuple(((x + 1) / 2, w / 2) for x, w in zip(nodes, weights))


def heavy_position(r, t, M, pos):
    if pos == 0:
        return M + r + r * t, r + r * t, M, mp.mpf(0)
    if pos == 1:
        return 1 + M * r + r * t, 1 + r * t, mp.mpf(1), mp.mpf(1)
    return 1 + r + M * r * t, 1 + r, mp.mpf(1), mp.mpf(1)


@pytest.mark.parametrize("M", (2, 3))
def test_two_light_mass_derivative_against_zero_dimensional_Gaussian(M):
    with mp.workdps(36):
        result = 0
        for pos in range(3):
            for r, wr in gauss_nodes():
                for t, wt in gauss_nodes():
                    B, C, _, _ = heavy_position(r, t, mp.mpf(M), pos)
                    # Gamma(3)=2, derivative power -3, two equal-light labels.
                    result -= 12 * wr * wt * r * C / B**4
        assert abs(result + mp.mpf(2) / M) < mp.mpf("1e-24")


@pytest.mark.parametrize("M", (2, 100, 100000))
@pytest.mark.parametrize("angle", (0, 1, 2, 3))
def test_explicit_corner_difference_and_split_bound(M, angle):
    with mp.workdps(34):
        M = mp.mpf(M)
        e = mp.exp(2j * mp.pi * angle / 4) / 16
        for pos in range(3):
            for t in (mp.mpf(".2"), mp.mpf(".8")):
                for r in (1 / (10 * M), 1 / (3 * M), mp.mpf(".4"), mp.mpf(1)):
                    B, C, B0, C0 = heavy_position(r, t, M, pos)
                    Fa = (1 - 2 * e) * C * (1 + t + r * t) ** (e - 2) * B ** (-2 * e)
                    F0 = (1 - 2 * e) * C0 * (1 + t) ** (e - 2) * B0 ** (-2 * e)
                    bound = min(10 * M * r, 5) * (3 * M) ** (mp.mpf(1) / 8)
                    assert abs(Fa - F0) < bound


@pytest.mark.parametrize("angle", (0, 1, 2, 3))
def test_full_corner_subtracted_mass_derivative_circle_enclosure(angle):
    with mp.workdps(34):
        e = mp.exp(2j * mp.pi * angle / 4) / 16
        M, mu, Q = mp.mpf(7), mp.mpf(2), 16 * mp.pi**2
        J = (2 ** (e - 1) - 1) / (e - 1)
        total = 0
        for pos in range(3):
            _, _, B0, C0 = heavy_position(0, mp.mpf(".5"), M, pos)
            total += 2 * (1 - 2 * e) * C0 * B0 ** (-2 * e) * J / e
            for r, wr in gauss_nodes():
                for t, wt in gauss_nodes():
                    B, C, B0, C0 = heavy_position(r, t, M, pos)
                    Fa = (1 - 2 * e) * C * (1 + t + r * t) ** (e - 2) * B ** (-2 * e)
                    F0 = (1 - 2 * e) * C0 * (1 + t) ** (e - 2) * B0 ** (-2 * e)
                    total += 2 * wr * wt * r ** (e - 1) * (Fa - F0)
        assert abs(total) < 700 * (3 * M) ** (mp.mpf(3) / 16)
        derivative = (
            mp.exp(2 * mp.euler * e)
            * mu ** (4 * e)
            * mp.gamma(2 * e - 1)
            * total
            / Q**2
        )
        assert (
            abs(derivative)
            < 28000 * mu ** (mp.mpf(1) / 4) * (3 * M) ** (mp.mpf(3) / 16) / Q**2
        )


def test_zero_cubic_limit_and_exact_power_caps():
    d = bounds.scalar_source_enclosure(2, 2, 1, 0, 144, 2, 2)
    assert d["complete_scalar_and_cubic_source_absolute_upper"] == 0
    assert bounds.field_reexpression_upper(0, 1, 1, 1, 144) == 0
    assert source.data()["full_first_parameter_source_reexpression"] != 0
    assert len(differentiation.data()["three_heavy_positions_each_twice"]) == 3


def test_complete_named_reference_scope_not_original_closure():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(ownership.rows()) == 8
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 15
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 10
    pending = {r["id"]: r["status"] for r in audit.matching()}
    assert pending[audit.TARGET] == audit.STATUS
    for name in (
        "remaining_scale_parameter_field_map_and_cross_terms",
        "finite_EFT_higher_order_truncation",
        "V_contour_and_cut_control",
        "finite_gravity_G",
        "common_parent_B",
    ):
        assert pending[name] == "OPEN"
    assert len(audit.residuals()) == 44
    assert len(audit.gates()) == 23
    assert audit.controls()["rejected_inputs"] == 138
