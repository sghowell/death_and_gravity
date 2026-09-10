"""Independent mixed-mass gaps, forest references and selected-order checks."""

from fractions import Fraction as F

import pytest
import sympy as sp
from p8_vacuum_fermion_spectral_pole import (
    audit,
    calibration,
    enclosure,
    kernel,
    quadratic,
    selected,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_explicit_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_all_unsupported_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("m,M", ((2, 1), (36, 2), (36, 10**200), (10**200, 10**197)))
@pytest.mark.parametrize(
    "Y,g,Q", ((0, 1, 144), (1, 0, 144), (F(1, 17), F(1, 31), 144), (1, 2, 256))
)
def test_family_bounds_with_independent_fractions(m, M, Y, g, Q):
    m, M, Y, g, Q = map(F, (m, M, Y, g, Q))
    T = 4 * m * m
    C = 12 * Y / Q
    d = enclosure.bound(m, Y, g, M, Q)
    slope = 4 * C * g / (Q * T)
    B = 4 * C * g / (3 * Q * T * T)
    assert d["finite_positive_outer_slope_upper"] == sp.Rational(slope)
    assert d["uniform_outer_OS_remainder_coefficient_upper"] == sp.Rational(B)
    assert B == Y * g / (Q * Q * m**4)
    assert slope == B * 3 * T
    assert d["strict_family_slope_and_spacelike_remainder_positive"] == (
        Y > 0 and g > 0
    )


@pytest.mark.parametrize("u,M", ((16, 1), (36, 2), (36, 10**200), (10**400, 10**197)))
@pytest.mark.parametrize("x", (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)))
@pytest.mark.parametrize("real,imag", ((-(10**100), 10**200), (1, 0), (4, 10**100)))
def test_unbounded_complex_parameter_gap_independently(u, M, x, real, imag):
    u, M = map(F, (u, M))
    A = x * (1 - x)
    reD = x * u + (1 - x) * M - A * real
    imD = -A * imag
    assert reD >= x * (u - 4) + (1 - x) * M > 0
    assert reD - x * (u - 4) - (1 - x) * M == A * (4 - real) + 4 * x * x
    assert A * A / (reD * reD + imD * imD) <= (1 - x) ** 2 / (u - 4) ** 2


@pytest.mark.parametrize(
    "u,M,x,t",
    ((16, 1, F(1, 4), 0), (36, 2, F(1, 2), 1), (10**400, 10**197, F(1, 4), 10**800)),
)
def test_global_spacelike_parameter_identity(u, M, x, t):
    u, M, x, t = map(F, (u, M, x, t))
    A = x * (1 - x)
    Delta1 = x * u + (1 - x) * M - A
    b = A / Delta1
    d = t + 1
    assert Delta1 > 0 and b > 0
    assert x * u + (1 - x) * M + A * t == Delta1 * (1 + d * b)


def test_mixed_trace_has_only_one_light_insertion():
    d = quadratic.data()
    assert d["placement_count"] == 1
    P = d["inserted_light_covariance"]
    phi = sp.Matrix(sp.symbols("phi1 phi2", real=True))
    G = sp.Symbol("cubic_G", real=True)
    a, b, c = sp.symbols("positive_K11 K12 positive_K22", real=True)
    H = sp.Matrix([[a, b], [b, c]]).inv()
    for i in range(2):
        for j in range(2):
            derivative = sp.diff(d["mixed_heavy_light_trace"], phi[i], phi[j])
            expected = -G * G * P[i, j] * H[i, j]
            assert sp.factor(derivative - expected) == 0
            assert sp.factor(derivative - 2 * expected) != 0


def test_fixed_H_reference_cancels_only_stationary_mass_piece():
    G, L, M, T = sp.symbols("G L M T")
    J = -G * T / 2
    original = (L - G * G / M) * T / 2
    counter = -G * J / M
    assert sp.factor(original + counter - L * T / 2) == 0
    assert sp.factor(original - L * T / 2) != 0


@pytest.mark.parametrize("constant,slope", ((0, 0), (3, 7), (10**800, -(10**400))))
def test_whole_outer_affine_counterterms_cancel(constant, slope):
    s = sp.symbols("s")
    local = sp.Integer(constant) + sp.Integer(slope) * s
    subtraction = local - local.subs(s, 1) - (s - 1) * sp.diff(local, s).subs(s, 1)
    assert sp.expand(subtraction) == 0


def test_finite_kinetic_coefficient_contains_Taylor_half_once():
    C, g, Q, T = sp.symbols("C g Q T")
    upper_second = (g / Q) * (16 * C / 3) * (1 / (2 * T * T))
    expected = 4 * C * g / (3 * Q * T * T)
    assert sp.factor(upper_second / 2 - expected) == 0
    assert sp.factor(upper_second - expected) != 0


def test_selected_reference_identity_with_independent_symbols():
    h, s, r, fp, r2, mF, mS, m2, RF, RS, R2 = sp.symbols(
        "h s r fp r2 mF mS m2 RF RS R2"
    )
    k = 1 + h * (r - fp) + h * h * r2
    mref = 1 + h * (mS - mF) + h * h * m2
    raw = (
        mref
        - s
        + h * (mF + (s - 1) * fp + RF - mS - (s - 1) * r - RS)
        - h * h * (m2 + (s - 1) * r2 + R2)
    )
    expected = 1 - s + (h * (RF - RS) - h * h * R2) / k
    assert sp.factor(raw / k - expected) == 0
    second = sp.diff(expected, h, 2).subs(h, 0) / 2
    assert sp.factor(second + R2 + (r - fp) * (RF - RS)) == 0
    assert sp.factor(second + R2) != 0


def test_selected_spacelike_defects_all_reduce_inverse():
    t, r, fp, r2, Ferm, Sc, Spec = sp.symbols("t r fp r2 Ferm Sc Spec")
    k = 1 + r - fp + r2
    inverse = 1 + t - (Ferm + Sc + Spec) / k
    assert sp.factor((1 - inverse / (1 + t)) - (Ferm + Sc + Spec) / (k * (1 + t))) == 0


def test_outer_OS_kernel_preserves_mass_and_residue():
    d = kernel.data()
    s = d["symbols"]["s"]
    R = d["anchored_analytic_OS_kernel"]
    assert R.subs(s, 1) == 0
    assert sp.diff(R, s).subs(s, 1) == 0
    assert sp.diff(R, s, 2).subs(s, 1) != 0


def test_actual_family_and_selected_bounds():
    d = calibration.data()
    family = d["quadratic_covariance_insertion_family_enclosure"]
    assert 0 < family["finite_positive_outer_slope_upper"] < sp.Rational(1, 10**616)
    assert (
        0
        < family["uniform_outer_OS_remainder_coefficient_upper"]
        < sp.Rational(1, 10**1017)
    )
    assert (
        0
        < d["selected_radius_two_OS_remainder_coefficient_upper"]
        < sp.Rational(2, 10**405)
    )
    assert d["selected_local_Phi_curvature_lower"] > 0
    assert d["selected_reference_window_energy_not_a_cutoff"] == 10**400
    assert d["selected_reference_window_inverse_enclosure"][
        "fractional_inverse_defect_upper"
    ] < sp.Rational(1, 10**202)


def test_zero_couplings_are_not_strict_positivity():
    assert (
        enclosure.bound(2, 0, 1, 1, 144)[
            "strict_family_slope_and_spacelike_remainder_positive"
        ]
        is False
    )
    assert (
        enclosure.bound(2, 1, 0, 1, 144)[
            "strict_family_slope_and_spacelike_remainder_positive"
        ]
        is False
    )


def test_scope_excludes_complete_two_loop_and_stable_heavy_claims():
    assert (
        "Not every enlarged-model two-loop quadratic graph" in quadratic.data()["scope"]
    )
    assert "not a declared exact stable particle" in kernel.data()["domain"]
    assert "not the complete new-model two-loop inverse" in selected.data()["scope"]
    assert audit.controls()["original_P8_not_closed"] is True


def test_exact_counts():
    assert len(audit.residuals()) == 44
    assert len(audit.gates()) == 28
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 68
