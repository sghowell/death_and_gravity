"""Independent Laurent forests, dimensional bubbles and forward coefficients."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_double_bubble_ms import (
    audit,
    bounds,
    calibration,
    conversion,
    forests,
    forward,
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


@pytest.mark.parametrize(
    "f0,f1,f2", ((2, 3, 5), (-1, 4, -3), (s.Rational(2, 3), s.Rational(5, 7), 0))
)
def test_independent_Laurent_projection_and_nested_forests(f0, f1, f2):
    e = s.symbols("e")
    P = s.Rational(2, 7) / e
    I = P + f0 + e * f1 + e**2 * f2

    def K(expr):
        expanded = s.series(expr, e, 0, 1).removeO().expand()
        return expanded - expanded.coeff(e, 0)

    raw, nested = K(I**2), K(P * I)
    assert s.expand(raw - P**2 - 2 * P * f0) == 0
    assert s.expand(nested - P**2 - P * f0) == 0
    six = I**2 - P * I - P * I - raw + nested + nested
    assert s.expand(six - (I - P) ** 2) == 0
    assert s.limit(six, e, 0) == f0**2
    # One incorrect raw-only projection leaves an uncancelled simple pole.
    wrong = I**2 - 2 * P * I - P**2 + 2 * nested
    assert s.expand(wrong).coeff(e, -1) == s.Rational(4, 7) * f0 != 0
    assert forests.poles(I**2, e) == raw


@pytest.mark.parametrize("index", range(24))
def test_each_actual_forest_has_only_its_allowed_topology(index):
    row = forests.data()["actual_MS_forest_rows"][index]
    e = forests.data()["symbols"]["epsilon"]
    terms = row["terms"]
    expected_count = {"single": 2, "disjoint": 4, "overlap": 6}[row["mode"]]
    assert len(terms) == expected_count
    result = s.expand(sum(t["term"] for t in terms))
    assert s.simplify(result - row["renormalized_form"]) == 0
    assert forests.poles(result, e) == 0
    if row["mode"] == "overlap":
        nested = [t for t in terms if len(t["forest"]) == 2]
        assert len(nested) == 2
        shared = set(nested[0]["forest"]) & set(nested[1]["forest"])
        assert len(shared) == 1
        forbidden = (set(nested[0]["forest"]) | set(nested[1]["forest"])) - shared
        assert not any(set(t["forest"]) == forbidden for t in terms)


@pytest.mark.parametrize("z", ("0", "2+.3j", "3"))
@pytest.mark.parametrize("ell", (0, 2))
def test_full_dimensional_bubble_finite_forest(z, ell):
    with mp.workdps(38):
        z = complex(z)
        z = mp.mpc(str(z.real), str(z.imag))
        IR = mp.quad(lambda x: -mp.log(1 - x * (1 - x) * z), [0, 1])
        f0 = ell + IR
        finite_square = 0
        finite_triangle = 0
        for k in range(24):
            e = mp.mpf(".04") * mp.exp(2j * mp.pi * k / 24)
            P = 1 / e
            I = (
                mp.exp(mp.euler * e + ell * e)
                * mp.gamma(1 + e)
                / e
                * mp.quad(lambda x, e=e: (1 - x * (1 - x) * z) ** (-e), [0, 1])
            )
            rawK = P * P + 2 * P * f0
            nestedK = P * P + P * f0
            six = I * I - P * I - P * I - rawK + nestedK + nestedK
            four = I * I - P * I - P * I + P * P
            assert abs(six - four) < mp.mpf("1e-30")
            finite_square += six / 24
            T = mp.mpf(".3") + e * mp.mpf(".7")
            finite_triangle += (I * T - P * T) / 24
        assert abs(finite_square - f0 * f0) < mp.mpf("1e-24")
        assert abs(finite_triangle - mp.mpf(".3") * f0) < mp.mpf("1e-24")


@pytest.mark.parametrize("M", (10, 40, 100))
def test_C_cubed_forward_coefficient_by_complex_contour(M):
    with mp.workdps(38):
        L, g = mp.mpf(2), mp.mpf(".5")
        D = M - 2
        coefficient = 0
        for k in range(32):
            nu = mp.exp(2j * mp.pi * k / 32)
            amplitude = (
                (-L + g / (D - nu)) ** 3 + (-L + g / (D + nu)) ** 3 + (-L + g / M) ** 3
            )
            coefficient += amplitude / nu**2 / 32
        tree = 2 * g / D**3
        target = tree * 3 * (-L + g / D) * (-L + 2 * g / D)
        assert abs(coefficient - target) < mp.mpf("1e-24")
        assert 0 < target / tree < 3 * L**2


@pytest.mark.parametrize(
    "r", (s.Rational(1, 100), s.Rational(1, 3), s.Rational(49, 100))
)
def test_quadratic_scale_exact_ratio_majorant(r):
    exact = 3 * (1 - r) * (1 - 2 * r)
    assert 0 < exact < 3
    assert s.expand(3 - exact - 3 * r * (3 - 2 * r)) == 0


@pytest.mark.parametrize("ell", (0, 1, 3))
def test_full_factorization_retains_both_scale_powers_and_triangles(ell):
    L, g, M, z, Q = s.symbols("L g M z Q")
    C = -L + g / (M - z)
    # Distinct analytic diagnostic factors; actual triangle bounds are tested separately.
    IR = 1 / (5 - z)
    TA = g / (M + 2 - z)
    TB = 2 * g / (M + 3 + z)
    d = s.Rational(ell) / Q
    raw = C * ((C * (IR + d) + TA) * (C * (IR + d) + TB) - TA * TB) / 4
    old = (C**3 * IR**2 + C**2 * IR * (TA + TB)) / 4
    linear = d * (2 * C**3 * IR + C**2 * (TA + TB)) / 4
    square = d * d * C**3 / 4
    assert s.factor(raw - old - linear - square) == 0
    if ell:
        assert s.diff(square, z, 2) != 0
        for missing in (TA, TB):
            assert s.diff(d * C * C * missing / 4, z, 2) != 0


@pytest.mark.parametrize("M", (40, 100, 1000))
def test_all_radius_heavy_triangle_majorant_integral(M):
    with mp.workdps(32):
        M = mp.mpf(M)
        exact = M * mp.log(4 * M) / (M - mp.mpf(".25")) ** 2 - 1 / (M - mp.mpf(".25"))
        integral = mp.quad(
            lambda y: y / ((y + mp.mpf(".25")) ** 2 * (y + M)), [0, 1, M, mp.inf]
        )
        assert abs(exact - integral) < mp.mpf("1e-25")
        assert 0 < integral < 2 * mp.log(4 * M) / M


def test_complete_actual_bound_and_matching_scope():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 10
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 3
    assert len(audit.gates()) == 23
    assert audit.controls()["rejected_inputs"] == 94
    atzero = bounds.enclosure(s.Rational(1, 100), 1, 0, 0)
    assert atzero["full_linear_scale_term_upper"] == 0
    assert atzero["quadratic_scale_b2_absolute_upper"] == 0
    assert conversion.data()["checks"]["zero_scale_change_recovers_old"] == 0
    assert forward.data()["checks"]["ordinary_disjoint_heavy_mass_product"] == 0
