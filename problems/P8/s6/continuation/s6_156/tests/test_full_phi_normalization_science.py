"""Independent full-regulator field/coordinate matching and local pole checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_full_phi_normalization import (
    audit,
    bounds,
    calibration,
    normalization,
    ownership,
    pole,
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
    "r,f,r1,f1",
    (
        (".2", ".1", ".3", "-.4"),
        (".1", ".3", "-.2", ".5"),
        (".2", ".2", ".3", ".1"),
    ),
)
def test_full_coordinate_shift_in_literal_field_ratio(r, f, r1, f1):
    with mp.workdps(36):
        r, f, r1, f1 = map(mp.mpf, (r, f, r1, f1))
        z, tH, z22, z21 = map(mp.mpf, ("-.7", ".08", ".3", "-.4"))
        k0 = r - f
        tMS = tH - k0 * r
        finite = 0
        for j in range(20):
            e = mp.mpf(".035") * mp.exp(2j * mp.pi * j / 20)
            kD = r - f + e * (r1 - f1)

            def residual(h, e=e, kD=kD):
                # g and Y start at one. Their regulated shifts carry
                # scalar weight two and fermionic weight one respectively.
                gstar = 1 - 2 * h * kD
                Ystar = 1 - h * kD
                kstar = gstar * (r + e * r1) - Ystar * (f + e * f1)
                ZH = (
                    1 + h * (z * Ystar / e - kstar) + h**2 * (z22 / e**2 + z21 / e - tH)
                )
                ZMS = 1 + h * z / e + h**2 * (z22 / e**2 + z21 / e)
                K = 1 + h * kD + h**2 * tMS
                return ZMS / K - ZH

            finite += mp.diff(residual, 0, 2) / 2 / 20
        assert abs(finite) < mp.mpf("1e-28")


@pytest.mark.parametrize("k1", (s.Rational(1, 3), s.Rational(-2, 7)))
def test_omitting_epsilon_coordinate_shift_leaves_finite_defect(k1):
    e, z, k0 = s.symbols("epsilon z k0")
    full = -z * (k0 + e * k1) / e
    truncated = -z * k0 / e
    assert s.expand(full - truncated).coeff(e, 0) == -z * k1
    assert s.expand(full - truncated).coeff(e, 0) != 0


@pytest.mark.parametrize("kind", ("L", "G", "M"))
@pytest.mark.parametrize("k1", (".3", "-.4"))
def test_literal_bare_vertex_match_including_first_pole_reexpression(kind, k1):
    with mp.workdps(36):
        L, G, Y, M, N, Q = map(mp.mpf, (2, 3, 1, 40, 6, 144))
        k0, k1, tMS = mp.mpf(".2"), mp.mpf(k1), mp.mpf(".05")
        w = {"L": 2, "G": 1, "M": 0}[kind]
        original = {"L": L, "G": G, "M": M}[kind]
        p = lambda l, g, y: {
            "L": (3 * l * l / 2 - 24 * N * y * y) / Q,
            "G": l * g / (2 * Q),
            "M": g * g / (2 * Q),
        }[kind]
        p0 = p(L, G, Y)
        Ep = {"L": (6 * L * L - 48 * N * Y * Y) / Q, "G": 3 * p0, "M": 2 * p0}[kind]
        comm = Ep - w * p0
        q2 = (w * (w + 1) * k0 * k0 / 2 - w * tMS) * original + k1 * comm
        double, simple = mp.mpf(".7"), mp.mpf("-.8")
        finite, pole1, pole2 = 0, 0, 0
        for j in range(20):
            e = mp.mpf(".04") * mp.exp(2j * mp.pi * j / 20)
            kD = k0 + e * k1

            def residual(h, e=e, kD=kD):
                ls = L - 2 * h * kD * L
                gs = G - h * kD * G
                ys = Y - h * kD * Y
                Xstar = original - h * w * kD * original + h * h * q2
                lhs = (
                    Xstar
                    + h * p(ls, gs, ys) / e
                    + h * h * (double / e**2 + (simple + k0 * comm) / e)
                )
                rhs = (original + h * p0 / e + h * h * (double / e**2 + simple / e)) / (
                    1 + h * kD + h * h * tMS
                ) ** w
                return rhs - lhs

            second = mp.diff(residual, 0, 2) / 2
            finite += second / 20
            pole1 += e * second / 20
            pole2 += e**2 * second / 20
        assert abs(finite) < mp.mpf("1e-27")
        assert abs(pole1) < mp.mpf("1e-29")
        assert abs(pole2) < mp.mpf("1e-30")
        if kind == "M":
            assert q2 == k1 * G**2 / Q != 0


@pytest.mark.parametrize("k1", (".3", "-.4"))
def test_full_forward_tree_map_has_G_square_and_induced_M_shift(k1):
    with mp.workdps(36):
        L, G, M, Q = map(mp.mpf, (2, 3, 40, 144))
        k0, k1, t = mp.mpf(".2"), mp.mpf(k1), mp.mpf(".05")
        G2 = (k0 * k0 - t) * G + k1 * L * G / Q
        L2 = (3 * k0 * k0 - 2 * t) * L + 3 * k1 * L * L / Q
        M2 = k1 * G * G / Q

        def forward(h, nu):
            gs = G - h * k0 * G + h * h * G2
            ls = L - 2 * h * k0 * L + h * h * L2
            ms = M + h * h * M2
            return -ls + gs * gs * (1 / (ms - 2 - nu) + 1 / (ms - 2 + nu) + 1 / ms)

        second = 0
        for j in range(32):
            nu = mp.exp(2j * mp.pi * j / 32)
            second += mp.diff(lambda h, nu=nu: forward(h, nu), 0, 2) / (2 * nu**2 * 32)
        tree = 2 * G**2 / (M - 2) ** 3
        expected = 3 * k0 * k0 - 2 * t + k1 * (2 * L - 3 * G**2 / (M - 2)) / Q
        assert abs(second / tree - expected) < mp.mpf("1e-27")
        without_M = 3 * k0 * k0 - 2 * t + 2 * k1 * L / Q
        assert abs(without_M - expected) > mp.mpf("1e-6")


@pytest.mark.parametrize("radius", (s.Rational(1, 4), s.Rational(1, 2), s.Integer(1)))
@pytest.mark.parametrize("degree", (2, 3, 8))
def test_Cauchy_extension_on_independent_holomorphic_tail_examples(radius, degree):
    with mp.workdps(32):
        rho = mp.mpf(int(radius.p)) / int(radius.q)
        E = mp.mpf(3)
        bound = mp.mpf(str(pole.tail_coefficient(3, radius).p)) / int(
            pole.tail_coefficient(3, radius).q
        )
        for j in range(12):
            v = rho * mp.exp(2j * mp.pi * j / 12)
            monomial = E * (v / 2) ** degree
            # This geometric function has supremum at most E on |v|<=2.
            geometric = (E / 2) / (1 - v / 4)
            geometric_OS = geometric - E / 2 - v * E / 8
            assert abs(monomial / v**2) <= bound + mp.mpf("1e-29")
            assert abs(geometric_OS / v**2) <= bound + mp.mpf("1e-29")


def test_second_coefficient_has_no_hidden_first_pole_product():
    d = normalization.data()
    expected = s.Symbol("t_H") - (s.Symbol("r") - s.Symbol("fp")) * s.Symbol("r")
    assert d["complete_MS_second_normalization_from_hybrid"] == expected
    assert not any(str(v) == "z11" for v in expected.free_symbols)
    assert "epsilon" in str(d["first_regulated_hybrid_parameter_direction"]["Y"])
    assert d["full_bare_second_scalar_parameter_map"]["M"] != 0


@pytest.mark.parametrize("B", (s.Rational(1, 1000), s.Rational(1, 10**18)))
def test_exact_pole_mass_residue_and_uniform_inverse_gap(B):
    v, h = s.symbols("v h")
    inverse = v + v * v * (B / 3 + h * B / 3 + h * h * B / 3)
    assert inverse.subs(v, 0) == 0
    assert s.diff(inverse, v).subs(v, 0) == 1
    assert s.limit(v / inverse, v, 0) == 1
    assert 1 - B > 0
    unprojected = inverse + s.Rational(1, 10)
    assert unprojected.subs(v, 0) != 0


def test_complete_owned_Phi_bound_is_not_amplitude_or_original_closure():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(d["complete_fixed_hybrid_two_loop_slope_group_uppers"]) == 3
    assert len(d["full_hybrid_OS_second_order_group_uppers"]) == 4
    assert len(ownership.rows()) == 3
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 14
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 7
    assert len(audit.residuals()) == 51
    assert len(audit.gates()) == 24
    assert audit.controls()["rejected_inputs"] == 96
    assert bounds.normalization_upper(1, 2, 3) == 7
    pending = {r["id"]: r["status"] for r in audit.matching()}
    assert pending["full_matched_two_loop_pole_and_amplitude_error"] == "OPEN"
    assert pending["common_parent_B"] == "OPEN"
