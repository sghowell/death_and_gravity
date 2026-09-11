"""Independent mapped-amplitude, coefficient and ownership checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_full_two_loop_amplitude import (
    audit,
    bounds,
    calibration,
    contact,
    matching,
    ownership,
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
    "k0,k1", ((".2", ".3"), ("-.1", ".3"), (".2", "-.4"), ("0", ".3"))
)
def test_literal_full_amplitude_map_by_independent_forward_circle(k0, k1):
    with mp.workdps(40):
        L, G, M, Y, Q = map(mp.mpf, (2, 3, 40, 1, 144))
        k0, k1, t = mp.mpf(k0), mp.mpf(k1), mp.mpf(".05")
        g = G * G
        G2 = (k0 * k0 - t) * G + k1 * L * G / Q
        L2 = (3 * k0 * k0 - 2 * t) * L + 3 * k1 * L * L / Q
        M2 = k1 * g / Q

        def scalar(l, g, nu):
            return (
                l * l / (9 - nu * nu) + l * g / (16 - nu * nu) + g * g / (25 - nu * nu)
            )

        def fermion(y, nu):
            return y * y / (36 - nu * nu)

        def mapped(h, nu):
            gs = G - h * k0 * G + h * h * G2
            ls = L - 2 * h * k0 * L + h * h * L2
            ys = Y - h * k0 * Y
            ms = M + h * h * M2
            tree = -ls + gs * gs * (1 / (ms - 2 - nu) + 1 / (ms - 2 + nu) + 1 / ms)
            first = scalar(ls, gs * gs, nu) + fermion(ys, nu)
            raw = (ls**3 + ys**3 + ys * ys * gs * gs) * nu * nu / 999
            return tree + h * first + h * h * raw

        coefficient1, coefficient2 = 0, 0
        for j in range(80):
            nu = mp.exp(2j * mp.pi * j / 80)
            coefficient1 += mp.diff(lambda h, nu=nu: mapped(h, nu), 0) / (nu * nu * 80)
            coefficient2 += mp.diff(lambda h, nu=nu: mapped(h, nu), 0, 2) / (
                2 * nu * nu * 80
            )
        tree = 2 * g / (M - 2) ** 3
        scalar_b2 = L * L / 81 + L * g / 256 + g * g / 625
        fermion_b2 = Y * Y / 1296
        raw_b2 = (L**3 + Y**3 + Y * Y * g) / 999
        first_expected = scalar_b2 + fermion_b2 - 2 * k0 * tree
        second_expected = (
            raw_b2
            - k0 * (4 * scalar_b2 + 2 * fermion_b2)
            + tree * (3 * k0 * k0 - 2 * t + k1 * (2 * L - 3 * g / (M - 2)) / Q)
        )
        assert abs(coefficient1 - first_expected) < mp.mpf("1e-32")
        assert abs(coefficient2 - second_expected) < mp.mpf("1e-32")
        missing_M = second_expected + tree * 3 * k1 * g / (Q * (M - 2))
        assert abs(coefficient2 - missing_M) > mp.mpf("1e-7")
        if k0:
            double_LSZ = second_expected - 2 * k0 * first_expected
            assert abs(coefficient2 - double_LSZ) > mp.mpf("1e-4")


@pytest.mark.parametrize("nu", (".3", ".7", "1.1"))
@pytest.mark.parametrize("drift", (".2", "-.4"))
def test_pointwise_contact_cancellation_after_noncommuting_coordinate_composition(
    nu, drift
):
    with mp.workdps(40):
        nu, drift = mp.mpf(nu), mp.mpf(drift)

        def residual(h):
            L = 2 - 4 * drift * h + mp.mpf(".3") * h * h
            g = 3 - 6 * drift * h + mp.mpf(".2") * h * h
            k = g / 7
            c = g * g / 11
            star = (L - h * c) / (1 - h * k)
            sigma = -k * star + c
            F = lambda l: (
                l * l / (9 - nu * nu) + l * g / (16 - nu * nu) + g * g / (25 - nu * nu)
            )
            dF = 2 * star / (9 - nu * nu) + g / (16 - nu * nu)
            old = -star + h * (F(star) - sigma) + h * h * sigma * dF
            direct = -L + h * F(L)
            return old - direct

        for n in range(3):
            assert abs(mp.diff(residual, 0, n)) < mp.mpf("1e-35")
        assert abs(mp.diff(residual, 0, 3)) > mp.mpf("1e-6")
        assert contact.data()["nonzero_direction_commutator_example"] != 0


@pytest.mark.parametrize("weight", (2, 4))
def test_epsilon_first_coordinate_shift_only_after_complete_finite_pair(weight):
    e, k0, k1, P, F0, F1 = s.symbols("e k0 k1 P F0 F1")
    raw = P / e + F0 + e * F1
    ct = -P / e
    actual = s.expand(-weight * (k0 + e * k1) * (raw + ct)).coeff(e, 0)
    assert actual == -weight * k0 * F0
    early_ct = s.expand(-weight * (k0 + e * k1) * raw - weight * k0 * ct).coeff(e, 0)
    assert s.expand(early_ct - actual) == -weight * k1 * P


@pytest.mark.parametrize(
    "h", (s.Rational(0), s.Rational(1, 7), s.Rational(1, 2), s.Integer(1))
)
def test_formal_positive_gap_is_uniform_in_loop_homotopy(h):
    d = bounds.assemble(10, 1, 1, 1, s.Rational(1, 10))
    E2 = d["second_absolute"]
    assert h + h * h * E2 <= 1 + E2
    assert 10 - h - h * h * E2 >= d["formal_uniform_lower"] > 0
    # These finite coefficients alone cannot bound a physical cubic tail.
    if h == 1:
        assert 10 - h - h * h * E2 - 100 * h**3 < 0


@pytest.mark.parametrize("value", (True, 1.0, "1", s.oo, s.I))
def test_assembly_rejects_inexact_or_nonnumeric_bounds(value):
    with pytest.raises((TypeError, ValueError)):
        bounds.assemble(1, 1, value, 1, 1)


def test_literal_heavy_source_changes_only_vacuum_and_mass():
    Phi, H, J, G, M = s.symbols("Phi H J G M", nonzero=True)
    potential = M * H * H / 2 + H * (J + G * Phi * Phi / 2)
    stationary = s.solve(s.diff(potential, H), H)[0]
    effective = s.expand(potential.subs(H, stationary))
    assert s.diff(effective, J, Phi, Phi, Phi, Phi) == 0
    assert s.diff(effective, J, Phi, Phi) == -G / M
    assert s.diff(effective, J, J) == -1 / M


def test_correct_fundamental_G_square_and_M2_are_both_nonzero():
    d = matching.data()
    assert d["full_second_bare_scalar_map"]["M"] != 0
    assert "k1" in str(d["tree_second_relative"])
    assert matching.data()["first_regulated_direction"]["M"] == 0
    assert bounds.first_variation_upper(1, 2, 3) == 14


def test_complete_actual_enclosure_advances_only_fixed_order_Phi_item():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(d["fixed_hybrid_two_loop_group_uppers"]) == 4
    assert len(ownership.rows()) == 8
    assert len(ownership.counterterms()) == 10
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 14
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 8
    pending = {r["id"]: r["status"] for r in audit.matching()}
    assert pending[audit.TARGET] == audit.STATUS
    for name in (
        "finite_EFT_higher_order_truncation",
        "V_contour_and_cut_control",
        "finite_gravity_G",
        "common_parent_B",
        "remaining_counterterm_and_vacuum_source_assembly",
    ):
        assert pending[name] == "OPEN"
    assert len(audit.residuals()) == 41
    assert len(audit.gates()) == 25
    assert audit.controls()["rejected_inputs"] == 90
