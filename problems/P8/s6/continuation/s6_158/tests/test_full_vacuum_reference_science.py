"""Independent Wick, regulator, simplex and complete-reference checks."""

from functools import cache
from itertools import permutations

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_full_vacuum_reference import (
    audit,
    bounds,
    calibration,
    forest,
    sector,
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


def pairings(labels):
    if not labels:
        yield ()
        return
    for j in range(1, len(labels)):
        for rest in pairings(labels[1:j] + labels[j + 1 :]):
            yield ((labels[0], labels[j]),) + rest


def test_full_species_labeled_cubic_Wick_pairings():
    vertices = (
        ("H", 0, 0),
        ("P", 0, 0),
        ("P", 0, 1),
        ("H", 1, 0),
        ("P", 1, 0),
        ("P", 1, 1),
    )
    allowed = [p for p in pairings(vertices) if all(a[0] == b[0] for a, b in p)]
    assert len(allowed) == 3
    cross = [p for p in allowed if all(a[1] != b[1] for a, b in p)]
    assert len(cross) == 2
    assert -s.Rational(1, 2) * s.Rational(1, 2) ** 2 * len(cross) == -s.Rational(1, 4)
    assert -s.Rational(1, 2) * s.Rational(1, 2) ** 2 * (
        len(allowed) - len(cross)
    ) == -s.Rational(1, 8)
    assert len(list(pairings(tuple(range(4))))) * s.Rational(1, 24) == s.Rational(1, 8)


@pytest.mark.parametrize("mu", ("2", "5", "13"))
def test_full_source_square_finite_part_by_complex_regulator_circle(mu):
    with mp.workdps(42):
        mu, G, M, Q = mp.mpf(mu), mp.mpf(".3"), mp.mpf(7), 16 * mp.pi**2
        finite = 0
        triple = 0
        for j in range(32):
            e = mp.mpf(".02") * mp.exp(2j * mp.pi * j / 32)
            Tad = mp.exp(mp.euler * e) * mu ** (2 * e) * mp.gamma(e - 1) / Q
            J = -G * Tad / 2
            square = -J * J / (2 * M)
            mass = -G * J * Tad / (2 * M)
            reducible = -G * G * Tad * Tad / (8 * M)
            finite += square / 32
            triple += (square + mass + reducible) / 32
        ell = 2 * mp.log(mu)
        expected = (
            -G * G * (2 * ell * ell + 4 * ell + 3 + mp.pi**2 / 6) / (8 * M * Q**2)
        )
        early = -G * G * (ell + 1) ** 2 / (8 * M * Q**2)
        assert abs(finite - expected) < mp.mpf("1e-34")
        assert abs(triple) < mp.mpf("1e-37")
        assert abs(finite - early) > mp.mpf("1e-7")


@pytest.mark.parametrize("values", ((2, 3, 5), (3, 7, 11)))
def test_literal_physical_OS_mass_and_slope_determine_scalar_trace(values):
    L, g, T = map(s.Integer, values)
    x, B0, B1, B2 = s.symbols("x B0 B1 B2")
    bubble = B0 + B1 * (x + 1) + B2 * (x + 1) ** 2
    inverse_loop = L * T / 2 - g * bubble
    deltaZ = -s.diff(inverse_loop, x).subs(x, -1)
    deltaM = -inverse_loop.subs(x, -1) + deltaZ
    ren = s.expand(inverse_loop + deltaZ * x + deltaM)
    assert ren.subs(x, -1) == 0
    assert s.diff(ren, x).subs(x, -1) == 0
    divided = s.cancel((deltaZ * x + deltaM) / (x + 1))
    assert s.cancel(divided - g * B1 - (-L * T / 2 + g * B0) / (x + 1)) == 0


@pytest.mark.parametrize(
    "r,t", ((s.Rational(1, 5), s.Rational(2, 3)), (s.Rational(3, 7), s.Rational(1, 4)))
)
@pytest.mark.parametrize("e", ("0.1+0.03j", "0.7-0.2j"))
def test_independent_simplex_change_of_variables_pointwise(r, t, e):
    with mp.workdps(36):
        rr, tt = mp.mpf(int(r.p)) / int(r.q), mp.mpf(int(t.p)) / int(t.q)
        ee = mp.mpc(complex(e))
        a, b, c = map(mp.mpf, (1, 2, 3))
        denom = 1 + rr + rr * tt
        xx, yy, zz = 1 / denom, rr / denom, rr * tt / denom
        U = xx * yy + xx * zz + yy * zz
        X = a * xx + b * yy + c * zz
        jac = mp.det(
            mp.matrix(
                [
                    [
                        mp.diff(lambda u: u / (1 + u + u * tt), rr),
                        mp.diff(lambda v: rr / (1 + rr + rr * v), tt),
                    ],
                    [
                        mp.diff(lambda u: u * tt / (1 + u + u * tt), rr),
                        mp.diff(lambda v: rr * v / (1 + rr + rr * v), tt),
                    ],
                ]
            )
        )
        original = U ** (ee - 2) * X ** (1 - 2 * ee) * jac
        transformed = (
            rr ** (ee - 1)
            * (1 + tt + rr * tt) ** (ee - 2)
            * (a + rr * b + rr * tt * c) ** (1 - 2 * ee)
        )
        assert abs(original - transformed) < mp.mpf("1e-31")


@cache
def gauss_nodes():
    with mp.workdps(42):
        nodes, weights = mp.gauss_quadrature(40, "legendre")
        return tuple(((x + 1) / 2, w / 2) for x, w in zip(nodes, weights))


@pytest.mark.parametrize("masses", ((1, 1, 2), (1, 2, 3)))
def test_sunset_zero_dimensional_Gaussian_normalization_independently(masses):
    # At D=0 (epsilon=2), the two momentum integrals are absent.
    # After removing the displayed MS measure factor, S=1/(a b c).
    with mp.workdps(36):
        total = 0
        for a, b, c in permutations(masses):
            for r, wr in gauss_nodes():
                for t, wt in gauss_nodes():
                    total += 2 * wr * wt * r / (a + r * b + r * t * c) ** 3
        assert abs(total - 1 / mp.mpf(masses[0] * masses[1] * masses[2])) < mp.mpf(
            "1e-24"
        )


@pytest.mark.parametrize("angle", (0, 1, 2, 3))
def test_corner_subtracted_sunset_and_Gamma_circle_majorants(angle):
    with mp.workdps(34):
        e = mp.exp(2j * mp.pi * angle / 4) / 16
        M, mu, Q = mp.mpf(7), mp.mpf(2), 16 * mp.pi**2
        A = 0
        J = (2 ** (e - 1) - 1) / (e - 1)
        for a, b, c in permutations((mp.mpf(1), mp.mpf(1), M)):
            A += a ** (1 - 2 * e) * J / e
            for r, wr in gauss_nodes():
                for t, wt in gauss_nodes():
                    F = (1 + t + r * t) ** (e - 2) * (a + r * b + r * t * c) ** (
                        1 - 2 * e
                    )
                    F0 = (1 + t) ** (e - 2) * a ** (1 - 2 * e)
                    A += wr * wt * r ** (e - 1) * (F - F0)
        assert abs(A) < 120 * (3 * M) ** (mp.mpf(9) / 8)
        gamma = mp.exp(2 * mp.euler * e) * mp.gamma(2 * e - 1)
        assert abs(gamma) < 40
        S = gamma * mu ** (4 * e) * A / Q**2
        assert abs(S) < 4800 * mu ** (mp.mpf(1) / 4) * (3 * M) ** (mp.mpf(9) / 8) / Q**2
        assert abs(mp.exp(mp.euler * e) * mp.gamma(e - 1)) < 70
        assert abs(mp.exp(mp.euler * e) * mp.gamma(e)) < 64


@pytest.mark.parametrize("mu,M", ((2, 2), (3, 5), (4, 7)))
def test_heavy_MS_mass_reference_retains_tadpole_epsilon(mu, M):
    with mp.workdps(42):
        finite = 0
        for j in range(32):
            e = mp.mpf(".02") * mp.exp(2j * mp.pi * j / 32)
            # Remove common g/(4Q^2), retain the extra pole.
            finite += (
                mp.exp(mp.euler * e)
                * mp.mpf(mu) ** (2 * e)
                * mp.gamma(e - 1)
                * mp.mpf(M) ** (1 - e)
                / e
                / 32
            )
        ell = mp.log(mp.mpf(mu) ** 2 / M)
        expected = -M * (1 + ell + ell * ell / 2 + mp.pi**2 / 12)
        assert abs(finite - expected) < mp.mpf("1e-31")


def test_explicit_power_caps_are_checked_exactly():
    out = bounds.scalar_enclosure(2, 2, 1, 1, 144, 2, 2, 2)
    assert out["complete_scalar_vacuum_absolute_upper"] == sum(
        out["group_uppers"].values()
    )
    assert all(isinstance(v, s.Rational) for v in out["group_uppers"].values())
    assert source.data()["early_finite_part_source_square_defect"] != 0
    assert sector.data()["regulator_circle_radius"] == s.Rational(1, 16)


def test_complete_vacuum_does_not_close_full_source_or_original_P8():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(forest.rows()) == 8
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 15
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 9
    pending = {r["id"]: r["status"] for r in audit.matching()}
    for name in (
        "remaining_counterterm_and_vacuum_source_assembly",
        "finite_EFT_higher_order_truncation",
        "V_contour_and_cut_control",
        "finite_gravity_G",
        "common_parent_B",
    ):
        assert pending[name] == "OPEN"
    assert len(audit.residuals()) == 50
    assert len(audit.gates()) == 25
    assert audit.controls()["rejected_inputs"] == 147
