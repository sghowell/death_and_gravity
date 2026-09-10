"""Independent contact integrals, bare maps and full insertion cancellation."""

from fractions import Fraction

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_finite_contact_conversion import (
    amplitude,
    audit,
    calibration,
    contact,
    conversion,
    ownership,
)


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_invalid_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("M", (s.Rational(3, 2), 2, 5, 32, 100, 1000))
def test_independent_four_dimensional_radial_integrals(M):
    with mp.workdps(45):
        m = mp.mpf(str(s.N(M, 47)))
        d = contact.data()
        for n, key in ((1, "J1"), (2, "J2")):
            radial = mp.quad(
                lambda y, n=n: y / ((y + 1) ** 2 * (y + m) ** n), [0, 1, m, mp.inf]
            )
            exact = mp.mpf(str(s.N(d[key].subs(d["symbols"]["M"], M), 47)))
            assert abs(radial - exact) < mp.mpf("1e-38") * max(1, abs(exact))
        J1 = mp.quad(lambda y: y / ((y + 1) ** 2 * (y + m)), [0, 1, m, mp.inf])
        assert 0 < J1 <= (mp.log(m + 1) + 1 / (m + 1)) / m <= mp.log(4 * m) / m


@pytest.mark.parametrize("e", ("-0.25", "0.2", "0.6", "1.25"))
@pytest.mark.parametrize("M", (2, 7, 30))
def test_dimensional_radial_integral_against_independent_Feynman_parameters(e, M):
    with mp.workdps(45):
        e, m, mu = mp.mpf(e), mp.mpf(M), mp.mpf("1.5")
        pref = mp.exp(mp.euler * e) * mu ** (2 * e)
        for n in (1, 2):
            radial = (
                pref
                / mp.gamma(2 - e)
                * mp.quad(
                    lambda y, n=n: y ** (1 - e) / ((y + 1) ** 2 * (y + m) ** n),
                    [0, 1, m, mp.inf],
                )
            )
            if n == 1:
                param = (
                    pref
                    * mp.gamma(1 + e)
                    * mp.quad(lambda x: x / (x + (1 - x) * m) ** (1 + e), [0, 1])
                )
            else:
                param = (
                    pref
                    * mp.gamma(2 + e)
                    * mp.quad(
                        lambda x: x * (1 - x) / (x + (1 - x) * m) ** (2 + e), [0, 1]
                    )
                )
            assert abs(radial - param) < mp.mpf("1e-30") * max(1, abs(param))


@pytest.mark.parametrize("M", (2, 4, 10, 100))
def test_nonzero_regulator_mass_and_quartic_derivatives(M):
    with mp.workdps(40):
        e, Q, g = mp.mpf(".3"), mp.mpf(150), mp.mpf(".5")
        pref = mp.exp(mp.euler * e) / mp.gamma(2 - e)

        def j(n, m):
            return pref * mp.quad(
                lambda y, n=n: y ** (1 - e) / ((y + 1) ** 2 * (y + m) ** n),
                [0, 1, m, mp.inf],
            )

        j1, j2 = j(1, M), j(2, M)
        assert abs(mp.diff(lambda m: j(1, m), M) + j2) < mp.mpf("1e-32")

        def sigma(L):
            return 6 * ((-L + g / M) * g * j1 + g * g * j2) / Q

        assert abs(mp.diff(sigma, 2) + 6 * g * j1 / Q) < mp.mpf("1e-32")
        assert sigma(2) < 0


@pytest.mark.parametrize(
    "L,k,c",
    (
        (2, s.Rational(1, 5), s.Rational(1, 7)),
        (3, s.Rational(1, 4), 2),
        (s.Rational(2, 5), s.Rational(1, 8), s.Rational(1, 6)),
    ),
)
def test_literal_bare_reference_expansion(L, k, c):
    h = s.symbols("h")
    sig = lambda ell: c - k * ell
    star = (L - h * c) / (1 - h * k)
    G, I0, fermion = s.Rational(3, 7), s.Rational(11, 13), s.Rational(5, 17)
    # Construct directly, independently of conversion.data().
    oldL = (
        star
        + h * (sig(star) + 3 * star**2 * I0 + fermion)
        + h * h * 6 * sig(star) * star * I0
    )
    oldG = G + h * star * G * I0 + h * h * sig(star) * G * I0
    oldM = 7 + h * G * G * I0
    for old, target in (
        (oldL, L + h * (3 * L * L * I0 + fermion)),
        (oldG, G + h * L * G * I0),
        (oldM, 7 + h * G * G * I0),
    ):
        assert s.series(old - target, h, 0, 3).removeO().expand() == 0
    for ell in (s.Rational(1, 3), 1):
        assert s.factor((star + h * sig(star) - L).subs(h, ell)) == 0


@pytest.mark.parametrize(
    "L,k,c", ((2, s.Rational(1, 5), s.Rational(1, 7)), (3, s.Rational(1, 4), 2))
)
def test_literal_complete_three_channel_amplitude(L, k, c):
    h, ell, z = s.symbols("h ell z")
    sigma = c - k * ell
    star = (L - h * c) / (1 - h * k)
    # Independent nonconstant channel functions, distinct heavy triangles.
    gs, M = s.Rational(2, 7), 11
    invariants = (z, 1 - z, s.Rational(3, 2))
    F = 0
    for j, inv in enumerate(invariants):
        C = -ell + gs / (M - inv)
        IR = (j + 1) / (3 - inv)
        TA, TB = gs / (7 - inv), gs * (j + 2) / (13 + inv)
        B = gs * gs / (17 + inv * inv)
        F += C * C * IR / 2 + C * (TA + TB) / 2 + B / 2
    tree = -ell + gs * sum(1 / (M - inv) for inv in invariants)
    old = (tree + h * (F - sigma) + h * h * sigma * s.diff(F, ell)).subs(ell, star)
    target = (tree + h * F).subs(ell, L)
    assert s.factor(s.diff(old - target, h, 2).subs(h, 0)) == 0
    assert s.factor(s.diff(old - target, h).subs(h, 0)) == 0
    # Omitting either integrated triangle leaves a nonconstant residual.
    for triangle in (gs / (7 - z), 2 * gs / (13 + z)):
        defect = -sigma * triangle / 2
        assert s.diff(defect, z, 2) != 0


@pytest.mark.parametrize("eps1", (1, 3, -2, s.Rational(7, 5)))
def test_early_regulator_truncation_leaves_finite_counterterm_defect(eps1):
    e, ell, G = s.symbols("e ell G")
    sigmaD, sigma0 = 2 + e * eps1, s.Integer(2)
    I0 = 1 / e + s.Rational(7, 3) + e * s.Rational(5, 11)
    for coeff in (6 * ell, G, 2 * G * G):
        full = sigmaD * coeff * I0 - sigmaD * coeff * I0
        wrong = sigmaD * coeff * I0 - sigma0 * coeff * I0
        assert s.expand(full) == 0
        assert s.limit(wrong, e, 0) == eps1 * coeff != 0


@pytest.mark.parametrize(
    "L", (s.Rational(1, 1000), s.Rational(1, 100), s.Rational(1, 20))
)
@pytest.mark.parametrize("h", (0, s.Rational(1, 3), 1))
def test_exact_coordinate_tail_bound_and_sign(L, h):
    bounds = calibration.enclosure(L, h)
    k = bounds["k_upper"] / 2
    sigma = -bounds["sigma_absolute_upper"] / 3
    c = sigma + k * L
    exact = (L - h * c) / (1 - h * k)
    trunc = L - h * sigma - h * h * k * sigma
    assert 0 <= exact - trunc <= bounds["coordinate_tail_after_order_two_upper"]
    assert exact >= L
    assert 1 - h * k > s.Rational(1, 2)


def test_nonconstant_sigma_requires_second_inverse_term():
    h = s.symbols("h")
    L, k, c = 2, s.Rational(1, 5), s.Rational(1, 7)
    sigma = c - k * L
    wrong = L - h * sigma
    defect = wrong + h * (c - k * wrong) - L
    assert s.expand(defect).coeff(h, 2) == k * sigma != 0


def test_bare_contact_has_nonzero_loop_insertion_b2():
    with mp.workdps(40):
        L, sigma = mp.mpf(2), mp.mpf("-.3")

        def inserted(z):
            return sigma * mp.quad(
                lambda x: -L * mp.log(1 - x * (1 - x) * z) / (16 * mp.pi**2), [0, 1]
            )

        value = mp.diff(inserted, 0, 2) / 2
        assert abs(value - sigma * L / (960 * mp.pi**2)) < mp.mpf("1e-34")
        assert value != 0


def test_quadratic_field_and_source_owners_from_explicit_functions():
    L, z, g, M, G, sigma = s.symbols("L z g M G sigma")
    f1 = L * s.Rational(3, 5) / 2 - g / (M - z)
    mass_inserted = sigma * s.diff(f1, L)
    assert s.diff(mass_inserted, z) == 0
    assert mass_inserted - mass_inserted.subs(z, 1) == 0
    source = -G * s.Rational(3, 5) / 2
    assert s.diff(source, L) == 0
    # The scalar bubble slope varies with g,M but not L.
    assert s.diff(s.diff(f1, z), L) == 0
    assert ownership.data()["checks"]["contact_zero_background_quadratic_hessian"] == 0


def test_actual_bounds_and_no_frontier_overclaim():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert len(audit.frontier()) == 9
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 1
    assert len(audit.gates()) == 23
    assert audit.controls()["rejected_inputs"] == 49
    assert d["matched_sigma_nonlocal_order_two_exact"] == 0
    assert calibration.enclosure(Fraction(1, 100), Fraction(1, 2))[
        "k_upper"
    ] == s.Rational(1, 12)
    assert conversion.data()["checks"]["same_heavy_mass_linear_reference"] == 0
    assert amplitude.data()["checks"]["order_two_full_nonlocal_cancellation"] == 0
