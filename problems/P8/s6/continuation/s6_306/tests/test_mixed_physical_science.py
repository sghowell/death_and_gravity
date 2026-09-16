"""Independent projective, physical-sheet, finite-term and all-angle checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_mixed_gravity_physical_rate import (
    assembly,
    audit,
    bounds,
    masters,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_all_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def lower_log(value):
    if mp.im(value) == 0 and mp.re(value) < 0:
        return mp.log(-mp.re(value)) - mp.j * mp.pi
    return mp.log(value)


def independent_light_integrals(channel, heavy):
    channel, heavy = map(mp.mpf, (channel, heavy))

    def arguments(x):
        if channel <= 4:
            return 1 - channel * x * (1 - x), mp.mpf(1)
        f = x * (1 - x) * (1 - 2 * x)
        den = (
            1
            - channel * (x * (1 - x) + f * f)
            - mp.j * channel * x * (1 - x) * (1 - 2 * x) ** 2
        )
        return den, 1 + mp.j * (1 - 6 * x + 6 * x * x)

    def radial(den):
        rr = den / heavy
        al = (1 + mp.sqrt(1 - 4 * rr)) / 2
        be = rr / al
        return (-lower_log(be) / al + mp.log(al) / be) / (heavy * (al - be))

    def integral(fun):
        return mp.quad(
            lambda x: fun(*arguments(x)),
            [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1],
        )

    return (
        integral(lambda d, v: v / d),
        integral(lambda d, v: v * lower_log(d) / d),
        integral(lambda d, v: v * lower_log(d)),
        -integral(lambda d, v: v * radial(d)),
    )


def independent_Box16_finite_remainder(channel, partner, heavy, j1):
    channel, partner, heavy = map(mp.mpf, (channel, partner, heavy))
    beta = mp.sqrt(1 - 4 / partner)
    x = (beta - 1) / (beta + 1)
    h = 2 / (mp.sqrt(heavy) + mp.sqrt(heavy - 4))
    bracket = (
        2 * mp.log(x) * (mp.log((heavy - channel) / mp.sqrt(heavy)) + mp.log(1 - x * x))
        + 2 * mp.log(h) ** 2
        + mp.polylog(2, x * x)
        + mp.polylog(2, 1 - x * h * h)
        + 2 * mp.polylog(2, 1 - x)
        + mp.polylog(2, 1 - x / (h * h))
        - mp.pi**2 / 6
    )
    return -x * bracket / (1 - x * x) - j1 / 2


@pytest.mark.parametrize(
    "partner,heavy", ((-3, 128), (2, 128), (mp.mpf("6.25"), 128), (9, 1024), (16, 128))
)
def test_whole_Box16_normalization_and_physical_sheet(partner, heavy):
    with mp.workdps(45):
        j0, j1, ll, tri = independent_light_integrals(partner, heavy)
        for channel in (-7, 3):
            actual = independent_Box16_finite_remainder(channel, partner, heavy, j1)
            expected = tri + j0 * mp.log(1 - mp.mpf(channel) / heavy)
            assert abs(actual - expected) < mp.mpf("1e-35")
            assert abs(actual + j1 / 2 - expected) > mp.mpf("1e-3")
        beta = mp.sqrt(1 - 4 / mp.mpf(partner))
        x = (beta - 1) / (beta + 1)
        assert abs(j0 + 2 * x * mp.log(x) / (1 - x * x)) < mp.mpf("1e-35")
        assert abs(j0) < 16 and abs(j1) < 128 and abs(ll) < 16
        assert abs(tri) < 4 * (mp.log(heavy) + 9) / heavy
        if partner > 4:
            assert mp.im(j0) > 0 and mp.im(tri) < 0
            assert abs(j0 - mp.conj(j0)) > mp.mpf(".1")


@pytest.mark.parametrize("channel,den,heavy", ((2, "1.5", 16), (-3, ".75", 32)))
def test_literal_original_two_parameter_K_independent_of_closed_formula(
    channel, den, heavy
):
    with mp.workdps(45):
        den = mp.mpf(den)

        def inner(u):
            if u == 0 or u == 1:
                return mp.mpf(0)
            maximum = 1 - u
            corner = u * u * den / heavy
            knots = sorted(
                {mp.mpf(0), min(corner, maximum), min(10 * corner, maximum), maximum}
            )
            return mp.quad(
                lambda h: (
                    u
                    * (u + 2 * h)
                    / (u * u * den + h * (heavy - channel * (1 - u - h))) ** 2
                ),
                knots,
            )

        value = mp.quad(inner, [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1])
        expected = -mp.log(1 - mp.mpf(channel) / heavy) / (channel * den)
        assert abs(value - expected) < mp.mpf("1e-25")


@pytest.mark.parametrize("order", (1, 2, 3))
@pytest.mark.parametrize("den", (2, 3 - 2j, -1 - 2j))
def test_full_endpoint_radial_moments_against_original_u_integral(order, den):
    with mp.workdps(60):
        den = mp.mpc(den)
        fn = s.lambdify(
            (masters.A, masters.N),
            masters.radial_moment(order, masters.A, masters.N),
            "mpmath",
        )
        exact = fn(den, 128)
        numeric = mp.quad(
            lambda u: u**order / (128 * (1 - u) + den * u * u),
            [0, mp.mpf(".5"), mp.mpf(".9"), 1],
        )
        assert abs(exact - numeric) < mp.mpf("1e-45")


@pytest.mark.parametrize("order", (0, 1, 2))
def test_Q_removable_zero(order):
    assert masters.Q_moment(order, 0) == s.Rational(1, order + 1)


@pytest.mark.parametrize("channel", (0, s.Rational(1, 10**50), -12, 16))
def test_exact_box_small_channel_and_heavy_gap(channel):
    n = s.Integer(128)
    if channel == 0:
        assert masters.box_multiplier(channel, n) == 1 / n
    else:
        assert masters.box_multiplier(channel, n) > 0


@pytest.mark.parametrize(
    "energy,transfer",
    (
        (s.Rational(25, 4), -s.Rational(9, 8)),
        (9, -1),
        (16, -6),
        (9, -s.Rational(1, 10**400)),
        (9, -5 + s.Rational(1, 10**400)),
    ),
)
@pytest.mark.parametrize("reference", (s.Rational(1, 4), s.Rational(1, 2), 1))
def test_original_uniform_full_Born_rate_domain(energy, transfer, reference):
    assert bounds.original_known_rate_bound(energy, transfer, reference) == s.Rational(
        4, 10**204
    )
    cosine = 1 + 2 * transfer / (energy - 4)
    am = source.born.matter_born(energy, cosine, source.HEAVY_MASS2, source.CUBIC)
    ag = source.born.gravity_born(energy, cosine, source.KAPPA)
    assert am > 4 * source.CUBIC**2 / source.HEAVY_MASS2**3 and ag > 0


@pytest.mark.parametrize(
    "value", (True, False, 1.0, s.Float(1), "1", None, s.I, s.oo, s.nan, s.Symbol("x"))
)
def test_inexact_or_unsupported_physical_coordinates_rejected(value):
    with pytest.raises((ValueError, TypeError)):
        bounds.require_physical(9, value)


@pytest.mark.parametrize(
    "energy,transfer", ((4, -1), (17, -1), (9, 0), (9, -5), (9, 1))
)
def test_physical_endpoints_and_wrong_compact_window_rejected(energy, transfer):
    with pytest.raises(ValueError):
        bounds.require_physical(energy, transfer)


@pytest.mark.parametrize("reference", (0, s.Rational(1, 8), 2))
def test_auxiliary_reference_window_not_silently_enlarged(reference):
    with pytest.raises(ValueError):
        bounds.original_known_rate_bound(9, -1, reference)


def test_general_bound_requires_the_quartic_premise():
    with pytest.raises(ValueError):
        bounds.general_amplitude_majorant(128, 1, 1, 1000, 1)
    with pytest.raises(ValueError):
        bounds.general_amplitude_majorant(64, 1, 0, 1000, 1)
    assert bounds.general_amplitude_majorant(128, 1, 0, 1000, 1) > 0


def test_no_default_physical_matching_is_assigned():
    with pytest.raises(TypeError):
        assembly.minimal_unassigned_matching(9, -1, source.HEAVY_MASS2, source.KAPPA)
    value = assembly.minimal_unassigned_matching(
        9,
        -1,
        source.HEAVY_MASS2,
        source.KAPPA,
        assembly.ELL,
        assembly.RH,
        assembly.LOCAL,
    )
    assert all(
        s.diff(value, x) != 0 for x in (assembly.ELL, assembly.RH, assembly.LOCAL)
    )


def test_uniform_known_bound_not_higher_loop_or_Regge_claim():
    data = bounds.data()
    assert all(x > 0 for x in data["whole_positive_arithmetic_margins"].values())
    assert data["whole_forward_normalized_known_limit"] == 0
    assert "loop squares" in data["whole_original_rate_proof"]
    assert "Regge" in data["whole_endpoint_boundary"]


def test_all_original_parameters_and_frontiers_unchanged():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.frontier() == audit.previous.frontier()
    assert audit.qualifications() == audit.previous.qualifications()
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
