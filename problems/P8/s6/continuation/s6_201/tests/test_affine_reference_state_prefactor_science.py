"""Independent occupation, full contact, radial and regulator controls."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_reference_state_prefactor import contact, limit, memory, prefactor


@pytest.mark.parametrize("pair", tuple(product(range(3), repeat=2)))
def test_complete_two_leg_phase_cancellation_with_independent_time_readouts(pair):
    r, c = pair
    b = s.Rational(r + 1, 16)
    d = s.Rational(c + 1, 20)
    phase1 = (3 + 4 * s.I) / 5
    phase2 = (5 - 12 * s.I) / 13
    a = s.sqrt(1 + b) * phase1
    z = s.sqrt(1 + d) * phase2
    v = s.Matrix([s.Rational(j + 1, 11) + (r - j) * s.I / 17 for j in range(10)])
    w = s.Matrix([s.Rational(j + 2, 13) + (c + j) * s.I / 19 for j in range(10)])
    x = s.Matrix([s.Rational(j + 3, 7) + (2 * r + j) * s.I / 23 for j in range(10)])
    y = s.Matrix([s.Rational(j + 1, 5) + (c - 2 * j) * s.I / 29 for j in range(10)])
    D = s.diag(1, -2, 1)
    G = s.Matrix([[2, 1, -1], [1, -1, 2], [-1, 2, -1]]) / 7
    MD = s.diag(-D, D, D, 0)
    MG = s.diag(-G, G, G, 0)
    left = ((a * v).T * MD * (z * w))[0]
    right = ((a * x).T * MG * (z * y))[0]
    unit_left = (v.T * MD * w)[0]
    unit_right = (x.T * MG * y)[0]
    actual = s.conjugate(left) * right - s.conjugate(unit_left) * unit_right
    expected = (b + d + b * d) * s.conjugate(unit_left) * unit_right
    assert s.simplify(actual - expected) == 0
    assert b * d > 0


@pytest.mark.parametrize(
    "occupation", (s.Rational(1, 100), s.Rational(1, 3), s.Rational(7, 8))
)
@pytest.mark.parametrize("fixture", range(3))
def test_literal_complete_ten_feature_contact_occupation_correction(
    occupation, fixture
):
    D = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]])
    G = s.Matrix([[2, -1, 1], [-1, -3, 2], [1, 2, 1]]) / 7
    H = (D * G + G * D) / 2
    M = s.diag(H, H, H, 0)
    v = s.Matrix([s.Rational(j + 1, 13) + (fixture + j) * s.I / 17 for j in range(10)])
    a = s.sqrt(1 + occupation) * (3 + 4 * s.I) / 5
    direct = ((a * v).conjugate().T * M * (a * v))[0] - (v.conjugate().T * M * v)[0]
    assert s.simplify(direct - occupation * (v.conjugate().T * M * v)[0]) == 0
    assert D * G != G * D and M[:3, :3] == H


@pytest.mark.parametrize("nu", (1000, 10000, 10**6))
@pytest.mark.parametrize("mu", (1000, 10000, 10**6))
def test_actual_full_occupation_increment_majorant(nu, mu):
    B = prefactor.BETA
    x = B**2 / s.Integer(nu) ** 12
    y = B**2 / s.Integer(mu) ** 12
    assert 0 < x < 1 and 0 < y < 1
    assert x + y + x * y < 2 * B**2 * (s.Integer(nu) ** -12 + s.Integer(mu) ** -12)
    source = prefactor.constants()["actual_initial_beta_envelope"]
    assert 0 < source / s.Integer(nu) ** 6 < B / s.Integer(nu) ** 6


@pytest.mark.parametrize("energy", (2, 5, 13))
@pytest.mark.parametrize("time", (0, s.Rational(1, 8), s.Rational(1, 4)))
def test_unit_WKB_Wronskian_is_restored_without_setting_physical_beta_zero(
    energy, time
):
    t = s.Symbol("t", real=True)
    W = energy + t + t * t
    phase = energy * t + t * t / 2 + t**3 / 3
    f = s.exp(-s.I * phase) / s.sqrt(2 * W)
    d = 1 + t
    p = s.diff(f, t) - d * f
    wr = s.simplify(f * s.conjugate(p) - p * s.conjugate(f))
    occupation = s.Rational(1, 7)
    alpha = s.sqrt(1 + occupation) * s.exp(s.I * s.Rational(2, 3))
    weighted = s.simplify(
        alpha * f * s.conjugate(alpha * p) - alpha * p * s.conjugate(alpha * f)
    )
    assert s.simplify(wr.subs(t, time) - s.I) == 0
    assert s.simplify(weighted.subs(t, time) - s.I * (1 + occupation)) == 0


def test_unit_Wronskian_is_not_an_exact_mode_equation():
    t = s.Symbol("t", real=True)
    W = 2 + t * t
    f = s.exp(-s.I * (2 * t + t**3 / 3)) / s.sqrt(2 * W)
    assert (
        s.simplify(
            (f * s.conjugate(s.diff(f, t)) - s.diff(f, t) * s.conjugate(f)).subs(t, 0)
        )
        == s.I
    )
    assert s.simplify((s.diff(f, t, 2) + 4 * f).subs(t, 0)) == -s.Rational(1, 4)


@pytest.mark.parametrize(
    "power,expected", ((10, 5 * s.pi / 256), (11, s.Rational(16, 315)))
)
def test_independent_complete_massive_radial_constants(power, expected):
    with mp.workdps(65):
        value = mp.quad(
            lambda theta: mp.sin(theta) ** 2 * mp.cos(theta) ** (power - 4),
            [0, mp.pi / 2],
        )
        assert abs(value - mp.mpf(str(s.N(expected, 70)))) < mp.mpf("1e-60")


def radial_tail(power, K, mass=1000):
    A = mp.mpf(25) / 16
    m = mp.mpf(mass)
    start = mp.atan(mp.mpf(K) / (A * m))
    return (
        A**3
        * m ** (3 - power)
        / (2 * mp.pi**2)
        * mp.quad(
            lambda theta: mp.sin(theta) ** 2 * mp.cos(theta) ** (power - 4),
            [start, mp.pi / 2],
        )
    )


@pytest.mark.parametrize("K", (1000, 3000, 10000, 10**6))
def test_independent_distinct_memory_and_contact_band_tails(K):
    with mp.workdps(65):
        A = mp.mpf(25) / 16
        assert radial_tail(10, mp.mpf(K) / 2) < 100 / mp.mpf(K) ** 7
        assert radial_tail(11, K) < A**11 / (144 * mp.mpf(K) ** 8)
        assert radial_tail(11, K) < radial_tail(11, mp.mpf(K) / 2)


def angular_mean(y, n, cut=None):
    if cut is None or y > cut:
        edge = mp.mpf(1)
    elif n == 0:
        return mp.mpf(0)
    elif y == 0:
        return mp.sqrt(1 + n * n) if n > cut else mp.mpf(0)
    else:
        edge = min(
            mp.mpf(1), max(mp.mpf(-1), (y * y + n * n - cut * cut) / (2 * y * n))
        )
    if edge == -1:
        return mp.mpf(0)
    upper = mp.sqrt(1 + (y + n) ** 2)
    lower = mp.sqrt(1 + y * y + n * n - 2 * y * n * edge)
    return (
        (1 + edge)
        * (upper * upper + upper * lower + lower * lower)
        / (3 * (upper + lower))
    )


def complete_pair_radial(P, K=None):
    A = mp.mpf(25) / 16
    m = mp.mpf(1000)
    n = mp.mpf(P) / (A * m)
    cut = None if K is None else mp.mpf(K) / (A * m)
    split = [mp.mpf(0), mp.pi / 2]
    if cut is not None:
        split += [mp.atan(abs(cut - n)), mp.atan(cut), mp.atan(cut + n)]
    split = sorted(set(split))

    def integrand(theta):
        y = mp.tan(theta)
        return mp.sin(theta) ** 2 * mp.cos(theta) ** 7 * angular_mean(y, n, cut)

    return 2 * A**3 / (2 * mp.pi**2 * m**7) * mp.quad(integrand, split)


@pytest.mark.parametrize("P", (0, 100, 1000, 10000, 10**6))
def test_independent_full_internal_pair_weight_at_arbitrary_transfer(P):
    with mp.workdps(65):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        J10 = 5 * A**3 / (512 * mp.pi * m**7)
        L = 1 + mp.mpf(P) / (A * m)
        actual = complete_pair_radial(P)
        assert actual <= 2 * L * J10 * (1 + mp.mpf("1e-60"))
        if P == 0:
            assert abs(actual / (2 * J10) - 1) < mp.mpf("1e-60")


@pytest.mark.parametrize("K", (1000, 10000))
@pytest.mark.parametrize("ratio", ("0", "0.25", "0.5", "1", "3"))
def test_independent_complete_removed_two_leg_union_weight(K, ratio):
    with mp.workdps(65):
        P = mp.mpf(K) * mp.mpf(ratio)
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        L = 1 + P / (A * m)
        actual = complete_pair_radial(P, K)
        assert actual > 0
        if P <= mp.mpf(K) / 2:
            assert actual <= 2 * L * radial_tail(10, mp.mpf(K) / 2)
        else:
            J10 = 5 * A**3 / (512 * mp.pi * m**7)
            assert actual / (1 + P * P) < 2 * J10 / (100 * K)


@pytest.mark.parametrize("module", (prefactor, memory, contact, limit))
def test_all_scientific_packets(module):
    d = module.data()
    for value in d["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(x) == 0 for x in entries)
    assert all(d["gates"].values())


from p8_vacuum_affine_reference_state_prefactor import audit


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_exact_residuals(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.cancel(x) == 0 for x in entries)


@pytest.mark.parametrize("case", audit.bad_cases(), ids=lambda case: case[0])
def test_all_scope_mutations_rejected(case):
    _, call, args = case
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_counts_and_guards():
    assert len(audit.residuals()) == 28 and audit.scalar_entry_count() == 135
    assert len(audit.gates()) == 36 and all(v is True for v in audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 123


def test_actual_state_and_unclosed_reference_boundary():
    assert "actual state is not changed" in audit.observable()["domain"]
    assert "not an exact bisolution" in audit.observable()["boundary"]
    assert "No momentum analyticity" in prefactor.data()["analyticity_boundary"]
    assert (
        "unknown contact-plus-five-endpoint sector is not included"
        in limit.data()["complete_known_piece"]
    )
    assert len(audit.frontier()) == 9
