"""Independent density, derivative, partition and spectral-transform fixtures."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_kernel_norm import audit, estimate, tail, threshold


def density_x(x):
    z = x * (x + 4) / (x + 2) ** 2
    v = mp.sqrt(z)
    C = 3 - 2 * z + 3 * z * z
    D = mp.mpf(16) / 15 + z - 3 * z * z + v * C * mp.atanh(v)
    U = v * C / 2
    return U / (D * D + mp.pi**2 * U * U)


def density_r(r):
    v = mp.tanh(r)
    A = mp.mpf(16) / 15 + v * v - 3 * v**4
    B = v * (3 - 2 * v * v + 3 * v**4)
    D = A + B * r
    return (B / 2) / (D * D + mp.pi**2 * B * B / 4)


def chi(y):
    if y <= 1:
        return mp.mpf(1)
    if y >= 2:
        return mp.mpf(0)
    t = y - 1
    return 1 - 10 * t**3 + 15 * t**4 - 6 * t**5


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_identity(name):
    assert audit.residuals()[name] == 0, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", list(audit.packets()))
def test_every_quantitative_packet(name):
    assert all(audit.packets()[name]["gates"].values())


@pytest.mark.parametrize(
    "text", ("1e-20", "1e-12", "0.0001", "0.1", "0.5", "1", "1.5", "2")
)
def test_literal_density_two_derivatives_at_threshold(text):
    with mp.workdps(100):
        x = mp.mpf(text)
        values = [abs(x**j * mp.diff(density_x, x, j)) / mp.sqrt(x) for j in range(3)]
        bounds = threshold.data()["weighted_a_envelopes_over_sqrt_x"]
        for value, bound in zip(values, bounds):
            assert value <= mp.mpf(str(bound))
        z = x * (x + 4) / (x + 2) ** 2
        C = 3 - 2 * z + 3 * z * z
        D = mp.mpf(16) / 15 + z - 3 * z * z + mp.sqrt(z) * C * mp.atanh(mp.sqrt(z))
        assert D >= mp.mpf(16) / 15


@pytest.mark.parametrize("text", ("1", "2", "10", "100", "1e6", "1e12", "1e20"))
def test_literal_density_derivatives_across_infinite_tail_coordinates(text):
    with mp.workdps(110):
        x = mp.mpf(text)
        r = mp.acosh(1 + x / 2)
        assert mp.almosteq(density_x(x), density_r(r), rel_eps=mp.mpf("1e-60"))
        bounds = tail.data()["weighted_x_derivative_envelopes_times_r_squared"]
        for j, bound in enumerate(bounds):
            actual = abs(x**j * mp.diff(density_x, x, j)) * r * r
            assert actual <= mp.mpf(str(bound))
        assert r >= mp.log(1 + x) > mp.mpf("0.5")


@pytest.mark.parametrize(
    "z", (s.Rational(1, 10000), s.Rational(1, 4), s.Rational(1, 2), s.Rational(3, 4))
)
def test_independent_threshold_integral_derivatives(z):
    with mp.workdps(60):
        zz = mp.mpf(str(z.p)) / int(z.q)
        for j, bound in enumerate((4, 16, 128)):
            direct = mp.diff(lambda q: mp.atanh(mp.sqrt(q)) / mp.sqrt(q), zz, j)
            integral = mp.quad(
                lambda t, j=j: (
                    mp.factorial(j) * t ** (2 * j) / (1 - zz * t * t) ** (j + 1)
                ),
                [0, 1],
            )
            assert mp.almosteq(direct, integral, rel_eps=mp.mpf("1e-50"))
            assert 0 < direct <= bound


@pytest.mark.parametrize("i", range(11))
def test_literal_C2_transition_derivative_envelopes(i):
    with mp.workdps(60):
        t = mp.mpf(i) / 10
        first = -30 * t * t * (1 - t) ** 2
        second = -60 * t + 180 * t * t - 120 * t**3
        assert -mp.mpf(15) / 8 <= first <= 0
        assert abs(second) < 6
        if i in (0, 10):
            assert first == 0 and second == 0


@pytest.mark.parametrize("power", (-30, -10, -1, 0, 1, 5, 10, 30))
def test_partition_telescopes_at_widely_separated_frequencies(power):
    with mp.workdps(60):
        x = mp.mpf(2) ** power * mp.mpf("1.3")
        values = [
            chi(x / mp.mpf(2) ** j) - chi(2 * x / mp.mpf(2) ** j)
            for j in range(power - 5, power + 6)
        ]
        assert all(0 <= y <= 1 for y in values)
        assert abs(sum(values) - 1) < mp.mpf("1e-55")
        assert sum(y != 0 for y in values) <= 2


@pytest.mark.parametrize("h", (s.Rational(1, 16), s.S.One, s.Integer(8)))
@pytest.mark.parametrize(
    "fraction", (s.Rational(3, 4), s.Rational(5, 4), s.Rational(7, 4))
)
def test_literal_partition_derivative_scaling(h, fraction):
    with mp.workdps(60):
        hh = mp.mpf(int(h.p)) / int(h.q)
        x = hh * mp.mpf(int(fraction.p)) / int(fraction.q)
        phi = lambda z: chi(z / hh) - chi(2 * z / hh)
        assert abs(mp.diff(phi, x)) * hh <= mp.mpf(15) / 4
        assert abs(mp.diff(phi, x, 2)) * hh**2 <= 24


@pytest.mark.parametrize("G0,G2", ((1, 1), (2, 9), (100, 1), (1, 100)))
def test_exact_half_line_Fourier_minimum_integral(G0, G2):
    t = s.symbols("positive_t", positive=True)
    T = s.sqrt(s.Rational(G2, G0))
    value = 2 * (G0 * T + s.integrate(G2 / t**2, (t, T, s.oo)))
    assert s.simplify(value - 4 * s.sqrt(G0 * G2)) == 0
    assert s.sqrt(1170) < 35


@pytest.mark.parametrize("p", (0, 1, 2, 10))
def test_same_spectral_Laplace_transform_independent_quadrature(p):
    with mp.workdps(55):
        pp = mp.mpf(p)
        transform = mp.quad(
            lambda r: (
                -2 * density_r(r) * mp.tanh(r) / (1 + (pp / (2 * mp.cosh(r))) ** 2)
            ),
            [0, 1, 3, 10, mp.inf],
        )
        H = 4 + mp.quad(
            lambda y: (
                y
                * y
                * (3 - 2 * y * y + 3 * y**4)
                * pp
                * pp
                / (4 + pp * pp * (1 - y * y))
            ),
            [0, 1],
        )
        assert abs(transform + 1 / H) < mp.mpf("1e-45")
        if p == 0:
            assert abs(transform + mp.mpf(1) / 4) < mp.mpf("1e-45")


@pytest.mark.parametrize("start", (5, 10, 20))
def test_nonzero_spectral_tail_is_not_deleted(start):
    with mp.workdps(45):
        value = mp.quad(lambda r: 2 * density_r(r) * mp.tanh(r), [start, mp.inf])
        assert value > 0
        assert value < mp.mpf(2) * int(tail.C_TAIL) / start


def test_pi_bound_has_positive_integral_proof():
    with mp.workdps(60):
        value = mp.quad(lambda t: t**4 * (1 - t) ** 4 / (1 + t * t), [0, 1])
        assert value > 0
        assert abs(value - (mp.mpf(22) / 7 - mp.pi)) < mp.mpf("1e-55")
        assert mp.pi**2 < 10


def test_signed_static_moment_is_not_an_absolute_kernel_bound():
    t = s.symbols("t", positive=True)
    f = s.exp(-t) * (1 - t)
    assert s.integrate(f, (t, 0, s.oo)) == 0
    norm = s.integrate(f, (t, 0, 1)) - s.integrate(f, (t, 1, s.oo))
    assert s.simplify(norm - 2 / s.E) == 0 and norm > 0


@pytest.mark.parametrize(
    "mass", (s.Rational(1, 1000), s.S.One, s.Integer(1000), s.Integer(10) ** 20)
)
def test_mass_scope_is_scalar_scaling_not_changed_current_parent(mass):
    assert audit.require_scope(mass) == (mass, 1, 1)
    assert estimate.data()["rounded_full_half_line_kernel_L1_upper"] == 300000000000
    assert "mass independent" in estimate.data()["mass_scope"]


def test_exact_total_and_only_one_current_parameter_replaced():
    d = estimate.data()
    assert d["unrounded_full_half_line_kernel_L1_upper"] == 35 * (
        24 * 200000000 + 48 * 20000000
    )
    assert (
        d["unrounded_full_half_line_kernel_L1_upper"]
        < d["rounded_full_half_line_kernel_L1_upper"]
    )
    assert "still not numerically bounded" in d["remaining_unevaluated_constant"]
    assert "C1<15000" in d["local_lapse_recovery"]
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 37
    assert len(audit.residuals()) == 83 and audit.scalar_entry_count() == 83
    assert len(audit.gates()) == 27 and audit.rejected_inputs() == 90
