"""Independent matrix, dimensional-integral and full-potential diagnostics."""

import mpmath as mp
import pytest
import sympy as s
from p8_affine import verify as exact
from p8_vacuum_flat_local_energy import (
    audit,
    energy,
    extension,
    frames,
    reference,
)
from p8_vacuum_superadiabatic_state_energy import audit as previous


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_invalid_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "angles",
    [
        ("0.1", "0.02", "-0.003", "0.0004"),
        ("-0.2", "0.1", "0.05", "-0.01"),
        ("0.003", "-0.001", "0.00002", "0.000004"),
    ],
)
def test_full_Pauli_product_against_Bloch_formula(angles):
    with mp.workdps(70):
        a, b, c, d = map(mp.mpf, angles)
        sigma1 = mp.matrix([[0, 1], [1, 0]])
        sigma2 = mp.matrix([[0, -1j], [1j, 0]])
        sigma3 = mp.matrix([[1, 0], [0, -1]])
        ux = lambda t: mp.cos(t / 2) * mp.eye(2) + 1j * mp.sin(t / 2) * sigma1
        uy = lambda t: mp.cos(t / 2) * mp.eye(2) - 1j * mp.sin(t / 2) * sigma2
        U = ux(a) * uy(b) * ux(c) * uy(d)
        C = U * mp.matrix([[0, 0], [0, 1]]) * U.H
        observed = -sum((sigma3 * C)[j, j] for j in range(2))
        expected = (
            mp.cos(a) * mp.cos(b) * mp.cos(c) * mp.cos(d)
            - mp.cos(a) * mp.sin(b) * mp.sin(d)
            - mp.sin(a) * mp.sin(c) * mp.cos(d)
        )
        assert abs(observed - expected) < mp.mpf("1e-65")
        bound = abs(a * c) + b * b + c * c + d * d
        assert abs(observed - mp.cos(a)) <= bound
        diagonal_product = mp.cos(a) * mp.cos(b) * mp.cos(c) * mp.cos(d)
        assert abs(observed - diagonal_product) > mp.mpf("1e-12")


@pytest.mark.parametrize("ratio", ["0.001", "0.01", "0.1"])
@pytest.mark.parametrize("scale", ["0.001", "0.05", "0.2"])
def test_Bloch_angle_remainder_without_erasing_mixed_terms(ratio, scale):
    with mp.workdps(65):
        r, A = mp.mpf(ratio), mp.mpf(scale)
        for signs in [(1, 1, 1, 1), (1, -1, 1, -1), (-1, 1, 1, -1)]:
            a, b, c, d = [signs[j] * A * r**j for j in range(4)]
            z = (
                mp.cos(a) * mp.cos(b) * mp.cos(c) * mp.cos(d)
                - mp.cos(a) * mp.sin(b) * mp.sin(d)
                - mp.sin(a) * mp.sin(c) * mp.cos(d)
            )
            assert abs(z - mp.cos(a)) < 3 * A * A * r * r


@pytest.mark.parametrize("q", ["0", "0.00001", "0.01", "0.3", "1", "3"])
def test_full_first_frame_Taylor_integral_remainder(q):
    with mp.workdps(65):
        q = mp.mpf(q)
        u = q * q
        residual = 1 / mp.sqrt(1 + u) - 1 + u / 2
        integral = (
            mp.quad(
                lambda t: (1 - t) * mp.mpf(3) / 4 / (1 + t * u) ** mp.mpf("2.5"), [0, 1]
            )
            * u
            * u
        )
        assert abs(residual - integral) < mp.mpf("1e-58")
        assert 0 <= residual <= mp.mpf(3) / 8 * u * u


@pytest.mark.parametrize("mass", ["0.5", "1", "3", "10"])
def test_complete_radial_factors_independent_quadrature(mass):
    with mp.workdps(65):
        m = mp.mpf(mass)
        cases = [
            (lambda p: p**3 / (p * p + m * m) ** mp.mpf("2.5"), mp.mpf(2) / (3 * m)),
            (lambda p: p**4 / (p * p + m * m) ** mp.mpf("3.5"), 1 / (5 * m * m)),
            (
                lambda p: p**6 / (p * p + m * m) ** mp.mpf("5.5"),
                mp.mpf(2) / (63 * m**4),
            ),
        ]
        for f, target in cases:
            assert abs(mp.quad(f, [0, m, mp.inf]) - target) < mp.mpf("1e-56")


@pytest.mark.parametrize("eps", ["0.125", "0.25", "0.375"])
def test_dimensionally_continued_derivative_integral_directly(eps):
    with mp.workdps(65):
        e = mp.mpf(eps)
        dim = 3 - 2 * e
        m = mp.mpf("1.7")

        def radial(u):
            if u >= 0:
                return mp.exp(-2 * e * u) / (1 + mp.exp(-2 * u)) ** mp.mpf("2.5")
            return mp.exp((dim + 2) * u) / (1 + mp.exp(2 * u)) ** mp.mpf("2.5")

        measure = 2 / (4 * mp.pi) ** (dim / 2) / mp.gamma(dim / 2)
        numeric = (
            measure
            * m ** (-2 * e)
            * mp.quad(radial, [-mp.inf, 0, 1, 5, 50, 500, mp.inf])
        )
        formula = (
            m ** (-2 * e)
            / (4 * mp.pi) ** (dim / 2)
            * (dim / 2)
            * mp.gamma(e)
            / mp.gamma(mp.mpf("2.5"))
        )
        assert abs(numeric - formula) < mp.mpf("1e-55")


@pytest.mark.parametrize("mass_ratio", ["0.997", "1", "1.003", "2"])
def test_dimensional_finite_parts_keep_evanescent_angular_factor(mass_ratio):
    with mp.workdps(80):
        ell = 2 * mp.log(mp.mpf(mass_ratio))
        e = mp.mpf("1e-25")
        common = mp.exp((mp.euler - ell) * e) * mp.gamma(1 + e)
        zero = 2 * common / ((e - 1) * (e - 2))
        second = (1 - 2 * e / 3) * common
        finite0 = (zero - 1) / e
        finite2 = (second - 1) / e
        assert abs(finite0 - (mp.mpf("1.5") - ell)) < mp.mpf("1e-22")
        assert abs(finite2 - (-mp.mpf(2) / 3 - ell)) < mp.mpf("1e-22")
        wrong = (common - 1) / e
        assert abs(wrong - finite2) > mp.mpf("0.6")


@pytest.mark.parametrize("mass", ["36", "50", "200"])
@pytest.mark.parametrize("ratio", ["-0.5", "-0.003", "0.003", "0.2", "0.5"])
def test_full_potential_and_fixed_anchor_against_independent_integrals(mass, ratio):
    with mp.workdps(90):
        m = mp.mpf(mass)
        r = mp.mpf(ratio)
        N = 6
        Q = 16 * mp.pi**2
        y = mp.mpf("0.02")
        chi = r * m / y
        J = mp.quad(lambda x: -mp.log1p(-x * (1 - x) / (m * m)), [0, mp.mpf("0.5"), 1])
        f0 = 4 * N * y * y * m * m / Q
        f1 = 2 * N * y * y / Q * (2 * m * m + (4 * m * m - 1) * J)
        increment = f1 - f0
        assert 0 < increment < 4 * N * y * y / (3 * Q)
        F = lambda u: u**4 * (2 * mp.log(u) - mp.mpf("1.5"))
        direct = -N * m**4 / (2 * Q) * (F(1 + r) + F(1 - r) + 3) - f1 * chi * chi / 2
        low = -increment * chi * chi / 2 - mp.mpf(8) * N * y**4 * chi**4 / (3 * Q)
        high = direct - low
        bound = N * m**4 * r**6 / (15 * Q * (1 - r * r))
        assert 0 < high <= bound
        series = (
            mp.fsum(
                48 * N * mp.mpf(r) ** k / (k * (k - 1) * (k - 2) * (k - 3) * (k - 4))
                for k in range(6, 242, 2)
            )
            * m**4
            / Q
        )
        assert abs(high - series) < mp.mpf("1e-65") * max(1, abs(high))
        assert (
            abs(direct)
            <= increment * chi * chi / 2
            + mp.mpf(8) * N * y**4 * chi**4 / (3 * Q)
            + bound
        )


def test_actual_huge_mass_anchor_cancellation_is_resolved():
    with mp.workdps(600):
        m = mp.mpf(10) ** 200
        y = 3 * mp.mpf(10) ** -103
        chi = mp.mpf(10) ** 300
        N = 6
        Q = 16 * mp.pi**2
        r = y * chi / m
        J = mp.quad(lambda x: -mp.log1p(-x * (1 - x) / (m * m)), [0, mp.mpf("0.5"), 1])
        f0 = 4 * N * y * y * m * m / Q
        f1 = 2 * N * y * y / Q * (2 * m * m + (4 * m * m - 1) * J)
        normalized = (f1 - f0) / (2 * N * y * y / Q)
        assert normalized > 0
        assert abs(normalized - mp.mpf(2) / 3) < mp.mpf("1e-180")
        F = lambda u: u**4 * (2 * mp.log(u) - mp.mpf("1.5"))
        direct = -N * m**4 / (2 * Q) * (F(1 + r) + F(1 - r) + 3) - f1 * chi * chi / 2
        quartic = -mp.mpf(8) * N * y**4 * chi**4 / (3 * Q)
        assert abs(direct - quartic) / abs(quartic) < mp.mpf("1e-6")
        assert abs(direct) < mp.mpf(10) ** 789


def test_complete_potential_tail_coefficients_are_positive_and_decrease():
    for n in range(6, 102, 2):
        a = s.Rational(48 * 6, n * (n - 1) * (n - 2) * (n - 3) * (n - 4))
        nxt = s.Rational(48 * 6, (n + 2) * (n + 1) * n * (n - 1) * (n - 2))
        assert 0 < nxt < a <= s.Rational(6, 15)


@pytest.mark.parametrize("time", ["-5", "-1", "0", "0.5", "2", "20"])
def test_both_full_profile_log_and_derivative_bounds(time):
    with mp.workdps(70):
        x = mp.mpf(time)
        m = mp.mpf(10) ** 200
        delta = 3 * mp.mpf(10) ** 197
        tau = mp.mpf(10) ** -100
        profile = x / (1 + x**8) ** (mp.mpf(1) / 8)
        derivative = (1 + x**8) ** (-mp.mpf(9) / 8)
        for sign in (-1, 1):
            M = m + sign * delta * profile
            Md = sign * delta / tau * derivative
            assert abs(mp.log(M * M / (m * m))) < mp.mpf("0.01")
            assert abs(Md) <= delta / tau
            assert abs(
                Md * Md * (mp.log(M * M / (m * m)) + mp.mpf(2) / 3)
            ) <= delta * delta / (tau * tau)


def test_halfline_unitary_bound_does_not_replace_projector_by_exact_state():
    with mp.workdps(70):
        s2 = mp.matrix([[0, -1j], [1j, 0]])
        s3 = mp.matrix([[1, 0], [0, -1]])
        U = mp.eye(2)
        bound = mp.mpf(0)
        for e, g, t in [
            ("3", "0.003", "0.4"),
            ("2", "-0.007", "0.6"),
            ("4", "0.002", "0.2"),
        ]:
            e, g, t = map(mp.mpf, (e, g, t))
            w = mp.sqrt(e * e + g * g)
            H = e * s3 + g * s2
            step = mp.cos(w * t) * mp.eye(2) - 1j * mp.sin(w * t) / w * H
            U = step * U
            bound += abs(g) * t
        assert 0 < abs(U[0, 1]) <= bound
        projector = mp.matrix([[0, 0], [0, 1]])
        evolved = U * projector * U.H
        assert mp.norm(evolved - projector) > 0
        assert mp.norm(U.H * U - mp.eye(2)) < mp.mpf("1e-60")


def test_curvature_improvement_zero_energy_not_zero_pressure():
    t = s.Symbol("t", real=True)
    F = s.Function("F")(t)
    box = s.diff(F, t, 2)
    energy_component = box - s.diff(F, t, 2)
    pressure_component = -box
    assert energy_component == 0
    assert pressure_component.subs(F, t * t).doit() == -2


def test_named_mass_reference_extension_and_original_unsaturated_control():
    phi, R = s.symbols("Phi R", positive=True)
    chi2 = phi * phi / (1 + (phi / R) ** 8) ** s.Rational(1, 4)
    f1 = s.Symbol("fixed_positive_mass_anchor", positive=True)
    assert s.limit(chi2, phi, s.oo) == R * R
    assert s.limit(-f1 * phi * phi / 2, phi, s.oo) == -s.oo
    assert s.limit(-f1 * chi2 / 2, phi, s.oo) == -f1 * R * R / 2
    assert extension.loop_floor(4) == 4
    assert extension.loop_floor(2) == 5
    assert extension.loop_floor(0) == 6
    assert extension.loop_floor(10) == 1
    assert "not fitted" in extension.data()["finite_functional_definition"]


def test_zero_variation_and_actual_complete_scope():
    m = s.Integer(10) ** 200
    tau = s.Rational(1, 10**100)
    assert all(v == 0 for v in energy.enclosures(m, 0, tau).values())
    actual = energy.data()
    bounds = actual["actual_uniform_exact_rational_enclosures"]
    assert 0 < bounds["complete_absolute_one_loop_energy_upper"] < 10**789
    assert 0 < actual["complete_exact_state_subtracted_remainder_upper"] < 10**411
    assert actual["absolute_one_loop_energy_over_named_kappa"] < s.Rational(1, 10**11)
    assert "not silently the entire" in reference.data()["reference_coordinate"]
    assert "not a relative error" in actual["reference_ratio_warning"]
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert "NOT_FULL_CURVED_INTERACTING_PARENT_STRESS" in audit.matching()[-1]["status"]


def test_evolution_generator_connection_is_not_the_physical_energy():
    with mp.workdps(70):
        omega = mp.mpf(3)
        for q in map(mp.mpf, ("0.01", "0.1", "0.3")):
            physical = -omega / mp.sqrt(1 + q * q)
            generator = -omega * mp.sqrt(1 + q * q)
            connection = -omega * q * q / mp.sqrt(1 + q * q)
            assert abs(generator - physical - connection) < mp.mpf("1e-65")
            assert abs(generator - physical) > 0
    assert "not to that physical energy" in frames.data()["physical_energy_operator"]


def test_all_own_report_payloads_use_exact_serialization():
    for mod in audit.MODULES:
        exact.serialize({k: v for k, v in mod.data().items() if k != "checks"})
    assert all(
        type(v) is str for v in energy.data()["decimal_diagnostics_only"].values()
    )
    assert all(
        v.is_Rational
        for v in energy.data()["actual_uniform_exact_rational_enclosures"].values()
    )
    with pytest.raises(ValueError, match="Inexact or nonfinite"):
        exact.serialize(s.Float(1))
