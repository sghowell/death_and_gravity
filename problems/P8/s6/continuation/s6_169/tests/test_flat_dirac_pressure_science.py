"""Independent pressure, dimensional-regulator and projector diagnostics."""

import runpy
from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_affine import verify as exact
from p8_vacuum_flat_dirac_pressure import audit, pressure
from p8_vacuum_flat_local_energy import audit as previous


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_each_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_bad_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@cache
def independent_frames():
    # Frozen parent's independent full-function Taylor-jet evaluator.
    return runpy.run_path(
        "problems/P8/s6/continuation/s6_167/tests/test_superadiabatic_state_energy_science.py"
    )["frames"]


@pytest.mark.parametrize("time", ["-2", "0.5", "2"])
@pytest.mark.parametrize("momentum", ["0.25", "1", "3"])
def test_full_four_frame_physical_pressure_matrix(time, momentum):
    with mp.workdps(80):
        x, p = mp.mpf(time), mp.mpf(momentum)
        m, delta, tau = mp.mpf(2), mp.mpf("0.01"), mp.mpf(3)
        M = m + delta * x / (1 + x**8) ** (mp.mpf(1) / 8)
        omega = mp.sqrt(p * p + M * M)
        s1 = mp.matrix([[0, 1], [1, 0]])
        s2 = mp.matrix([[0, -1j], [1j, 0]])
        s3 = mp.matrix([[1, 0], [0, -1]])
        ux = lambda t: mp.cos(t / 2) * mp.eye(2) + 1j * mp.sin(t / 2) * s1
        uy = lambda t: mp.cos(t / 2) * mp.eye(2) - 1j * mp.sin(t / 2) * s2
        chain = independent_frames()(x, p, m, delta, tau)
        angles = [mp.atan(g / e) for e, g in chain[:4]]
        U = ux(angles[0]) * uy(angles[1]) * ux(angles[2]) * uy(angles[3])
        P = mp.matrix([[0, 0], [0, 1]])
        mass_rotation = uy(mp.atan(p / M))
        Cphysical = mass_rotation * U * P * U.H * mass_rotation.H
        observed = sum((p * s1 * Cphysical / 3)[j, j] for j in range(2))
        a, b, c, d = angles
        nx = mp.cos(b) * mp.sin(d) + mp.sin(b) * mp.cos(c) * mp.cos(d)
        nz = (
            mp.cos(a) * mp.cos(b) * mp.cos(c) * mp.cos(d)
            - mp.cos(a) * mp.sin(b) * mp.sin(d)
            - mp.sin(a) * mp.sin(c) * mp.cos(d)
        )
        expected = -p * p * nz / (3 * omega) - p * M * nx / (3 * omega)
        assert abs(observed - expected) < mp.mpf("1e-70")
        assert abs(observed + p * p * nz / (3 * omega)) > mp.mpf("1e-15")
        Pi = p * p * s3 / (3 * omega) + p * M * s1 / (3 * omega)
        assert mp.norm(Pi * Pi - p * p * mp.eye(2) / 9) < mp.mpf("1e-70")
        # These moderate parameters test identities, not the large-gap enclosure.


@pytest.mark.parametrize("A", ["0.001", "0.01", "0.1"])
@pytest.mark.parametrize("ratio", ["0.001", "0.01", "0.1"])
def test_complete_x_remainder_all_angle_signs(A, ratio):
    with mp.workdps(80):
        A, r = mp.mpf(A), mp.mpf(ratio)
        for bits in range(16):
            signs = [1 if (bits >> j) & 1 else -1 for j in range(4)]
            q0 = signs[0] * A / 4
            q1 = signs[1] * A * r
            b = mp.atan(q1)
            c = mp.atan(signs[2] * A * r * r)
            d = mp.atan(signs[3] * A * r**3)
            nx = mp.cos(b) * mp.sin(d) + mp.sin(b) * mp.cos(c) * mp.cos(d)
            leading = q1 * (1 + q0 * q0) ** mp.mpf("1.5")
            bound = A * r**3 + 3 * A**3 * r
            assert abs(nx - leading) <= bound
            assert abs(leading - q1) <= 2 * abs(q1) * q0 * q0
            assert abs(mp.atan(q1) - q1) <= abs(q1) ** 3 / 3
            assert abs(mp.sin(b) - b) <= abs(b) ** 3 / 6


def test_tiny_angle_remainder_is_resolved_at_actual_parameter_scales():
    with mp.workdps(600):
        A = mp.mpf("1e-103")
        r = mp.mpf("1e-98")
        q0 = A / 4
        q1 = A * r
        b = mp.atan(q1)
        c = mp.atan(A * r * r)
        d = mp.atan(A * r**3)
        nx = mp.cos(b) * mp.sin(d) + mp.sin(b) * mp.cos(c) * mp.cos(d)
        leading = q1 * (1 + q0 * q0) ** mp.mpf("1.5")
        difference = nx - leading
        assert difference != 0
        assert abs(difference) <= A * r**3 + 3 * A**3 * r


@pytest.mark.parametrize("mass", ["0.5", "1", "3", "10"])
def test_full_off_diagonal_radial_integrals(mass):
    with mp.workdps(70):
        m = mp.mpf(mass)
        a = mp.quad(lambda p: p**3 / (p * p + m * m) ** 3, [0, m, mp.inf])
        b = mp.quad(lambda p: p**5 / (p * p + m * m) ** 5, [0, m, mp.inf])
        assert abs(a - 1 / (4 * m * m)) < mp.mpf("1e-60")
        assert abs(b - 1 / (24 * m**4)) < mp.mpf("1e-60")


@pytest.mark.parametrize("eps", ["0.125", "0.25", "0.375"])
def test_dimensional_J1_radial_integral_with_full_measure(eps):
    with mp.workdps(70):
        e = mp.mpf(eps)
        dim = 3 - 2 * e
        m = mp.mpf("1.7")

        def integrand(u):
            if u >= 0:
                return mp.exp(-2 * e * u) / (1 + mp.exp(-2 * u)) ** mp.mpf("3.5")
            return mp.exp((dim + 4) * u) / (1 + mp.exp(2 * u)) ** mp.mpf("3.5")

        measure = 2 / (4 * mp.pi) ** (dim / 2) / mp.gamma(dim / 2)
        direct = (
            measure
            * m ** (-2 * e)
            * mp.quad(integrand, [-mp.inf, 0, 1, 5, 50, 500, mp.inf])
        )
        target = (
            m ** (-2 * e)
            / (4 * mp.pi) ** (dim / 2)
            * (dim / 2)
            * (dim / 2 + 1)
            * mp.gamma(e)
            / mp.gamma(mp.mpf("3.5"))
        )
        assert abs(direct - target) < mp.mpf("1e-60")
        Q = 16 * mp.pi**2
        J2 = mp.quad(
            lambda p: p**4 / (p * p + m * m) ** mp.mpf("3.5"), [0, m, mp.inf]
        ) / (2 * mp.pi**2)
        assert abs(J2 - 8 / (5 * Q * m * m)) < mp.mpf("1e-60")


@pytest.mark.parametrize("mass_ratio", ["0.997", "1", "1.003", "2"])
def test_isotropic_regulator_factor_changes_pressure_finite_terms(mass_ratio):
    with mp.workdps(90):
        ell = 2 * mp.log(mp.mpf(mass_ratio))
        e = mp.mpf("1e-25")
        common = mp.exp((mp.euler - ell) * e) * mp.gamma(1 + e)
        J0 = (1 - 2 * e / 3) * common
        J1 = (1 - 2 * e / 3) * (1 - 2 * e / 5) * common
        first = (-2 * J0 / (3 - 2 * e) + mp.mpf(2) / 3) / e
        second = (J1 / (3 - 2 * e) - mp.mpf(1) / 3) / e + mp.mpf(4) / 5
        assert abs(first - 2 * ell / 3) < mp.mpf("1e-22")
        assert abs(second - (mp.mpf(2) / 3 - ell / 3)) < mp.mpf("1e-22")
        wrong_first = (-2 * J0 / 3 + mp.mpf(2) / 3) / e
        wrong_second = (J1 / 3 - mp.mpf(1) / 3) / e + mp.mpf(4) / 5
        assert abs(wrong_first - first) > mp.mpf("0.44")
        assert abs(wrong_second - second) > mp.mpf("0.22")


@pytest.mark.parametrize("time", ["-10", "-2", "-0.5", "0", "0.5", "2", "10"])
def test_full_profile_pressure_finite_term_bound(time):
    with mp.workdps(75):
        x = mp.mpf(time)
        m = mp.mpf(10) ** 200
        delta = 3 * mp.mpf(10) ** 197
        tau = mp.mpf(10) ** -100
        f = lambda z: z / (1 + z**8) ** (mp.mpf(1) / 8)
        one = mp.diff(f, x)
        two = mp.diff(f, x, 2)
        assert abs(one) <= 1 and abs(two) <= 9
        L = 2 * delta / (m - delta)
        cap = (
            6 * (m + delta) * delta * L + delta * delta * (L / 3 + mp.mpf(2) / 3)
        ) / (tau * tau)
        for sign in (-1, 1):
            M = m + sign * delta * f(x)
            v = sign * delta / tau * one
            a = sign * delta / (tau * tau) * two
            ell = mp.log(M * M / (m * m))
            exact_term = 2 * M * a * ell / 3 - v * v * ell / 3 + 2 * v * v / 3
            assert abs(exact_term) <= cap


def test_curvature_improvement_has_zero_energy_but_nonzero_pressure():
    with mp.workdps(75):
        m, v, a, mu = mp.mpf(2), mp.mpf("0.1"), mp.mpf("0.03"), mp.mpf(3)
        M = lambda t: m + v * t + a * t * t / 2
        F = lambda t: M(t) ** 2 * (mp.log(M(t) ** 2 / (mu * mu)) - 1)
        Q = 16 * mp.pi**2
        ell = mp.log(m * m / (mu * mu))
        rho = -v * v * (ell + mp.mpf(2) / 3) / Q
        p = (2 * m * a * ell / 3 - v * v * ell / 3 + 2 * v * v / 3) / Q
        improvement = mp.diff(F, 0, 2) / (3 * Q)
        assert abs(p - rho - improvement) < mp.mpf("1e-65")
        assert abs(improvement) > mp.mpf("1e-7")
        assert p != rho
    t = s.Symbol("t", real=True)
    F = s.Function("F")(t)
    assert s.diff(F, t, 2) - s.diff(F, t, 2) == 0


def test_rotational_second_moments_and_rest_frame_tensor():
    directions = [
        s.Matrix(v)
        for v in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    ]
    assert sum(directions, s.zeros(3, 1)) == s.zeros(3, 1)
    assert sum((v * v.T for v in directions), s.zeros(3)) / 6 == s.eye(3) / 3
    assert "fixed Minkowski orthonormal" in pressure.data()["homogeneous_stress_tensor"]
    assert "not every Lorentz-boosted" in pressure.data()["homogeneous_stress_tensor"]


def test_zero_variation_and_exact_actual_stress_allowances():
    m = s.Integer(10) ** 200
    t = s.Rational(1, 10**100)
    assert all(value == 0 for value in pressure.enclosures(m, 0, t).values())
    data = pressure.data()
    values = data["all_actual_uniform_rational_allowances"]
    assert 0 < data["complete_subtracted_pressure_remainder_upper"] < 10**410
    assert 0 < values["restored_MS_derivative_pressure_upper"] < 10**595
    assert (
        values["inherited_absolute_energy_upper"]
        < values["complete_absolute_pressure_upper"]
        < 10**789
    )
    assert data["all_stress_components_over_named_kappa"] < s.Rational(1, 10**11)
    assert "not the full physical interacting parent stress" in data["scope"]
    assert (
        "insensitive to homogeneous flat curvature improvements"
        in data["curvature_choice"]
    )


def test_same_frontier_and_exact_serializer_boundaries():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert "NOT_FULL_INTERACTING_CURVED_PARENT_STRESS" in audit.matching()[-1]["status"]
    for mod in audit.MODULES:
        exact.serialize({k: v for k, v in mod.data().items() if k != "checks"})
    assert all(
        type(v) is str for v in pressure.data()["decimal_diagnostics_only"].values()
    )
    assert all(
        v.is_Rational
        for v in pressure.data()["all_actual_uniform_rational_allowances"].values()
    )
    with pytest.raises(ValueError, match="Inexact or nonfinite"):
        exact.serialize(s.Float(1))
