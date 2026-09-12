"""Independent covariant polarizations, tensor projection and cut checks."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_flat_tensor_cut import audit, curvature, dispersion
from p8_vacuum_affine_flat_tensor_cut import polarizations as pol
from p8_vacuum_affine_flat_tensor_cut import projectors as pro

ETA = s.diag(1, -1, -1, -1)


def covariant_modes(sign, m, k, E):
    momentum = s.Matrix([E, 0, 0, -sign * k])
    potentials = (
        s.Matrix([0, 1, 0, 0]),
        s.Matrix([0, 0, 1, 0]),
        s.Matrix([-k / m, 0, 0, sign * E / m]),
    )
    return momentum, potentials


def covariant_pair(k, A, l, B, m, metric=ETA):
    F = -s.I * (k * A.T - A * k.T)
    G = -s.I * (l * B.T - B * l.T)
    return (
        -(F * metric * G.T + G * metric * F.T) / 2
        + metric * s.trace(F * metric * G.T * metric) / 4
        + m * m * (A * B.T + B * A.T) / 2
        - metric * m * m * (A.T * metric * B)[0] / 2
    )


@pytest.mark.parametrize("r,c", tuple(product(range(3), repeat=2)))
def test_literal_covariant_four_tensor_reconstructs_every_physical_pair(r, c):
    k, m = pol.MOMENTUM, pol.MASS
    E = s.sqrt(k * k + m * m)
    K, A = covariant_modes(1, m, k, E)
    L, B = covariant_modes(-1, m, k, E)
    direct = covariant_pair(K, A[r], L, B[c], m)
    actual = pol.pairs()[3 * r + c].subs(pol.ENERGY, E)
    assert all(s.simplify(x) == 0 for x in direct - actual)
    assert (K.T * ETA * A[r])[0] == 0 and s.simplify((A[r].T * ETA * A[r])[0] + 1) == 0


@pytest.mark.parametrize("sign", (-1, 1))
def test_full_covariant_polarization_completeness(sign):
    k, m = pol.MOMENTUM, pol.MASS
    E = s.sqrt(k * k + m * m)
    K, A = covariant_modes(sign, m, k, E)
    result = sum((a * a.T for a in A), s.zeros(4))
    assert all(s.simplify(x) == 0 for x in result + ETA - K * K.T / (m * m))


@pytest.mark.parametrize("r,c", tuple(product(range(3), repeat=2)))
def test_noncollinear_Lorentz_boost_retains_the_full_stress_and_Ward_identity(r, c):
    m = s.Integer(3)
    k = s.Integer(4)
    E = s.Integer(5)
    K, A = covariant_modes(1, m, k, E)
    L, B = covariant_modes(-1, m, k, E)
    n = s.Matrix([s.Rational(1, 3), s.Rational(2, 3), s.Rational(2, 3)])
    gamma = s.Rational(5, 4)
    speed = s.Rational(3, 5)
    boost = s.eye(4)
    boost[0, 0] = gamma
    boost[0, 1:4] = (gamma * speed * n).T
    boost[1:4, 0] = gamma * speed * n
    boost[1:4, 1:4] = s.eye(3) + (gamma - 1) * n * n.T
    transform = ETA * boost * ETA
    assert boost.T * ETA * boost == ETA
    actual = covariant_pair(
        transform * K, transform * A[r], transform * L, transform * B[c], m
    )
    expected = transform * covariant_pair(K, A[r], L, B[c], m) * transform.T
    assert actual == expected
    total = transform * (K + L)
    assert total.T * ETA * actual == s.zeros(1, 4)


def mp_pair(k, A, l, B, m):
    metric = mp.diag([1, -1, -1, -1])
    F = -1j * (k * A.T - A * k.T)
    G = -1j * (l * B.T - B * l.T)
    trace = lambda M: sum(M[i, i] for i in range(M.rows))
    return (
        -(F * metric * G.T + G * metric * F.T) / 2
        + metric * trace(F * metric * G.T * metric) / 4
        + m * m * (A * B.T + B * A.T) / 2
        - metric * m * m * (A.T * metric * B)[0] / 2
    )


@pytest.mark.parametrize("ratio", ("0.01", "1", "100"))
def test_independent_full_angular_tensor_projection(ratio):
    with mp.workdps(65):
        m = mp.mpf(1)
        k = mp.mpf(ratio)
        E = mp.sqrt(m * m + k * k)
        points, weights = mp.gauss_quadrature(6, "legendre")
        tensor = mp.zeros(9)
        for z, w in zip(points, weights):
            radial = mp.sqrt(1 - z * z)
            for j in range(16):
                phi = 2 * mp.pi * j / 16
                n = mp.matrix([radial * mp.cos(phi), radial * mp.sin(phi), z])
                e1 = mp.matrix([z * mp.cos(phi), z * mp.sin(phi), -radial])
                e2 = mp.matrix([-mp.sin(phi), mp.cos(phi), 0])
                K = mp.matrix([E, *[-k * n[a] for a in range(3)]])
                L = mp.matrix([E, *[k * n[a] for a in range(3)]])
                A = [
                    mp.matrix([0, *e1]),
                    mp.matrix([0, *e2]),
                    mp.matrix([-k / m, *[E * n[a] / m for a in range(3)]]),
                ]
                B = [
                    mp.matrix([0, *e1]),
                    mp.matrix([0, *e2]),
                    mp.matrix([-k / m, *[-E * n[a] / m for a in range(3)]]),
                ]
                for r, c in product(range(3), repeat=2):
                    T = mp_pair(K, A[r], L, B[c], m)
                    v = mp.matrix([T[a + 1, b + 1] for a in range(3) for b in range(3)])
                    tensor += w * v * v.transpose_conj() / 32
        sigma = 4 * E * E
        a = (13 * sigma * sigma + 56 * m * m * sigma + 48 * m**4) / 120
        b = (sigma * sigma - 4 * m * m * sigma + 12 * m**4) / 12
        I = mp.eye(3)
        expected = mp.matrix(9)
        for u, v in product(range(9), repeat=2):
            i, j = u // 3, u % 3
            r, c = v // 3, v % 3
            P0 = I[i, j] * I[r, c] / 3
            P2 = (I[i, r] * I[j, c] + I[i, c] * I[j, r]) / 2 - P0
            expected[u, v] = a * P2 + b * P0
        assert mp.norm(tensor - expected) < mp.mpf("1e-55") * mp.norm(expected)


def test_independent_Lorentz_phase_space_Jacobian():
    r, m = s.symbols("r mass", positive=True)
    E = s.sqrt(r * r + m * m)
    angular = 4 * s.pi
    measure = angular * r * r / ((2 * s.pi) ** 2 * (2 * E) ** 2)
    jac = s.diff(2 * E, r)
    phase = s.simplify(measure / jac)
    assert s.simplify(phase - r / (8 * s.pi * E)) == 0


@pytest.mark.parametrize("E", (s.Integer(3), s.Integer(1000)))
def test_literal_oscillator_current_sign_and_spectral_quarter(E):
    t = s.Symbol("t", real=True)
    a = s.zeros(3)
    a[0, 1] = 1
    a[1, 2] = s.sqrt(2)
    q = (a * s.exp(-s.I * E * t) + a.T * s.exp(s.I * E * t)) / s.sqrt(2 * E)
    O = q * q
    O0 = O.subs(t, 0)
    current = s.expand_complex(s.I * (O * O0 - O0 * O)[0, 0] / 4).expand()
    assert s.trigsimp(current - s.sin(2 * E * t) / (4 * E * E)) == 0
    decay = s.Symbol("decay", positive=True)
    laplace = s.integrate(
        s.exp(-decay * t) * s.sin(2 * E * t) / (4 * E * E), (t, 0, s.oo)
    )
    assert s.simplify(laplace - 1 / (2 * E * (4 * E * E + decay * decay))) == 0
    # Connected <q^2(t)q^2(0)> gives2/(2E)^2. Converting its
    # positive-frequency delta to delta(s-4E^2), then dividing by4,
    # yields the same spectral coefficient1/(2E).
    assert (2 / (2 * E) ** 2) * (4 * E) / 4 == 1 / (2 * E)


@pytest.mark.parametrize("spin", (0, 2))
@pytest.mark.parametrize("ratio", (s.Rational(1, 2), 1, 3, 4, 5, 100))
def test_full_density_support_and_strict_positive_massive_cut(spin, ratio):
    value = pro.density(spin).subs(
        {pro.MASS: 1000, pro.S: s.Integer(1000) ** 2 * ratio}
    )
    if ratio <= 4:
        assert value == 0
    else:
        assert value > 0 and value.is_real is True


def dimensionless_density_polynomial(spin, v):
    x = (1 - v * v) / 4
    return (
        (13 + 56 * x + 48 * x * x, mp.mpf(3840))
        if spin == 2
        else (1 - 4 * x + 12 * x * x, mp.mpf(384))
    )


@pytest.mark.parametrize("spin", (0, 2))
@pytest.mark.parametrize("order", (3, 4))
def test_independent_complete_spectral_moments(spin, order):
    with mp.workdps(75):

        def integrand(v):
            polynomial, den = dimensionless_density_polynomial(spin, v)
            return (
                v
                * v
                * polynomial
                / (2 * den * mp.pi**2)
                * ((1 - v * v) / 4) ** (order - 3)
            )

        value = mp.quad(integrand, [0, 1])
        expected = mp.mpf(
            str(s.N(dispersion.moment(spin, order).subs(pro.MASS, 1), 70))
        )
        assert abs(value - expected) < mp.mpf("1e-65") * expected


@pytest.mark.parametrize("spin", (0, 2))
@pytest.mark.parametrize("point", ("1", "-1", "2", "-2", "complex"))
def test_independent_subtracted_dispersion_and_next_order_error(spin, point):
    with mp.workdps(75):
        z = mp.mpc(1, 1) if point == "complex" else mp.mpf(point)

        def integrand(v):
            polynomial, den = dimensionless_density_polynomial(spin, v)
            return 2 * v * v * polynomial / (den * mp.pi**2 * (4 - z * (1 - v * v)))

        actual = z**3 * mp.quad(integrand, [0, 1])
        I = mp.mpf(str(s.N(dispersion.moment(spin, 3).subs(pro.MASS, 1), 70)))
        J = mp.mpf(str(s.N(dispersion.moment(spin, 4).subs(pro.MASS, 1), 70)))
        assert 0 < abs(actual) < 2 * abs(z) ** 3 * I
        assert 0 < abs(actual - z**3 * I) < 2 * abs(z) ** 4 * J
        if point not in ("complex",):
            assert mp.sign(mp.re(actual)) == mp.sign(z)


@pytest.mark.parametrize("spin", (0, 2))
@pytest.mark.parametrize("lower", (4, 8, 100))
def test_independent_discarded_spectral_tail_bound(spin, lower):
    with mp.workdps(70):
        L = mp.mpf(lower)
        z = mp.mpf(2)
        start = mp.sqrt(1 - 4 / L)

        def integrand(v):
            polynomial, den = dimensionless_density_polynomial(spin, v)
            return 2 * v * v * polynomial / (den * mp.pi**2 * (4 - z * (1 - v * v)))

        tail_value = z**3 * mp.quad(integrand, [start, 1])
        den = 128 if spin == 2 else 384
        assert 0 < tail_value < 2 * abs(z) ** 3 / (den * mp.pi**2 * L)


@pytest.mark.parametrize("fixture", (0, 1, 2))
def test_independent_four_dimensional_linear_curvature_tensor_Hessians(fixture):
    w = s.Symbol("frequency", positive=True)
    p = s.Matrix([w, 0, 0, 0])
    metric = ETA
    D = s.Matrix([[1 + fixture, 2, -1], [2, -3, 1], [-1, 1, 4 - fixture]])
    G = s.Matrix([[2, -1, 3], [-1, 1 + fixture, 2], [3, 2, -2]]) / 7

    def linear_curvature(H):
        h = s.zeros(4)
        h[1:4, 1:4] = H
        R = {}
        for a, b, c, d in product(range(4), repeat=4):
            R[a, b, c, d] = (
                -(
                    p[c] * p[b] * h[a, d]
                    + p[d] * p[a] * h[b, c]
                    - p[d] * p[b] * h[a, c]
                    - p[c] * p[a] * h[b, d]
                )
                / 2
            )
        Ric = s.Matrix(
            4, 4, lambda a, b: sum(metric[r, r] * R[r, a, r, b] for r in range(4))
        )
        scalar = s.trace(metric * Ric)
        return R, Ric, scalar

    A, RA, sa = linear_curvature(D)
    B, RB, sb = linear_curvature(G)
    riem = sum(
        metric[a, a]
        * metric[b, b]
        * metric[c, c]
        * metric[d, d]
        * A[a, b, c, d]
        * B[a, b, c, d]
        for a, b, c, d in product(range(4), repeat=4)
    )
    ric = sum(
        metric[a, a] * metric[b, b] * RA[a, b] * RB[a, b]
        for a, b in product(range(4), repeat=2)
    )
    Weyl_mixed = 2 * (riem - 2 * ric + sa * sb / 3)
    assert (
        s.expand(Weyl_mixed - w**4 * (s.trace(D * G) - s.trace(D) * s.trace(G) / 3))
        == 0
    )
    assert s.expand(2 * sa * sb - 2 * w**4 * s.trace(D) * s.trace(G)) == 0


@pytest.mark.parametrize("bad", (True, False, 1, 3, 1.0, s.Float(2), "2", None))
def test_nonphysical_spin_selector_rejected(bad):
    with pytest.raises(ValueError):
        pro.density(bad)
    with pytest.raises(ValueError):
        dispersion.moment(bad, 3)


@pytest.mark.parametrize("bad", (True, 3.0, s.Float(3), 2, 0, -1, "3", None))
def test_nonconvergent_or_nonexact_moment_order_rejected(bad):
    with pytest.raises(ValueError):
        dispersion.moment(2, bad)


@pytest.mark.parametrize(
    "name", ("polarizations", "projectors", "dispersion", "curvature")
)
def test_exact_scientific_packet(name):
    p = {
        "polarizations": pol,
        "projectors": pro,
        "dispersion": dispersion,
        "curvature": curvature,
    }[name].data()
    for v in p["checks"].values():
        entries = list(v) if isinstance(v, s.MatrixBase) else [v]
        assert all(s.cancel(x) == 0 for x in entries)
    assert all(p["gates"].values())


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_exact_residuals(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.cancel(entry) == 0 for entry in entries)


@pytest.mark.parametrize("case", audit.bad_cases(), ids=lambda case: case[0])
def test_all_scope_mutations_rejected(case):
    _, call, args = case
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_counts_and_guards():
    assert len(audit.residuals()) == 48
    assert audit.scalar_entry_count() == 450
    assert len(audit.gates()) == 31
    assert all(value is True for value in audit.gates().values())
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 121


def test_flat_benchmark_does_not_change_state_or_close_matching():
    assert "actual CD state is not changed" in audit.observable()["domain"]
    assert (
        "not a dimensionally continued finite local action"
        in curvature.data()["comparison_boundary"]
    )
    assert (
        "not asserted to be a canonically reduced propagating scalar"
        in curvature.data()["canonical"]
    )
    assert "finite local counterterms" in dispersion.data()["boundary"]
    assert audit.ITEM["status"].endswith(
        "NOT_CURVED_MATCHING_FULL_PARENT_AMPLITUDE_OR_V_G_B"
    )
    assert len(audit.frontier()) == 9
