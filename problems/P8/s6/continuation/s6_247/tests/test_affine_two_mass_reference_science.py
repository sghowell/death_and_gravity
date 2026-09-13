"""Independent full scalar stress, curvature, two-threshold and causal inverse tests."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_two_mass_reference import audit, inverse, spectral
from p8_vacuum_affine_two_mass_reference import geometry as g


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_complete_exact_identity(name):
    value = audit.residuals()[name]
    assert all(x == 0 for x in value) if isinstance(value, s.MatrixBase) else value == 0


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_every_written_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_input_is_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_actual_parameters_frontier_and_scope_are_preserved():
    assert audit.require_parameters(g.PROCA_MASS2, g.HEAVY_MASS2, g.KAPPA) == (
        g.PROCA_MASS2,
        g.HEAVY_MASS2,
        g.KAPPA,
    )
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 103
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching()) is True
    assert all(v is True for k, v in audit.controls().items() if k != "rejected_inputs")
    assert len(audit.residuals()) == 103 and audit.scalar_entry_count() == 205
    assert len(audit.gates()) == 37 and audit.rejected_inputs() == 181


@pytest.mark.parametrize("case", range(6))
def test_independent_scalar_Lagrangian_stress_and_noncollinear_covariant_Ward(case):
    m, k, E = ((3, 4, 5), (5, 12, 13), (8, 15, 17))[case // 2]
    eta = s.diag(1, -1, -1, -1)
    direction = (
        s.Matrix([s.Rational(1, 3), s.Rational(2, 3), s.Rational(2, 3)])
        if case % 2 == 0
        else s.Matrix([0, s.Rational(3, 5), s.Rational(4, 5)])
    )
    gamma, speed = (
        (s.Rational(5, 4), s.Rational(3, 5))
        if case % 2 == 0
        else (s.Rational(13, 12), s.Rational(5, 13))
    )
    boost = s.eye(4)
    boost[0, 0] = gamma
    boost[0, 1:4] = (gamma * speed * direction).T
    boost[1:4, 0] = gamma * speed * direction
    boost[1:4, 1:4] = s.eye(3) + (gamma - 1) * direction * direction.T
    assert boost.T * eta * boost == eta
    transform = eta * boost * eta
    K, L = s.Matrix([E, 0, 0, -k]), s.Matrix([E, 0, 0, k])

    def lagrangian_pair(a, b):
        da, db = -s.I * a, -s.I * b
        return (da * db.T + db * da.T) / 2 - eta * ((da.T * eta * db)[0] - m * m) / 2

    original = lagrangian_pair(K, L)
    actual = lagrangian_pair(transform * K, transform * L)
    assert actual == transform * original * transform.T
    assert actual == g.stress_pair(transform * K, transform * L, m * m)
    assert ((transform * (K + L)).T * eta * actual) == s.zeros(1, 4)
    assert original[0, :] == s.zeros(1, 4)
    assert s.trace(original[1:4, 1:4]) == -2 * k * k - 3 * m * m


@pytest.mark.parametrize("ratio", ("0.01", "1", "100"))
def test_independent_full_nine_by_nine_scalar_angular_tensor(ratio):
    with mp.workdps(65):
        k = mp.mpf(ratio)
        m = mp.mpf(3)
        E = mp.sqrt(k * k + m * m)
        roots, weights_ = mp.gauss_quadrature(6, "legendre")
        result = mp.zeros(9)
        for z0, w0 in zip(roots, weights_):
            r = mp.sqrt(1 - z0 * z0)
            for j in range(12):
                phi = 2 * mp.pi * j / 12
                nhat = mp.matrix([r * mp.cos(phi), r * mp.sin(phi), z0])
                pair = k * k * nhat * nhat.T - E * E * mp.eye(3)
                v = mp.matrix([pair[a, b] for a, b in product(range(3), repeat=2)])
                result += w0 * v * v.T / 24
        sig = 4 * E * E
        a = (sig - 4 * m * m) ** 2 / 120
        b = (sig + 2 * m * m) ** 2 / 12
        expected = mp.zeros(9)
        I = mp.eye(3)
        for v, w in product(range(9), repeat=2):
            i, j = v // 3, v % 3
            r, c = w // 3, w % 3
            p0 = I[i, j] * I[r, c] / 3
            p2 = (I[i, r] * I[j, c] + I[i, c] * I[j, r]) / 2 - p0
            expected[v, w] = a * p2 + b * p0
        assert mp.norm(result - expected) < mp.mpf("1e-55") * mp.norm(expected)


@pytest.mark.parametrize("dimension", range(2, 7))
@pytest.mark.parametrize("case", (0, 1))
def test_independent_complete_dimension_sphere_not_tracefree_only(dimension, case):
    v = s.symbols("nhat0:" + str(dimension))
    D = s.diag(*range(1, dimension + 1))
    G = s.diag(*range(dimension, 0, -1))
    D[0, 1] = D[1, 0] = case + 1
    G[0, dimension - 1] = G[dimension - 1, 0] = s.Rational(1, case + 2)
    assert D * G != G * D
    vec = s.Matrix(v)
    f = (s.trace(D) - (vec.T * D * vec)[0]) * (s.trace(G) - (vec.T * G * vec)[0]) / 4
    actual = s.S.Zero
    for powers, coef in s.Poly(s.expand(f), *v).terms():
        if any(p % 2 for p in powers):
            continue
        half = sum(powers) // 2
        numerator = s.prod(s.factorial2(p - 1) for p in powers)
        denominator = s.prod(dimension + 2 * j for j in range(half))
        actual += coef * numerator / denominator
    C = s.Rational(1, 2 * dimension * (dimension + 2))
    B = s.Rational(dimension * dimension - 3, 4 * dimension * (dimension + 2))
    assert s.factor(actual - C * s.trace(D * G) - B * s.trace(D) * s.trace(G)) == 0


def independent_homogeneous_curvature(Q):
    eta = (1, -1, -1, -1)
    h = s.zeros(4)
    h[1:4, 1:4] = -Q
    derivative = lambda a, b, c, d: h[a, b] if c == d == 0 else s.S.Zero
    R = {}
    for a, b, c, d in product(range(4), repeat=4):
        R[a, b, c, d] = (
            derivative(a, d, b, c)
            + derivative(b, c, a, d)
            - derivative(a, c, b, d)
            - derivative(b, d, a, c)
        ) / 2
    Ric = s.Matrix(4, 4, lambda b, d: sum(eta[a] * R[a, b, a, d] for a in range(4)))
    scalar = sum(eta[a] * Ric[a, a] for a in range(4))
    return R, Ric, scalar


@pytest.mark.parametrize("case", range(3))
def test_independent_literal_curvature_Hessians_keep_action_factor_two(case):
    D = (s.diag(1, -1, 0), 2 * s.eye(3), s.Matrix([[2, 1, 0], [1, -1, 2], [0, 2, 3]]))[
        case
    ]
    G = (
        s.Matrix([[2, 1, 0], [1, 1, 0], [0, 0, -3]]),
        2 * s.eye(3),
        s.Matrix([[1, 0, 2], [0, 2, 1], [2, 1, -2]]),
    )[case]
    RD, PD, SD = independent_homogeneous_curvature(D)
    RG, PG, SG = independent_homogeneous_curvature(G)
    eta = (1, -1, -1, -1)
    riem = sum(s.prod(eta[i] for i in inds) * RD[inds] * RG[inds] for inds in RD)
    ric = sum(
        eta[a] * eta[b] * PD[a, b] * PG[a, b] for a, b in product(range(4), repeat=2)
    )
    weyl = riem - 2 * ric + SD * SG / 3
    euler = riem - 4 * ric + SD * SG
    assert riem == s.trace(D * G)
    assert euler == 0
    assert weyl == s.trace(D * G) / 2 - s.trace(D) * s.trace(G) / 6
    ell = s.Symbol("ell")
    Hessian = -2 * ell * (SD * SG / 72 + (riem - ric) / 180) * 2
    expected = -ell * weyl / 30 - ell * SD * SG / 18
    assert s.expand(Hessian - expected) == 0
    if case == 1:
        assert s.simplify(Hessian + 2 * ell) == 0


W = {
    ("Proca", 0): [0, 3, -2, 3],
    ("Proca", 2): [0, 1, -s.Rational(2, 3), s.Rational(1, 10)],
    ("heavy", 0): [0, 9, -6, 1],
    ("heavy", 2): [0, 0, 0, s.Rational(1, 30)],
}


def coefficients(key):
    return tuple(mp.mpf(str(s.numer(c))) / int(s.denom(c)) for c in W[key])


def weight(key, y):
    return sum(c * y ** (2 * j) for j, c in enumerate(coefficients(key)))


def moment(key, j):
    return sum(
        c * mp.beta(mp.mpf(k) + mp.mpf(".5"), j + 1) / 2
        for k, c in enumerate(coefficients(key))
    )


def parts(key, z):
    constant = mp.mpc(0)
    coefficient = mp.mpc(0)
    pol = mp.mpc(0)
    mon = mp.mpc(1)
    for k, c in enumerate(coefficients(key)):
        if k:
            pol = z * pol - mp.mpf(1) / (2 * k - 1)
            mon *= z
        constant += c * pol
        coefficient += c * mon
    return constant, coefficient


def small_integral(key, r):
    ans = mp.mpc(0)
    term = r
    M = sum(abs(c) / (2 * k + 1) for k, c in enumerate(coefficients(key)))
    for j in range(1000):
        ans += term * moment(key, j)
        tail = abs(r) ** (j + 2) * M / (1 - abs(r))
        if tail < mp.eps * max(abs(ans), mp.mpf("1e-1000")) / 100:
            return ans
        term *= -r
    raise AssertionError("series failed")


def closed(key, p, m2):
    if p == 0:
        return mp.mpc(0)
    r = p / (4 * m2)
    if abs(r) < mp.mpf(".2"):
        return small_integral(key, r)
    z = 1 + 1 / r
    pol, coef = parts(key, z)
    if not mp.im(p) and mp.re(p) < 0 and mp.re(p) > -4 * m2:
        a = mp.sqrt(-z)
        base = -mp.atan(1 / a) / a
    else:
        a = mp.sqrt(z)
        base = mp.log((a + 1) ** 2 * r) / (2 * a)
    return pol + coef * base


def direct(key, p, m2):
    if p == 0:
        return mp.mpc(0)
    r = p / (4 * m2)
    turn = max(mp.mpf(1), mp.log(1 + abs(r)))
    scale = r if abs(r) < 1 else mp.mpf(1)

    def f(v):
        u = mp.exp(-v)
        y = 1 - u
        c = u * (2 - u)
        return (r / scale) * weight(key, y) * u / (1 + r * c)

    knots = sorted({mp.mpf(0), mp.mpf(1), mp.mpf(4), turn, turn + 4, turn + 16})
    return scale * mp.quad(f, knots + [mp.inf])


def bank(key, u):
    if u < 0:
        return mp.re(closed(key, -mp.exp(u), mp.mpf(".25"))), mp.mpf(0)
    if u == 0:
        return mp.re(parts(key, 0)[0]), mp.mpf(0)
    z = -mp.expm1(-u)
    b = mp.sqrt(z)
    pol, coef = parts(key, z)
    real = pol + coef * (u / 2 + mp.log1p(b)) / b
    return mp.re(real), mp.re(coef / (2 * b))


def total(spin, p, mu2, n):
    ell = mp.log(n)
    c = 2 * (ell + 2) if spin == 0 else (ell + 2) / 60
    return c + closed(("Proca", spin), p, mu2) + closed(("heavy", spin), p, n)


def density_u(spin, u, mu2, n):
    ell = mp.log(n)
    shift = mp.log(n / mu2)
    c = 2 * (ell + 2) if spin == 0 else (ell + 2) / 60
    dp, up = bank(("Proca", spin), u)
    dh, uh = bank(("heavy", spin), u - shift)
    D, U = c + dp + dh, up + uh
    return U / (D * D + mp.pi**2 * U * U)


def dispersion(spin, p, mu2, n):
    shift = mp.log(n / mu2)
    knots = [
        mp.mpf(0),
        mp.mpf(1),
        mp.mpf(4),
        shift / 2,
        shift - 4,
        shift,
        shift + 1,
        shift + 4,
        shift + 32,
    ]
    knots = sorted({x for x in knots if x >= 0})
    integrand = lambda u: density_u(spin, u, mu2, n) / (1 + p / (4 * mu2) * mp.exp(-u))
    return mp.quad(integrand, knots + [mp.inf])


@pytest.mark.parametrize("species,spin", tuple(W))
@pytest.mark.parametrize("case", range(4))
def test_complete_radial_integral_against_independent_closed_function(
    species, spin, case
):
    with mp.workdps(70):
        actual_n = mp.mpf(10) ** 200 / 512 + 2
        m2, pvalue = (
            (mp.mpf(7), mp.mpc(3, 2)),
            (mp.mpf(7), mp.mpf(-5)),
            (actual_n, mp.mpc("1e6", "2e6")),
            (mp.mpf("1e6"), mp.mpc("1e205", "3e204")),
        )[case]
        a, b = closed((species, spin), pvalue, m2), direct((species, spin), pvalue, m2)
        assert abs(a - b) < mp.mpf("1e-55") * max(abs(b), mp.mpf("1e-1000"))


@pytest.mark.parametrize("mass_case", (0, 1))
@pytest.mark.parametrize("spin", (0, 2))
@pytest.mark.parametrize("case", range(3))
def test_entire_two_threshold_reciprocal_dispersion_including_actual_mass(
    mass_case, spin, case
):
    with mp.workdps(70):
        mu2, n = (
            (mp.mpf(1), mp.mpf("1e30"))
            if mass_case == 0
            else (mp.mpf("1e6"), mp.mpf(10) ** 200 / 512 + 2)
        )
        pvalue = (mp.mpf(0), mp.mpc(mu2, mu2 / 3), mp.mpc(n, n / 5))[case]
        target = 1 / total(spin, pvalue, mu2, n)
        actual = dispersion(spin, pvalue, mu2, n)
        assert abs(target - actual) < mp.mpf("1e-50") * abs(target)


@pytest.mark.parametrize("mass_case", (0, 1))
def test_second_trace_threshold_has_both_nonzero_one_sided_cusps(mass_case):
    with mp.workdps(80):
        mu2, n = (
            (mp.mpf(1), mp.mpf("1e30"))
            if mass_case == 0
            else (mp.mpf("1e6"), mp.mpf(10) ** 200 / 512 + 2)
        )
        h = mp.log(n / mu2)
        ell = mp.log(n)
        dp, up = bank(("Proca", 0), h)
        dh, uh = bank(("heavy", 0), 0)
        D, U = 2 * (ell + 2) + dp + dh, up + uh
        den = D * D + mp.pi**2 * U * U
        expected_minus = -9 * mp.pi * D * U / den**2
        expected_plus = mp.mpf(9) / 2 * (D * D - mp.pi**2 * U * U) / den**2
        eps = mp.mpf("1e-30")
        at = density_u(0, h, mu2, n)
        minus = (density_u(0, h - eps, mu2, n) - at) / mp.sqrt(eps)
        plus = (density_u(0, h + eps, mu2, n) - at) / mp.sqrt(eps)
        assert abs(minus - expected_minus) < mp.mpf("1e-10") * abs(expected_minus)
        assert abs(plus - expected_plus) < mp.mpf("1e-10") * abs(expected_plus)
        assert expected_minus < 0 < expected_plus


@pytest.mark.parametrize("case", range(4))
def test_complete_spectral_reconstruction_in_both_quotient_matrix_inverse_products(
    case,
):
    with mp.workdps(70):
        mu2 = mp.mpf("1e6")
        n = mp.mpf(10) ** 200 / 512 + 2
        L, qvalue = (
            (mp.mpc(1000, 100), mp.mpf(0)),
            (mp.mpc(1000, 200), mp.mpf("1e6")),
            (mp.mpc(3000, 100), mp.mpf("1e7")),
            (mp.sqrt(n) * mp.mpc(1, mp.mpf(".2")), n),
        )[case]
        pvalue = L * L + qvalue
        B = mp.matrix([[L * L + 2 * qvalue / 3, -L * L / 3], [-qvalue, -L * L]])
        forward = (
            B.T
            * mp.diag(
                [-total(0, pvalue, mu2, n), -mp.mpf(8) / 3 * total(2, pvalue, mu2, n)]
            )
            * B
        )
        Binv = mp.matrix(
            [
                [1 / pvalue, -1 / (3 * pvalue)],
                [1 / pvalue - 1 / L**2, -1 / (3 * pvalue) - mp.mpf(2) / (3 * L**2)],
            ]
        )
        recovered = (
            Binv
            * mp.diag(
                [
                    -dispersion(0, pvalue, mu2, n),
                    -mp.mpf(3) / 8 * dispersion(2, pvalue, mu2, n),
                ]
            )
            * Binv.T
        )
        assert mp.norm(forward * recovered - mp.eye(2)) < mp.mpf("1e-45")
        assert mp.norm(recovered * forward - mp.eye(2)) < mp.mpf("1e-45")


@pytest.mark.parametrize("spin", (0, 2))
def test_actual_lowest_threshold_and_old_pole_are_not_transferred(spin):
    with mp.workdps(75):
        mu2 = mp.mpf("1e6")
        n = mp.mpf(10) ** 200 / 512 + 2
        ell = mp.log(n)
        cp = 4 if spin == 0 else mp.mpf(1) / 30
        pvalue = -mp.mpf(".57985") * mu2
        P = cp + closed(("Proca", spin), pvalue, mu2)
        H = (2 * ell if spin == 0 else ell / 60) + closed(("heavy", spin), pvalue, n)
        actual = 1 / (P + H)
        assert mp.re(P + H) > 5
        assert abs(actual - (1 / P + 1 / H)) > mp.mpf(".001")
        # The independently specified old Proca shear factor still has its zero.
        if spin == 2:
            a = cp + closed(("Proca", spin), -mp.mpf(".5798") * mu2, mu2)
            b = cp + closed(("Proca", spin), -mp.mpf(".5799") * mu2, mu2)
            assert mp.re(a) > 0 > mp.re(b)


@pytest.mark.parametrize("qvalue", (0, 1, 16, 10000))
def test_independent_coordinate_oscillator_and_complete_initial_boundary(qvalue):
    t = s.Symbol("t", nonnegative=True)
    q = s.Integer(qvalue)
    wave = t if q == 0 else s.sin(s.sqrt(q) * t) / s.sqrt(q)
    R = s.Matrix([[wave, -wave / 3], [wave - t, -wave / 3 - 2 * t / 3]])
    assert R.subs(t, 0) == s.zeros(2)
    assert R.diff(t).subs(t, 0) == s.Matrix([[1, -s.Rational(1, 3)], [0, -1]])
    actual = inverse.coordinate_kernel(qvalue).subs(inverse.t, t)
    assert actual == R
    # Applying the full differential coordinate map includes the delta from R'(0).
    second = R.diff(t, 2)
    regular = s.Matrix(
        [
            [second[0, j] + 2 * q * R[0, j] / 3 - second[1, j] / 3 for j in range(2)],
            [-q * R[0, j] - second[1, j] for j in range(2)],
        ]
    )
    assert regular.applyfunc(s.simplify) == s.zeros(2)
    principal = s.Matrix([[1, -s.Rational(1, 3)], [0, -1]])
    assert principal * R.diff(t).subs(t, 0) == s.eye(2)
    assert R.diff(t).subs(t, 0) * principal == s.eye(2)


def test_normalized_reference_bound_does_not_drop_physical_force_or_metric_norm():
    assert inverse.L1_BOUND.subs(g.ell, 394) == s.Rational(125, 4224)
    assert s.Rational(125, 4224) < s.Rational(3, 100)
    expression = inverse.data()["physical_force_inverse_bound"]
    assert (
        s.factor(
            expression / (64 * s.pi**2 * g.KAPPA) - inverse.L1_BOUND * inverse.T**4
        )
        == 0
    )
    assert expression.subs({g.ell: 462, inverse.T: 1}) > g.KAPPA
    E = s.Matrix([[12, -4], [-4, 4]])
    assert (E - 2 * s.eye(2)).det() == 4 and (14 * s.eye(2) - E).det() == 4


@pytest.mark.parametrize("species,spin", tuple(W))
def test_independent_closed_functions_match_complete_symbolic_polynomial_parts(
    species, spin
):
    with mp.workdps(70):
        poly, coef = spectral.polynomial_parts(species, spin)
        fn = s.lambdify(spectral.z, poly, "mpmath")
        gn = s.lambdify(spectral.z, coef, "mpmath")
        independent = sum(c * g.y ** (2 * j) for j, c in enumerate(W[species, spin]))
        actual = (g.W_P if species == "Proca" else g.W_H)[spin]
        assert s.expand(independent - actual) == 0
        for zvalue in (mp.mpc("1.3", ".2"), mp.mpc("1.1", "-.6")):
            pp, cc = parts((species, spin), zvalue)
            assert abs(fn(zvalue) - pp) < mp.mpf("1e-60")
            assert abs(gn(zvalue) - cc) < mp.mpf("1e-60")
