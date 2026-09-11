"""Independent coordinate, dense Gaussian, full-function and pole controls."""

import copy
from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_affine_vacuum_domain import family
from p8_vacuum_affine_vector_schur import (
    audit,
    gaussian,
    homogeneous,
    lorentz,
    vacuum,
    verify,
)


@pytest.mark.parametrize("name,value", tuple(audit.residuals().items()))
def test_every_exact_residual(name, value):
    if isinstance(value, s.MatrixBase):
        assert value == s.zeros(*value.shape), name
    else:
        assert s.simplify(value) == 0, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_strict_rejected_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@cache
def literal_coordinate_euler():
    x = s.symbols("x0:4", real=True)
    w = s.symbols("W0:4", real=True)
    source = s.symbols("S0:4", real=True)
    eta = (1, -1, -1, -1)
    jets = s.Matrix(4, 4, lambda i, j: s.Symbol(f"q{i}{j}"))
    F = jets - jets.T
    z = s.Rational(1, 2000)
    action = (
        -z * sum(eta[i] * eta[j] * F[i, j] ** 2 for i in range(4) for j in range(4)) / 4
    )
    action += sum(eta[i] * (w[i] - source[i]) ** 2 for i in range(4)) / 2
    actual_fields = [
        x[i] ** 3 + x[(i + 1) % 4] * x[(i + 2) % 4] ** 2 + x[0] * x[1] * x[2] * x[3]
        for i in range(4)
    ]
    sources = [(i + 1) * x[(i + 1) % 4] ** 2 + x[i] ** 2 for i in range(4)]
    substitutions = {w[i]: actual_fields[i] for i in range(4)}
    substitutions.update({source[i]: sources[i] for i in range(4)})
    substitutions.update(
        {jets[i, j]: s.diff(actual_fields[j], x[i]) for i in range(4) for j in range(4)}
    )
    euler = s.Matrix(
        [
            s.expand(
                s.diff(action, w[j]).subs(substitutions)
                - sum(
                    s.diff(s.diff(action, jets[i, j]).subs(substitutions), x[i])
                    for i in range(4)
                )
            )
            for j in range(4)
        ]
    )
    expected = s.Matrix(
        [
            eta[j]
            * (
                actual_fields[j]
                - sources[j]
                + z
                * sum(
                    eta[i]
                    * (
                        s.diff(actual_fields[j], x[i], 2)
                        - s.diff(actual_fields[i], x[i], x[j])
                    )
                    for i in range(4)
                )
            )
            for j in range(4)
        ]
    )
    divergence = s.expand(sum(s.diff(euler[i], x[i]) for i in range(4)))
    expected_div = s.expand(
        sum(eta[i] * s.diff(actual_fields[i] - sources[i], x[i]) for i in range(4))
    )
    return euler, expected, divergence, expected_div


@pytest.mark.parametrize("component", range(4))
def test_original_flat_action_coordinate_Euler_and_sourced_constraint(component):
    actual, expected, div, expected_div = literal_coordinate_euler()
    assert actual[component] == expected[component]
    assert div == expected_div
    assert expected_div != 0


@pytest.mark.parametrize(
    "p",
    (
        (0, 0, 0, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (1, 1, 0, 0),
        (2, 1, 1, 0),
        (3, s.Rational(1, 2), s.Rational(1, 3), s.Rational(1, 4)),
    ),
)
def test_dense_inverse_source_contact_and_norm(p):
    eta = s.diag(1, -1, -1, -1)
    v = s.Matrix(p)
    z = vacuum.MAX_ZETA
    dot = (v.T * eta * v)[0]
    matrix = s.eye(4) + z * (-dot * s.eye(4) + v * v.T * eta)
    inverse = matrix.inv(method="DM")
    expected = (s.eye(4) - z * v * v.T * eta) / (1 - z * dot)
    assert inverse == expected
    src = s.Matrix(
        [s.Rational(1, 2), -s.Rational(2, 3), s.Rational(3, 5), s.Rational(4, 7)]
    )
    solved = inverse * src
    stationary = (
        (solved.T * eta * matrix * solved)[0]
        - 2 * (solved.T * eta * src)[0]
        + (src.T * eta * src)[0]
    ) / 2
    schur = (src.T * eta * (s.eye(4) - inverse) * src)[0] / 2
    assert stationary == schur
    omitted_contact = stationary - (src.T * eta * src)[0] / 2
    assert omitted_contact != schur
    assert v.T * eta * solved == v.T * eta * src
    numeric = mp.matrix((s.eye(4) - inverse).tolist())
    assert mp.norm(numeric * mp.matrix(list(src))) <= float(
        lorentz.remainder_upper(4, z)
    ) * mp.norm(mp.matrix(list(src))) + mp.mpf("1e-14")


@pytest.mark.parametrize("p", ((0, 0, 0, 0), (1, 2, 0, 0), (2, 1, 1, 1)))
def test_literal_Euclidean_gaussian_integrals(p):
    with mp.workdps(55):
        vec = mp.matrix(p)
        O = mp.eye(4) + mp.mpf(1) / 7 * (mp.fdot(vec, vec) * mp.eye(4) - vec * vec.T)
        eigen, basis = mp.eigsy(O)
        src = mp.matrix([mp.mpf(1) / 2, -mp.mpf(2) / 3, mp.mpf(3) / 5, mp.mpf(4) / 7])
        rotated = basis.T * src
        integral = mp.mpf(1)
        for i in range(4):
            f = lambda w, index=i: mp.exp(
                -eigen[index] * w * w / 2 + rotated[index] * w - rotated[index] ** 2 / 2
            )
            integral *= mp.quad(f, [-mp.inf, mp.inf])
        predicted = (
            (2 * mp.pi) ** 2
            / mp.sqrt(mp.det(O))
            * mp.exp(-mp.fdot(src, (mp.eye(4) - O**-1) * src) / 2)
        )
        assert abs(integral / predicted - 1) < mp.mpf("1e-45")
        shifted = O**-1 * src
        assert mp.fdot(src, (mp.eye(4) - O**-1) * src) >= -mp.mpf("1e-50")
        assert mp.fdot(src, shifted) <= mp.fdot(src, src) + mp.mpf("1e-50")


@pytest.mark.parametrize("power", (1, 2, 4, 8, 16))
def test_transverse_near_pole_control(power):
    z = vacuum.MAX_ZETA
    eps = s.Rational(1, 10**power)
    p = s.Matrix([s.sqrt(1 / z + eps), 0, 0, 0])
    eta = s.diag(1, -1, -1, -1)
    O = (1 - z * (p.T * eta * p)[0]) * s.eye(4) + z * p * p.T * eta
    transverse = s.Matrix([0, 1, 0, 0])
    response = O.inv() * transverse
    assert response == -transverse / (z * eps)
    assert abs(response[1]) > 1 / (2 * z * eps)
    # This is a pointwise multiplier control, not a claim about S(Phi)'s image.


def test_nonsymmetric_boundary_and_closed_source_controls():
    err = gaussian.variation()["retarded_single_branch_error"]
    assert err != s.zeros(2, 1)
    d = homogeneous.data()
    assert "compatible vector data" in d["homogeneous_lift"]
    assert "compact smooth" in d["causal_theorem_hypotheses"]
    z = s.Rational(1, 2000)
    p = s.Matrix([1, 2, 0, 0])
    O = s.eye(4) + z * (
        -(p.T * s.diag(1, -1, -1, -1) * p)[0] * s.eye(4)
        + p * p.T * s.diag(1, -1, -1, -1)
    )
    assert O * p == p
    assert (p.T * s.diag(1, -1, -1, -1) * p)[0] != 0
    # A gradient source has a nonzero divergence but exactly zero Schur term.
    assert O.inv() * p == p


@cache
def actual_target_functions():
    d = family.data()
    delta = d["R"] - 1
    # Keep the exact symbolic difference: numerical R-1 would lose the tiny
    # full source at the actual kappa=10^800 in finite precision.
    values = (delta, d["RX"], s.diff(delta, family.u))
    return s.lambdify((family.u, family.X), values, "mpmath", cse=True)


def exact_source_at(phi, grad, hess, kappa, euclidean=False):
    scale = mp.sqrt(kappa)
    u = mp.mpf(phi) / scale
    g = mp.matrix([mp.mpf(v) / scale for v in grad])
    hh = mp.matrix([[mp.mpf(v) / scale for v in row] for row in hess])
    eta = mp.eye(4) if euclidean else mp.diag([1, -1, -1, -1])
    x = mp.fdot(g, eta * g)
    box = sum(eta[i, i] * hh[i, i] for i in range(4))
    Z = mp.fdot(eta * g, hh * (eta * g))
    # The real Euclidean counterpart has the continued Lorentz contractions.
    if euclidean:
        x, box = -x, -box
    delta, RX, Ru = actual_target_functions()(u, x)
    R = 1 + delta
    r = delta / (x * x) if x else -mp.mpf(1024) / (1 + u * u) ** 3
    H = 4 * u / (1 + u * u)
    bracket = (
        delta * (-3 * H + 3 * Ru / (4 * R))
        + x * r * box
        + (-r + 3 * x * r * RX / (2 * R)) * Z
    )
    return g * bracket, (x, delta, RX, Ru, r)


@pytest.mark.parametrize("kappa_text", ("65536", "1e12", "1e800"))
@pytest.mark.parametrize(
    "grad",
    (
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (1, 1, 0, 0),
        (1, 1, 1, 1),
        (s.Rational(1, 2), -1, s.Rational(1, 3), 0),
    ),
)
@pytest.mark.parametrize("euclidean", (False, True))
def test_complete_actual_source_pointwise_majorant(kappa_text, grad, euclidean):
    with mp.workdps(100):
        kappa = mp.mpf(kappa_text)
        gg = [
            mp.mpf(int(v.p)) / int(v.q) if isinstance(v, s.Rational) else mp.mpf(v)
            for v in grad
        ]
        hh = [[mp.mpf((i + j) % 3 - 1) for j in range(4)] for i in range(4)]
        result, (x, delta, RX, Ru, r) = exact_source_at(
            mp.mpf("0.7"), gg, hh, kappa, euclidean
        )
        n = 1024
        u = mp.mpf("0.7") / mp.sqrt(kappa)
        h = (1 + u * u) ** 3
        # This saturated numerical boundary has explicit relative roundoff;
        # the exact uniform inequality is proved from the Fourier jet bound.
        assert abs(x) <= 4 / kappa * (1 + mp.mpf("1e-90"))
        assert abs(delta) <= 2 * n * x * x / h
        assert abs(r) <= 2 * n / h
        assert abs(RX) <= 4 * n * abs(x) / h
        assert abs(Ru + 6 * u * delta / (1 + u * u)) < mp.mpf("1e-90") * max(
            abs(Ru), mp.mpf("1e-2000")
        )
        bound = (
            64 * n / kappa**2
            + 384 * n / kappa**3
            + 6144 * n * n / kappa**4
            + 9216 * n * n / kappa**5
        )
        for i in range(4):
            assert abs(result[i]) <= bound * abs(gg[i]) + mp.mpf("1e-100") * bound
        assert mp.sqrt(mp.fdot(result, result)) <= 256 * n / kappa**2


@pytest.mark.parametrize(
    "x_text", ("-0.00006103515625", "-1e-12", "0", "1e-12", "0.00006103515625")
)
@pytest.mark.parametrize("u_text", ("0", "0.1", "3"))
def test_full_analytic_coefficient_envelopes_both_X_signs(x_text, u_text):
    with mp.workdps(100):
        u, x = mp.mpf(u_text), mp.mpf(x_text)
        delta, RX, Ru = actual_target_functions()(u, x)
        n = 1024
        h = (1 + u * u) ** 3
        assert abs(delta) <= 2 * n * x * x / h
        assert abs(RX) <= 4 * n * abs(x) / h
        assert 1 + delta > mp.mpf("0.5")
        assert abs(Ru + 6 * u * delta / (1 + u * u)) < mp.mpf("1e-90")
        if x:
            T = x**n / (x**n + (1 - x) ** n)
            w = n * x * x * mp.exp(-n * x * x)
            B = T + (1 - T) * w
            assert 0 <= T <= 1 and 0 <= w <= 1
            assert B <= (n + 1) * x * x
            assert abs(delta - B * (x - 1) / h) < mp.mpf("1e-90")


@pytest.mark.parametrize("scale_text", ("1e-6", "1e-10"))
def test_full_function_recovers_leading_vacuum_source(scale_text):
    with mp.workdps(100):
        eps = mp.mpf(scale_text)
        grad = [mp.mpf("0.3"), mp.mpf("0.7"), -mp.mpf("0.4"), mp.mpf("0.2")]
        hh = [[mp.mpf(i + j - 3) / 5 for j in range(4)] for i in range(4)]
        eta = mp.diag([1, -1, -1, -1])
        gv = mp.matrix(grad)
        X = mp.fdot(gv, eta * gv)
        Box = sum(eta[i, i] * hh[i][i] for i in range(4))
        Z = mp.fdot(eta * gv, mp.matrix(hh) * (eta * gv))
        expected = 1024 * gv * (Z - X * Box)
        actual, _ = exact_source_at(mp.mpf("0.4"), grad, hh, eps**-2)
        assert mp.norm(expected) > 0
        assert mp.norm(actual / eps**4 - expected) / mp.norm(expected) < 10000 * eps**2


def test_actual_null_limit_from_frozen_second_X_derivative():
    d = family.data()
    coefficient = s.diff(d["R"], family.X, 2).subs(family.X, 0) / 2
    assert s.factor(coefficient + 1024 / (1 + family.u**2) ** 3) == 0


@pytest.mark.parametrize("n", (4, 6, 16, 1024))
def test_uniform_component_envelope_at_worst_allowed_kappa(n):
    n, k, _ = vacuum.require_small_domain(n, 64 * n, vacuum.MAX_ZETA)
    assert vacuum.component_bound(n, k) < 128 * n / k**2
    rho = 1 / (16 * n)
    assert 4 * (2 * rho) ** (n - 2) <= 1
    assert 8 * (2 * rho) ** (n - 2) <= 1
    assert (n + 1) * (1 + rho) < 2 * n
    assert (n + 1) * rho + 3 * n * (1 + rho) < 4 * n


def test_exact_serialization_and_scope_not_promoted():
    for packet in (
        gaussian.data(),
        gaussian.variation(),
        homogeneous.data(),
        vacuum.data(),
        lorentz.data(),
    ):
        verify.serialize(verify.payload(packet))
    assert all(v is True for v in audit.gates().values())
    rows = copy.deepcopy(audit.matching())
    rows[-1]["status"] = "COMPLETE"
    with pytest.raises(ValueError):
        audit.validate_scope(audit.frontier(), rows)
    with pytest.raises(ValueError):
        verify.serialize(s.Float("0.1"))
