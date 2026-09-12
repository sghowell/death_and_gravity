"""Independent literal shear, complete primitive, weak-log and causal-order checks."""

import pytest
import sympy as s
from p8_vacuum_affine_homogeneous_shear_inverse import audit, local, matching, volterra
from p8_vacuum_affine_local_tensor_response import operator as loc


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_all_exact_residuals(name, value):
    expected = s.zeros(*value.shape) if isinstance(value, s.MatrixBase) else s.S.Zero
    assert value == expected, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_all_proof_gates(name, value):
    assert value is True, name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_unchanged_original_scope_and_specific_homogeneous_boundary():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 84
    assert audit.frontier() == audit.previous.frontier()
    assert "HOMOGENEOUS" in audit.ITEM["status"]
    assert "NOT_NONZERO_TRANSFER" in audit.ITEM["status"]
    assert all(row["status"] != "COMPLETE" for row in audit.frontier())


def test_independent_matching_local_tree_and_primitives():
    d = s.Symbol("d", positive=True)
    # Independent rotational moments, with tr(D)=0 and tr(D^2)=1.
    TT = 1 - 2 / d + (d - 1) / (2 * d * (d + 2))
    LL = 1 / (2 * d * (d + 2))
    C = s.factor(TT + LL)
    assert C == (2 * d * d + d - 8) / (2 * d * (d + 2))
    assert TT.subs(d, 3) == s.Rational(2, 5) and LL.subs(d, 3) == s.Rational(1, 30)
    assert C.subs(d, 3) == s.Rational(13, 30) and s.diff(C, d).subs(d, 3) == s.Rational(
        91, 450
    )
    radial = s.factor(32 * C.subs(d, 3))
    lag, eps = s.symbols("lag eps", positive=True)
    abel = s.factor(s.im(24 / (eps - 2 * s.I * lag) ** 5))
    assert s.limit(abel, eps, 0, dir="+") == s.Rational(3, 4) / lag**5
    assert s.factor(radial * s.Rational(3, 4) - s.Rational(52, 5)) == 0
    assert s.diff(s.Rational(13, 30) / lag, lag, 4) == s.Rational(52, 5) / lag**5
    a0, p = s.symbols("a0 p", positive=True)
    assert s.simplify(a0 ** (-d - 2) * (a0 * p) ** (d + 1) * a0 - p ** (d + 1)) == 0
    # Independent first lag jet, with time-varying scale and its original output volume.
    H = s.Symbol("H", real=True)
    a_t = a0 * (1 + H * lag)
    delta = lag / a0 * (1 - H * lag / 2)
    ratio = s.series(lag**5 / (a_t**4 * a0 * delta**5), lag, 0, 2).removeO()
    assert s.expand(ratio - 1 + s.Rational(3, 2) * H * lag) == 0
    # Original q0 finite local Euler coefficients and canonical tree normalization.
    t = loc.t
    f = s.Function("f")(t)
    A, H = loc.A, loc.H
    L = s.diff(f, t, 2) + 3 * H * s.diff(f, t)
    Dt = s.diff(f, t, 2) + H * s.diff(f, t)
    Dstar = lambda v: (
        s.diff(v, t, 2) + 5 * H * s.diff(v, t) + (2 * s.diff(H, t) + 6 * H**2) * v
    )
    Tloc = -2 * A * L - 2 * s.diff(A, t) * s.diff(f, t) - Dstar(Dt) / 30
    Ttree = -16 * s.pi**2 * loc.KAPPA * L
    total = s.expand(Tloc + Ttree)
    expected = {
        4: -s.Rational(1, 30),
        3: -H / 5,
        2: -2 * A - (4 * s.diff(H, t) + 11 * H * H) / 30 - 16 * s.pi**2 * loc.KAPPA,
        1: -6 * H * A
        - 2 * s.diff(A, t)
        - (s.diff(H, t, 2) + 7 * H * s.diff(H, t) + 6 * H**3) / 30
        - 48 * s.pi**2 * loc.KAPPA * H,
        0: s.S.Zero,
    }
    for j, v in expected.items():
        assert s.factor(total.coeff(s.diff(f, t, j)) - v) == 0
    Qloc = (A * L + s.diff(A, t) * s.diff(f, t) + Dstar(Dt) / 60) / (
        8 * s.pi**2 * loc.KAPPA
    )
    assert s.factor(Tloc + 16 * s.pi**2 * loc.KAPPA * Qloc) == 0
    # Independent polynomial integration for every lower variable coefficient term.
    x, y = s.symbols("x y", nonnegative=True)
    source = y**6 * (1 + y + y**2)
    coef = 1 + 2 * y + 3 * y * y + y**4
    for j in range(4):
        raw = s.integrate((x - y) ** 3 * s.diff(source, y, j) * coef / 6, (y, 0, x))
        ker = sum(
            (-1) ** r
            * s.binomial(j, r)
            * (x - y) ** (3 - j + r)
            * s.diff(coef, y, r)
            / s.factorial(3 - j + r)
            for r in range(j + 1)
        )
        rebuilt = s.integrate(ker * source, (y, 0, x))
        assert s.expand(raw - rebuilt) == 0
        if j:
            bad = s.integrate(
                (x - y) ** (3 - j) * coef * source / s.factorial(3 - j), (y, 0, x)
            )
            assert s.expand(raw - bad) != 0
    # Exact finite-weight contraction, with C including all tree and quantum terms.
    K, Cb = s.symbols("K Cb", positive=True)
    beta = K * Cb
    assert (
        s.factor(
            s.Rational(1, 2) - 2 * beta / (4 * beta + 1) - 1 / (2 * (4 * beta + 1))
        )
        == 0
    )
    assert s.factor(2 * beta / (4 * beta + 1)).is_positive


@pytest.mark.parametrize("dimension", (3, 4, 5, 6))
def test_independent_full_transverse_polarization_square(dimension):
    n = s.zeros(dimension, 1)
    n[0], n[1], n[2] = s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3)
    P = s.eye(dimension) - n * n.T
    D = s.diag(1, -1, *([0] * (dimension - 2)))
    G = s.diag(2, 1, -3, *([0] * (dimension - 3)))
    D[0, 1] = D[1, 0] = s.Rational(1, 3)
    G[0, 2] = G[2, 0] = s.Rational(2, 5)
    dn, gn = (n.T * D * n)[0], (n.T * G * n)[0]
    TD = P * D * P + dn * P / 2
    TG = P * G * P + gn * P / 2
    complete = s.trace(TD * TG) + dn * gn / 4
    independent = (
        s.trace(D * G) - 2 * (n.T * D * G * n)[0] + s.Rational(dimension, 4) * dn * gn
    )
    assert s.factor(complete - independent) == 0
    assert s.factor(dn * gn / 4) != 0
    assert P * P == P and s.trace(P) == dimension - 1


@pytest.mark.parametrize("power", (0, 5))
def test_weak_log_simultaneous_derivative_and_initial_boundary(power):
    x = s.Symbol("positive_output", positive=True)
    I = lambda n: x ** (n + 1) * (s.log(x) - s.harmonic(n + 1)) / (n + 1)
    original = (1 + x) * I(power) + I(power + 1)
    source_derivative = 0 if power == 0 else power * ((1 + x) * I(power - 1) + I(power))
    diagonal_derivative = 2 * I(power)
    initial = (1 + x) * s.log(x) if power == 0 else s.S.Zero
    assert (
        s.simplify(
            s.diff(original, x) - source_derivative - diagonal_derivative - initial
        )
        == 0
    )
    assert s.simplify(s.diff(original, x) - source_derivative - initial) != 0
    if power == 0:
        assert (
            s.simplify(s.diff(original, x) - source_derivative - diagonal_derivative)
            != 0
        )


@pytest.mark.parametrize("variant", (1, 2, 3))
def test_independent_full_noncommuting_causal_normal_form_inverse(variant):
    v = s.Integer(variant)
    K = s.Matrix(
        [
            [s.Rational(1, 5), 0, 0],
            [s.Rational(1, 7), s.Rational(2, 5), 0],
            [v / 9, s.Rational(1, 6), s.Rational(3, 5)],
        ]
    )
    V = s.Matrix(
        [
            [s.Rational(1, 3), 0, 0],
            [s.Rational(1, 8), -v / 4, 0],
            [s.Rational(1, 9), s.Rational(1, 11), s.Rational(1, 5)],
        ]
    )
    D = s.Matrix([[1, 0, 0], [-1, 1, 0], [0, -1, 1]])
    I = s.eye(3)
    I4 = D.inv() ** 4
    E = (I + K * V).inv() * K
    F = K.inv()
    T = D**4 * (F + V)
    inverse = E * I4
    assert (F + V) * E == I and E * (F + V) == I
    assert T * inverse == I and inverse * T == I
    assert I4 * T == F + V
    assert K * V != V * K
    assert K * (I + K * V).inv() != E
    assert I4 * E != inverse
    assert K * I4 != inverse


@pytest.mark.parametrize(
    "kernel,majorant",
    (
        (1, 1),
        (2, 3),
        (s.Rational(1, 7), s.Rational(3, 5)),
        (10**4, 10**8),
        (10**3, 10**810),
    ),
)
def test_finite_weight_keeps_arbitrary_actual_majorant_scale(kernel, majorant):
    beta = (
        s.Integer(kernel) * majorant if isinstance(kernel, int) else kernel * majorant
    )
    weight = (4 * beta + 1) ** 2
    upper = 2 * beta / (4 * beta + 1)
    assert weight > 1 and 0 < upper < s.Rational(1, 2)
    assert s.Rational(1, 2) - upper == 1 / (2 * (4 * beta + 1))


@pytest.mark.parametrize("module", (matching, local, volterra))
def test_module_checks_are_original_exact_arithmetic(module):
    for value in module.data()["checks"].values():
        assert not value.atoms(s.Float)
        value = (
            value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
        )
        assert (
            all(v == 0 for v in value)
            if isinstance(value, s.MatrixBase)
            else value == 0
        )
