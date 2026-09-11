"""Independent all-momentum tail, Fock projection and support checks."""

from functools import cache
from itertools import pairwise, product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_separated_response import (
    audit,
    projection,
    response,
    support,
    tail,
)


def trace(M):
    return sum(M[i, i] for i in range(M.rows))


def quad_gauss(f, a, b, n=48):
    if a == b:
        return mp.mpf(0)
    nodes, weights = mp.gauss_quadrature(n, "legendre")
    return (
        (b - a)
        / 2
        * sum(w * f((a + b + (b - a) * x) / 2) for x, w in zip(nodes, weights))
    )


@pytest.mark.parametrize("power", (4, 10))
@pytest.mark.parametrize("K", (1000, 10**6, 10**16))
def test_independent_complete_massive_radial_tail(power, K):
    with mp.workdps(70):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        R = mp.mpf(K) if power == 4 else mp.mpf(K) / 2
        alpha = A * m / R
        integral = (
            A**power
            * R ** (3 - power)
            * mp.quad(
                lambda y: (
                    y ** (power - 4) / (1 + (alpha * y) ** 2) ** (mp.mpf(power) / 2)
                ),
                [0, 1],
            )
            / (2 * mp.pi**2)
        )
        exact_majorant = A**power * R ** (3 - power) / (2 * mp.pi**2 * (power - 3))
        assert 0 < integral < exact_majorant
        if power == 10:
            assert integral < 100 / mp.mpf(K) ** 7


def angular_pair_integral(r, p, K, kind):
    if r == 0:
        return mp.mpf(0)
    nu = mp.sqrt(1 + r * r)
    if p == 0:
        if r <= K:
            return mp.mpf(0)
        return 1 / (32 * nu**4) if kind == "reference" else 8 / nu**10
    top = mp.mpf(1) if r > K else min(mp.mpf(1), (p * p + r * r - K * K) / (2 * p * r))
    if top <= -1:
        return mp.mpf(0)
    alpha = 1 + p * p + r * r
    beta = 2 * p * r
    upper = mp.sqrt(alpha + beta)
    lower = mp.sqrt(alpha - beta * top)
    if kind == "reference":
        primitive = lambda z: (
            -1 / (3 * (1 + z) ** 3) + 1 / (2 * (1 + z) ** 4) - 1 / (5 * (1 + z) ** 5)
        )
        return 2 / (beta * nu**2) * (primitive(upper / nu) - primitive(lower / nu))
    moment = lambda exponent: (
        (upper ** (exponent + 2) - lower ** (exponent + 2))
        / (beta * (mp.mpf(exponent) / 2 + 1))
    )
    return nu**-11 * moment(1) + 2 * nu**-5 * moment(-5) + nu * moment(-11)


@pytest.mark.parametrize("p", ("0", "0.5", "2", "6"))
@pytest.mark.parametrize("kind", ("reference", "error"))
def test_independent_nonzero_external_momentum_union_tail(p, kind):
    # Dimensionless radial/angular integral; A^3/m or A^3/m^7
    # cancels from both sides of the corresponding physical inequality.
    with mp.workdps(65):
        P = mp.mpf(p)
        K = mp.mpf(2)
        cuts = sorted({mp.mpf(0), K, max(mp.mpf(0), min(K, abs(P - K)))})
        value = mp.mpf(0)
        for a, b in pairwise(cuts):
            value += quad_gauss(
                lambda r: r * r * angular_pair_integral(r, P, K, kind), a, b
            )
        value += quad_gauss(
            lambda x: K**3 / x**4 * angular_pair_integral(K / x, P, K, kind),
            mp.mpf(0),
            mp.mpf(1),
        )
        value /= 4 * mp.pi**2
        if kind == "reference":
            radial = mp.quad(lambda y: 1 / (1 + (y / K) ** 2) ** 2, [0, 1]) / (
                2 * mp.pi**2 * K
            )
            bound = radial / 2
        elif P <= K / 2:
            R = K / 2
            radial = mp.quad(lambda y: y**6 / (1 + (y / R) ** 2) ** 5, [0, 1]) / (
                2 * mp.pi**2 * R**7
            )
            bound = 4 * (1 + P) * radial
        else:
            bound = 4 * (1 + P) * 5 / (512 * mp.pi)
        assert value > 0 and value < bound


@pytest.mark.parametrize("K", (1000, 2000, 10**16))
@pytest.mark.parametrize("ratio", (s.Rational(501, 1000), 1, 10))
def test_high_external_spatial_H1_weight_without_a_band(K, ratio):
    P = s.Integer(K) * ratio
    A = s.Rational(25, 16)
    m = s.Integer(1000)
    L = 1 + P / (A * m)
    assert L / (1 + P * P) <= 4 / s.Integer(K) ** 2 + 2 / (A * m * K)
    assert 4 / s.Integer(K) ** 2 + 2 / (A * m * K) < s.Rational(1, 100 * K)


def test_low_external_both_leg_implication_and_its_missing_hypothesis():
    K = s.Integer(1000)
    for k, l in ((1500, -1100), (1100, -600), (1001, -501)):
        assert abs(k) > K and abs(k + l) <= K / 2
        assert abs(k) > K / 2 and abs(l) > K / 2
    k, l = 2 * K, 0
    assert max(abs(k), abs(l)) > K and abs(k + l) > K / 2
    assert not (abs(k) > K / 2 and abs(l) > K / 2)


@cache
def fock_fixture():
    count = 3
    states = list(product(range(3), repeat=count))
    dim = len(states)
    a = s.zeros(3)
    a[0, 1] = 1
    a[1, 2] = s.sqrt(2)
    annihilators = []
    for mode in range(count):
        annihilators.append(
            s.kronecker_product(*[a if j == mode else s.eye(3) for j in range(count)])
        )
    creators = [A.T for A in annihilators]
    B = (
        s.Matrix(
            [
                [1 + s.I, 2 - s.I, -1 + 2 * s.I],
                [2 - s.I, -1 + s.I, 3],
                [-1 + 2 * s.I, 3, 2 - s.I],
            ]
        )
        / 17
    )
    C = (
        s.Matrix(
            [
                [2 - s.I, -1 + 2 * s.I, 1],
                [-1 + 2 * s.I, 1 + s.I, -2 - s.I],
                [1, -2 - s.I, 3 + 2 * s.I],
            ]
        )
        / 19
    )

    def observable(pair):
        result = s.zeros(dim)
        for i in range(count):
            result += (i + 1) * creators[i] * annihilators[i] / 7
            for j in range(count):
                result += (
                    pair[i, j] * creators[i] * creators[j]
                    + s.conjugate(pair[i, j]) * annihilators[i] * annihilators[j]
                )
        return result.applyfunc(s.expand)

    O, P = observable(B), observable(C)
    vacuum = s.eye(dim)[:, 0]
    return states, B, C, O, P, vacuum


@pytest.mark.parametrize("retained", ((0,), (0, 2), (0, 1, 2)))
def test_literal_Fock_current_and_common_two_particle_projection(retained):
    states, B, C, O, P, vac = fock_fixture()
    x, y = O * vac, P * vac
    inner = lambda u, v: s.expand((u.conjugate().T * v)[0])
    projector = s.diag(
        *[
            int(
                sum(state) == 2
                and all(state[j] == 0 for j in range(3) if j not in retained)
            )
            for state in states
        ]
    )
    z = inner(x, y)
    pair_inner = 2 * sum(
        s.conjugate(B[i, j]) * C[i, j] for i in range(3) for j in range(3)
    )
    assert s.expand(z - pair_inner) == 0
    assert s.expand(inner(x, x) - 2 * s.trace(B.conjugate().T * B)) == 0
    current = s.I * (vac.T * (O * P - P * O) * vac)[0] / 4
    assert s.simplify(current + s.im(z) / 2) == 0
    remainder = z - inner(projector * x, projector * y)
    assert (
        s.expand(
            remainder
            - inner(
                (s.eye(len(states)) - projector) * x,
                (s.eye(len(states)) - projector) * y,
            )
        )
        == 0
    )
    assert s.expand(inner(projector * x, (s.eye(len(states)) - projector) * y)) == 0
    if len(retained) < 3:
        assert (
            inner(
                (s.eye(len(states)) - projector) * x,
                (s.eye(len(states)) - projector) * x,
            )
            > 0
        )


def test_pair_exchange_factor_and_second_leg_projection_cannot_be_omitted():
    _, B, C, O, _P, vac = fock_fixture()
    norm = s.expand(((O * vac).conjugate().T * (O * vac))[0])
    assert s.simplify(norm - s.trace(B.conjugate().T * B)) > 0
    P = s.diag(1, 0, 1)
    assert P * B != P * B * P
    correct = 2 * s.trace((P * B * P).conjugate().T * (P * C * P))
    one_leg = 2 * s.trace((P * B).conjugate().T * (P * C))
    assert s.expand(correct - one_leg) != 0


@pytest.mark.parametrize(
    "intervals",
    (
        (-s.Rational(3, 8), -s.Rational(1, 8), s.Rational(1, 8), s.Rational(3, 8)),
        (
            -s.Rational(49, 100),
            -s.Rational(48, 100),
            s.Rational(48, 100),
            s.Rational(49, 100),
        ),
        (
            -s.Rational(1, 4),
            -s.Rational(1, 1000),
            s.Rational(1, 1000),
            s.Rational(1, 4),
        ),
    ),
)
def test_strict_time_ordered_domain(intervals):
    assert support.require_ordered(*intervals) == intervals


@pytest.mark.parametrize(
    "intervals",
    (
        (0, s.Rational(1, 4), s.Rational(1, 8), s.Rational(3, 8)),
        (-s.Rational(1, 4), 0, 0, s.Rational(1, 4)),
        (-s.Rational(1, 2), -s.Rational(1, 4), 0, s.Rational(1, 4)),
        (0, s.Rational(1, 8), s.Rational(1, 4), s.Rational(1, 2)),
    ),
)
def test_overlap_touching_and_slab_endpoints_are_rejected(intervals):
    with pytest.raises(ValueError, match="strict"):
        support.require_ordered(*intervals)


@pytest.mark.parametrize("values", ((3, 1, 1), (5, 2, 1), (s.Rational(5, 2), 1, 1)))
def test_actual_clock_spacelike_enclosure(values):
    assert support.require_spacelike_enclosure(*values) > 0


@pytest.mark.parametrize(
    "values", ((2, 1, 1), (1, 1, 1), (-1, 0, 0), (3, -1, 1), (3, 1, -1))
)
def test_nonspacelike_or_invalid_enclosures_rejected(values):
    with pytest.raises(ValueError, match="spacelike"):
        support.require_spacelike_enclosure(*values)


def test_overlapping_time_order_cannot_be_replaced_by_full_commutator():
    with mp.workdps(60):
        T = mp.mpf(1) / 4
        bump = lambda x: mp.exp(-1 / (x * (1 - x))) if 0 < x < 1 else mp.mpf(0)
        # For a single oscillator q^2 in its vacuum, the i/4 commutator
        # kernel is sin(2(t-s))/4. Its full square pairing with the same
        # real test is zero by antisymmetry, but its retarded triangle
        # pairing is strictly positive on this short time interval.
        nodes, weights = mp.gauss_quadrature(20, "legendre")
        nodes = [(x + 1) / 2 for x in nodes]
        weights = [w / 2 for w in weights]
        square = (
            T
            * T
            * sum(
                w * v * bump(x) * bump(y) * mp.sin(2 * T * (x - y)) / 4
                for x, w in zip(nodes, weights)
                for y, v in zip(nodes, weights)
            )
        )
        triangle = (
            T
            * T
            * sum(
                w * v * x * bump(x) * bump(x * y) * mp.sin(2 * T * x * (1 - y)) / 4
                for x, w in zip(nodes, weights)
                for y, v in zip(nodes, weights)
            )
        )
        assert abs(square) < mp.mpf("1e-55")
        assert triangle > mp.mpf("1e-12")
        assert 2 * T < mp.pi


@pytest.mark.parametrize("K", (1000, 10**6, 10**16))
def test_complete_tail_current_and_two_canonical_factor_displays(K):
    variance = tail.variance_tail_display(K)
    assert variance == s.Integer(10) ** 52 / K
    assert 4 * (variance / 2) / s.Integer(10) ** 800 == 2 * s.Rational(1, 10) ** 748 / K
    assert response.CANONICAL == 2 * s.Rational(1, 10) ** 750


@pytest.mark.parametrize("name", ("tail", "projection", "support", "response"))
def test_exact_scientific_packet(name):
    p = {
        "tail": tail,
        "projection": projection,
        "support": support,
        "response": response,
    }[name].data()
    for value in p["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(x) == 0 for x in entries)
    assert all(p["gates"].values())


@pytest.mark.parametrize(
    "name,value",
    audit.residuals().items(),
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_all_exact_residuals(name, value):
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.cancel(x) == 0 for x in entries), name


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_all_scope_mutations_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_complete_counts_and_unclosed_frontier():
    assert len(audit.residuals()) == 32 and audit.scalar_entry_count() == 81
    assert len(audit.gates()) == 37 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 118
    assert len(audit.frontier()) == 9 and audit.matching()[-1] == audit.ITEM


def test_weak_domain_and_coincident_boundary_are_explicit():
    assert "operator-norm differentiability" in projection.data()["boundary"]
    assert "No coincident-support extension" in support.data()["boundary"]
    assert (
        "not a fully reduced mixed scalar/metric norm"
        in response.data()["normalization"]
    )
    assert "No external momentum cutoff" in tail.data()["high_external"]
