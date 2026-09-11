"""Independent radial, literal retarded-operator and contact comparisons."""

from functools import cache
from itertools import pairwise

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_retarded_state_remainder import (
    audit,
    contact,
    memory,
    response,
    tail,
)


@cache
def gauss(n):
    nodes, weights = mp.gauss_quadrature(n, "legendre")
    return tuple((x + 1) / 2 for x in nodes), tuple(w / 2 for w in weights)


def quad(f, a, b, n=48):
    if a == b:
        return mp.mpf(0)
    nodes, weights = gauss(n)
    return (b - a) * sum(w * f(a + (b - a) * x) for x, w in zip(nodes, weights))


def angular(r, p, K, kind, removed):
    if r == 0:
        return mp.mpf(0)
    nu = mp.sqrt(1 + r * r)
    if p == 0:
        if removed and r <= K:
            return mp.mpf(0)
        return 4 / nu**4 if kind == "linear" else 8 / nu**10
    top = mp.mpf(1)
    if removed and r <= K:
        top = min(top, (p * p + r * r - K * K) / (2 * p * r))
    if top <= -1:
        return mp.mpf(0)
    alpha = 1 + p * p + r * r
    beta = 2 * p * r
    upper = mp.sqrt(alpha + beta)
    lower = mp.sqrt(alpha - beta * top)

    def moment(power):
        return (upper ** (power + 2) - lower ** (power + 2)) / (
            beta * (mp.mpf(power) / 2 + 1)
        )

    if kind == "linear":
        return nu**-5 * moment(1) + nu * moment(-5)
    return nu**-11 * moment(1) + 2 * nu**-5 * moment(-5) + nu * moment(-11)


@pytest.mark.parametrize("p", ("0", "0.5", "2", "6"))
@pytest.mark.parametrize("kind", ("linear", "quadratic"))
@pytest.mark.parametrize("removed", (False, True))
def test_independent_two_momentum_integrals_with_all_external_transfers(
    p, kind, removed
):
    with mp.workdps(70):
        gauss.cache_clear()
        P = mp.mpf(p)
        K = mp.mpf(2)
        cuts = sorted({mp.mpf(0), K, max(mp.mpf(0), min(K, abs(P - K)))})
        value = sum(
            quad(lambda r: r * r * angular(r, P, K, kind, removed), a, b)
            for a, b in pairwise(cuts)
        )
        value += quad(
            lambda x: K**3 / x**4 * angular(K / x, P, K, kind, removed),
            mp.mpf(0),
            mp.mpf(1),
        )
        value /= 4 * mp.pi**2
        power = 4 if kind == "linear" else 10
        factor = 2 if kind == "linear" else 4
        if removed and P <= K / 2:
            R = K / 2
            radial = mp.quad(
                lambda y: y ** (power - 4) / (1 + (y / R) ** 2) ** (mp.mpf(power) / 2),
                [0, 1],
            ) / (2 * mp.pi**2 * R ** (power - 3))
        else:
            radial = 1 / (8 * mp.pi) if kind == "linear" else 5 / (512 * mp.pi)
        assert 0 < value < factor * (1 + P) * radial


@pytest.mark.parametrize("K", (1000, 10**6, 10**16))
def test_independent_complete_and_removed_contact_radial(K):
    with mp.workdps(80):
        A = mp.mpf(25) / 16
        m = mp.mpf(1000)
        R = mp.mpf(K)
        complete = (
            A**3
            * mp.quad(lambda y: y * y / (1 + y * y) ** mp.mpf("2.5"), [0, 1, mp.inf])
            / (2 * mp.pi**2 * m * m)
        )
        exact = A**3 / (6 * mp.pi**2 * m * m)
        assert abs(complete - exact) < mp.mpf("1e-70")
        removed = (
            A**5
            * mp.quad(lambda y: y / (1 + (A * m * y / R) ** 2) ** mp.mpf("2.5"), [0, 1])
            / (2 * mp.pi**2 * R * R)
        )
        assert 0 < removed < A**5 / (4 * mp.pi**2 * R * R)
        assert removed * mp.mpf(str(contact.MODE)) < mp.mpf("1e14") / R**2


def dagger(A):
    return A.transpose_conj()


def tr(A):
    return sum(A[i, i] for i in range(A.rows))


def kron(A, B):
    return mp.matrix(
        [
            [
                A[i // B.rows, j // B.cols] * B[i % B.rows, j % B.cols]
                for j in range(A.cols * B.cols)
            ]
            for i in range(A.rows * B.rows)
        ]
    )


def literal_fixture(beta):
    a = mp.zeros(3)
    a[0, 1] = 1
    a[1, 2] = mp.sqrt(2)
    annihilators = [kron(a, mp.eye(3)), kron(mp.eye(3), a)]
    D = mp.matrix([[1, mp.mpf(2) / 3], [mp.mpf(2) / 3, -1]])
    G = mp.matrix([[mp.mpf(1) / 2, -mp.mpf(1) / 3], [-mp.mpf(1) / 3, -mp.mpf(1) / 2]])
    frequencies = [mp.mpf(3), mp.mpf(7)]
    betas = [beta, -beta / 2]

    def values(t, reference):
        out = []
        for w, b in zip(frequencies, betas):
            alpha = mp.sqrt(1 + b * b)
            f = mp.exp(-1j * w * t) / mp.sqrt(2 * w)
            out.append(alpha * f if reference else alpha * f + b * mp.conj(f))
        return out

    def observable(t, M, reference):
        u = values(t, reference)
        q = [
            u[j] * annihilators[j] + mp.conj(u[j]) * dagger(annihilators[j])
            for j in range(2)
        ]
        T = mp.zeros(9)
        for i in range(2):
            for j in range(2):
                T += M[i, j] * q[i] * q[j]
        return T

    return D, G, values, observable


@pytest.mark.parametrize("beta", ("0.01", "0.2", "0.5"))
def test_literal_retarded_Fock_comparison_keeps_both_mixed_terms_and_error_square(beta):
    with mp.workdps(65):
        gauss.cache_clear()
        D, G, values, observable = literal_fixture(mp.mpf(beta))
        T = mp.mpf(1) / 4
        nodes, weights = gauss(8)
        literal = mp.mpf(0)
        pair = mp.mpf(0)
        mixed_only = mp.mpf(0)
        for x, wx in zip(nodes, weights):
            t = T * x
            for y, wy in zip(nodes, weights):
                tau = t * y
                w = wx * wy * T * T * x * x * (1 - x) * (x * y) * (1 - x * y)
                od = observable(t, D, False)
                og = observable(tau, G, False)
                rd = observable(t, D, True)
                rg = observable(tau, G, True)
                literal += (
                    w * 1j * ((od * og - og * od)[0, 0] - (rd * rg - rg * rd)[0, 0]) / 4
                )
                ud, ug = values(t, False), values(tau, False)
                vd, vg = values(t, True), values(tau, True)
                for i in range(2):
                    for j in range(2):
                        # The creation coefficient is conjugate(u_i u_j).
                        ad = D[i, j] * mp.conj(ud[i] * ud[j])
                        ag = G[i, j] * mp.conj(ug[i] * ug[j])
                        bd = D[i, j] * mp.conj(vd[i] * vd[j])
                        bg = G[i, j] * mp.conj(vg[i] * vg[j])
                        ed, eg = ad - bd, ag - bg
                        diff = mp.conj(bd) * eg + mp.conj(ed) * bg + mp.conj(ed) * eg
                        pair -= w * mp.im(diff)
                        mixed_only -= w * mp.im(mp.conj(bd) * eg + mp.conj(ed) * bg)
        assert abs(literal - pair) < mp.mpf("1e-55")
        assert abs(literal) > mp.mpf("1e-12")
        assert abs(pair - mixed_only) > mp.mpf("1e-15")


@pytest.mark.parametrize("nu", (1000, 10**6, 10**16))
@pytest.mark.parametrize("fixture", (0, 1, 2))
def test_literal_ten_field_three_mode_contact_difference(nu, fixture):
    with mp.workdps(170):
        N = mp.mpf(nu)
        F = mp.mpf(4000)
        R = mp.mpf(2) * 10**6
        D = mp.matrix(
            [
                [1, mp.mpf(2) / 3, -mp.mpf(1) / 5],
                [mp.mpf(2) / 3, -2, mp.mpf(1) / 7],
                [-mp.mpf(1) / 5, mp.mpf(1) / 7, 1],
            ]
        )
        G = mp.matrix(
            [
                [2, -mp.mpf(1) / 3, mp.mpf(1) / 2],
                [-mp.mpf(1) / 3, -3, mp.mpf(2) / 5],
                [mp.mpf(1) / 2, mp.mpf(2) / 5, 1],
            ]
        )
        H = (D * G + G * D) / 2
        S = mp.zeros(10)
        for block in range(3):
            for i in range(3):
                for j in range(3):
                    S[3 * block + i, 3 * block + j] = H[i, j]
        direct = mp.mpf(0)
        expanded = mp.mpf(0)
        for mode in range(3):
            v = mp.matrix(
                [
                    mp.mpc(
                        (j + 1) * (mode + 1) + fixture, (-1) ** j * (j + 2 + fixture)
                    )
                    for j in range(10)
                ]
            )
            e = mp.matrix(
                [
                    mp.mpc((-1) ** (j + mode) * (j + 2), 2 * j + 1 + fixture)
                    for j in range(10)
                ]
            )
            v *= F * mp.sqrt(N) / mp.sqrt((dagger(v) * v)[0])
            e *= F * R / N ** mp.mpf("5.5") / mp.sqrt((dagger(e) * e)[0])
            direct += ((dagger(v + e) * S * (v + e))[0] - (dagger(v) * S * v)[0]) / 2
            expanded += (
                (dagger(v) * S * e)[0] + (dagger(e) * S * v)[0] + (dagger(e) * S * e)[0]
            ) / 2
            # V0 is retained in the readout but has no direct shear contact.
            assert v[9] != 0 and e[9] != 0 and S[9, 9] == 0
        scale = mp.mpf(str(contact.MODE)) / N**5 * mp.norm(D) * mp.norm(G)
        assert abs(direct - expanded) < scale * mp.mpf("1e-65")
        assert 0 < abs(direct) < scale
        assert D * G != G * D


@pytest.mark.parametrize(
    "beta", (s.Rational(1, 100), s.Rational(1, 5), s.Rational(1, 2))
)
def test_constant_alpha_reference_is_not_a_new_CCR_normalized_state(beta):
    t = s.Symbol("t", real=True)
    w = s.Integer(3)
    alpha = s.sqrt(1 + beta**2)
    f = s.exp(-s.I * w * t) / s.sqrt(2 * w)
    actual = alpha * f + beta * s.conjugate(f)
    ref = alpha * f
    wronskian = lambda u: s.simplify(
        u * s.conjugate(s.diff(u, t)) - s.diff(u, t) * s.conjugate(u)
    )
    assert wronskian(actual) == s.I
    assert wronskian(ref) == s.I * (1 + beta**2)
    assert wronskian(ref) != s.I


@pytest.mark.parametrize("K", (1000, 10**6, 10**16))
def test_complete_time_ordered_and_contact_tail_displays(K):
    c = tail.constants()
    assert c["all_external_memory_tail_numerator"] / K + c[
        "contact_K_squared_tail_numerator"
    ] / K**2 < tail.error_bound(K)
    assert (
        4 * tail.error_bound(K) / s.Integer(10) ** 800
        == 4 * s.Rational(1, 10) ** 773 / K
    )
    assert 4 * response.DISPLAY / s.Integer(10) ** 800 == response.CANONICAL


def test_low_external_tail_requires_both_legs_and_keeps_high_transfer_contact():
    K = s.Integer(1000)
    k = 2 * K
    l = 0
    assert max(abs(k), abs(l)) > K and abs(k + l) > K / 2 and abs(l) <= K / 2
    # Pointwise Gamma=D*cos(px) gives zero mean first harmonic but
    # a nonzero second contact with spatial average D^2/2.
    D = s.diag(1, -1, 0)
    x = s.Symbol("x", real=True)
    H = D * D * s.cos(x) ** 2
    mean = H.applyfunc(lambda v: s.integrate(v, (x, 0, 2 * s.pi)) / (2 * s.pi))
    assert mean == D * D / 2 and s.trace(mean) > 0


def test_reference_remainder_does_not_claim_the_singular_extension():
    assert "not a new physical state" in memory.data()["comparison"]
    assert "not CCR-normalized" in response.data()["reference_not_a_state"]
    assert (
        "ORIGINAL covariant retarded prescription"
        in response.data()["matching_boundary"]
    )
    assert "theta(t-s)" in response.data()["full_comparison"]


@pytest.mark.parametrize("name", ("memory", "contact", "tail", "response"))
def test_exact_scientific_packet(name):
    p = {"memory": memory, "contact": contact, "tail": tail, "response": response}[
        name
    ].data()
    for v in p["checks"].values():
        entries = list(v) if isinstance(v, s.MatrixBase) else [v]
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
    "name,call,args", audit.bad_cases(), ids=lambda x: x if isinstance(x, str) else None
)
def test_all_scope_mutations_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_complete_counts_and_unclosed_frontier():
    assert len(audit.residuals()) == 28 and audit.scalar_entry_count() == 51
    assert len(audit.gates()) == 38 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 119
    assert len(audit.frontier()) == 9 and audit.matching()[-1] == audit.ITEM
