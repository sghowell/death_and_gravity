"""Independent full graph, dimensional, cut and analytic-majorant tests."""

import math

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from numpy.polynomial.legendre import leggauss
from p8_vacuum_affine_whole_mixed_heavy_gravity_sector import (
    audit,
    cuts,
    forward,
    graphs,
    source,
)
from scipy.special import gamma

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_scope_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


class IndependentWhole:
    def __init__(self, N, mu, heavy):
        nodes, weights = leggauss(N)
        self.x = (nodes + 1) / 2
        self.w = weights / 2
        self.mu = float(mu)
        self.n = float(heavy)
        self.cache = {}

    def moments(self, a, e):
        key = (a, e)
        if key in self.cache:
            return self.cache[key]
        x, w, mu, n = self.x, self.w, self.mu, self.n
        A = mu - a * x * (1 - x)
        Mll = np.dot(w, A**e)
        J = np.dot(w, A ** (e - 1))
        J1 = np.dot(w, np.log(A) / A)
        Lll = np.dot(w, np.log(A))
        mix = n * x + mu * (1 - x) ** 2
        Mmix = np.dot(w, mix**e)
        Lmix = np.dot(w, np.log(mix))
        # v=x^4 resolves the logarithmic massless bubble endpoint.
        v = x**4
        Mh = np.dot(w, 4 * x**3 * (v * (n - a * (1 - v))) ** e)
        Lh = -1 + np.dot(w, np.log(n - a * (1 - x)))
        # The active-light massive triangle, including its full e dependence.
        u = x[:, None]
        vv = x[None, :]
        A2 = mu - a * vv * (1 - vv)
        C = u * u * A2 + n * (1 - u)
        I = np.sum(w[:, None] * w[None, :] * u * C ** (e - 1))
        # Both parameters fourth-powered for the integrable off-shell cusp corner.
        r = x[:, None] ** 4
        vv = x[None, :] ** 4
        den = (n - a) * vv + r * (a * vv + mu * (1 - vv) ** 2)
        Uint = np.sum(
            w[:, None]
            * w[None, :]
            * 16
            * x[:, None] ** 3
            * x[None, :] ** 3
            * r**e
            * den ** (e - 1)
        )
        out = {
            "Mll": Mll,
            "J": J,
            "J1": J1,
            "Lll": Lll,
            "Mmix": Mmix,
            "Lmix": Lmix,
            "Mh": Mh,
            "Lh": Lh,
            "I": I,
            "U": Uint,
        }
        self.cache[key] = out
        return out

    def box(self, a, b, e):
        key = ("box", a, b, e)
        if key in self.cache:
            return self.cache[key]
        x, w, mu, n = self.x, self.w, self.mu, self.n
        v = x[:, None]
        ww = x[None, :]
        vw = w[:, None] * w[None, :]
        A = mu - b * v * (1 - v)
        q = n - a
        totalK = 0.0
        direct_finite = 0.0
        for xx, weight in zip(x, w):
            u = xx * xx
            jac = 2 * xx
            C = u * u * A
            zmax = q * (1 - u) / (C + q * (1 - u))
            z = zmax * ww
            h = C * z / (q * (1 - z))
            Rminus = a * u * z / q + a * C * z * z / (q * q * (1 - z))
            R = 1 + Rminus
            integrand = (
                u
                * (u + 2 * h)
                * C ** (e - 1)
                / q
                * (1 - z) ** (-e)
                * R ** (e - 2)
                * zmax
            )
            totalK += weight * jac * np.sum(vw * integrand)
            # Independently subtract only the literal box's u=0 endpoint value.
            difference = np.expm1((e - 2) * np.log1p(Rminus))
            inner = np.sum(w[None, :] * (1 - z) ** (-e) * difference * zmax, axis=1)
            end = (1 - zmax[:, 0]) ** (1 - e) / (1 - e)
            direct_finite += (
                weight
                * jac
                * u ** (2 * e - 1)
                * np.dot(w, A[:, 0] ** (e - 1) * (inner - end))
            )
        out = (totalK, direct_finite)
        self.cache[key] = out
        return out

    def whole_raw(self, s, t, e, nu2):
        mu, n = self.mu, self.n
        rows = (s, t, 4 * mu - s - t)
        dim = 4 + 2 * e
        pref = gamma(1 - e) * (4 * math.pi * nu2) ** (-e)
        Bon = -pref * mu**e / (e * (1 + 2 * e))
        Am = -pref * mu ** (1 + e) / (e * (1 + e))
        An = -pref * n ** (1 + e) / (e * (1 + e))
        Sprime = 2 * mu * (dim - 1) * Bon / (dim - 2)
        total = 0.0
        max_box_difference = 0.0
        for a in rows:
            m = self.moments(a, e)
            Va = (a - 2 * mu) ** 2 - 4 * mu * mu / (dim - 2)
            LL = -pref * m["J"] / (2 * e)
            U = -pref * m["U"]
            Bll = -pref * m["Mll"] / e
            Bmix = -pref * m["Mmix"] / e
            Bh = -pref * m["Mh"] / e
            P = (
                -Va * LL
                - 2 * n * (a - 4 * mu / (dim - 2)) * U
                - (a - 2 * mu) * Bll
                + 2 * n * Bmix
                - 4 * mu * (dim - 4) * Bon / (dim - 2)
                + (-a - n + 4 * n / (dim - 2)) * Bh
                + 2 * Am
                + An
            )
            SH = (4 * n * a - 4 * n * n / (dim - 2)) * Bh - 2 * n * An
            total += (
                (2 * P - 2 * Sprime) / (n - a)
                + SH / (n - a) ** 2
                - 2 * dim * Bh / (dim - 2)
            )
        for i, a in enumerate(rows):
            for j, b in enumerate(rows):
                if i == j:
                    continue
                mb = self.moments(b, e)
                ma = self.moments(a, e)
                K, direct = self.box(a, b, e)
                direct_finite = (1 - e) * direct
                ibp_finite = -mb["I"] - (1 - e) * a * K
                max_box_difference = max(
                    max_box_difference, abs(direct_finite - ibp_finite)
                )
                Dbox = pref / (n - a) * (mb["J"] / (2 * e) + direct_finite)
                Vb = (b - 2 * mu) ** 2 - 4 * mu * mu / (dim - 2)
                T = -pref * mb["I"]
                U = -pref * ma["U"]
                total += 2 * (Vb * Dbox + (b - 2 * mu) * (T - 2 * U))
        return total, max_box_difference

    def finite(self, s, t, E2):
        mu, n = self.mu, self.n
        rows = (s, t, 4 * mu - s - t)
        H = sum(1 / (n - a) for a in rows)
        softsum = (
            sum(
                ((a - 2 * mu) ** 2 - 2 * mu * mu) * self.moments(a, 0.0)["J"]
                for a in rows
            )
            - 2 * mu
        )
        c = np.euler_gamma - math.log(4 * math.pi * E2)
        total = H * (
            sum(
                ((a - 2 * mu) ** 2 - 2 * mu * mu) * self.moments(a, 0.0)["J1"]
                + 2 * mu * mu * self.moments(a, 0.0)["J"]
                for a in rows
            )
            + c * softsum
        )
        for a in rows:
            m = self.moments(a, 0.0)
            q = n - a
            K1 = (
                -2 * (a - 2 * mu) * m["Lll"] / q
                + 4 * n * m["Lmix"] / q
                - 2 * mu * (math.log(mu) - 1) / q
                + (
                    2 * a * a
                    + 4 * (-a * a / 2 + 2 * n * a - n * n) * m["Lh"]
                    - 2 * a * n * (math.log(n) - 1)
                )
                / q**2
            )
            T = -m["I"]
            U = -m["U"]
            total += -K1 + 4 * (a - 2 * mu) * T - 4 * (a * a - 2 * n * mu) * U / q
        for i, a in enumerate(rows):
            for j, b in enumerate(rows):
                if i == j:
                    continue
                m = self.moments(b, 0.0)
                K, _ = self.box(a, b, 0.0)
                V = (b - 2 * mu) ** 2 - 2 * mu * mu
                total += 2 * V * (-m["I"] - a * K) / (n - a)
        return total, H * softsum


@pytest.mark.parametrize(
    "mu,n,ss,tt,nu,E2",
    (
        (1, 9, 2.2, 0.0, 1.0, 1.0),
        (2, 17, 4.3, 0.0, 3.0, 0.5),
        (1, 100, 1.7, 0.5, 2.0, 0.25),
    ),
)
def test_independent_whole_D_raw_graphs_and_finite_implementation(
    mu, n, ss, tt, nu, E2
):
    previous = None
    for nodes in (48, 96):
        whole = IndependentWhole(nodes, mu, n)
        target, pole = whole.finite(ss, tt, E2)
        errors = []
        for ep in (0.001, 0.0001):
            raw, boxdiff = whole.whole_raw(ss, tt, ep, nu)
            value = raw - pole * (E2 / nu) ** ep / ep
            errors.append(abs(value - target))
            assert boxdiff < 0.0002
        assert errors[1] < errors[0] / 5
        if previous is not None:
            assert abs(target - previous) < 0.0002
        previous = target
    mm, nn, aa, ttt, res = (s.Rational(str(x)) for x in (mu, n, ss, tt, E2))
    rows = graphs.channels(aa, ttt, mm)
    replacements = {}
    for a in rows:
        m = whole.moments(float(a), 0.0)
        for expression, value in (
            (forward.J0(a, mm), m["J"]),
            (forward.J1(a, mm), m["J1"]),
            (forward.Lll(a, mm), m["Lll"]),
            (forward.Lmix(mm, nn), m["Lmix"]),
            (forward.Lh(a, nn) + 1, m["Lh"] + 1),
            (-graphs.massive_T(a, mm, nn), m["I"]),
            (-graphs.offshell_U(a, mm, nn), m["U"]),
        ):
            replacements[expression] = s.Float(value, 17)
    for i, a in enumerate(rows):
        for j, b in enumerate(rows):
            if i != j:
                replacements[graphs.finite_box_K(a, b, mm, nn)] = s.Float(
                    whole.box(float(a), float(b), 0.0)[0], 17
                )
    library = forward.finite_bracket(aa, ttt, mm, nn, res).xreplace(replacements)
    assert not library.atoms(s.Integral), library.atoms(s.Integral)
    assert abs(float(library) - target) < 1e-10


@pytest.mark.parametrize(
    "mv,nv,sv,nu,E2",
    ((1, 9, 6, 1, 1), (2, 17, 10, 3, ".5"), (1, 9, 20, 2, ".25"), (1, 9, 50, 1, 1)),
)
def test_independent_full_D_light_cut_and_entire_Coulomb_subtraction(
    mv, nv, sv, nu, E2
):
    with mp.workdps(55):
        mv, nv, sv, nu, E2 = map(mp.mpf, (mv, nv, sv, nu, E2))
        q = sv - 4 * mv
        beta = mp.sqrt(q / sv)
        r = q / (2 * nv + q)
        kh = 4 / (2 * nv + q)
        R = 1 / (nv - sv)
        H = lambda x: R + kh / (1 - r * r * x * x)
        V0 = (sv - 2 * mv) ** 2 - 2 * mv * mv
        Air = H(1) * V0 / (8 * mp.pi * sv * beta)
        target = cuts.hard_light_cut(
            *(s.Rational(str(x)) for x in (sv, mv, nv)), 1, 1, s.Rational(str(E2))
        ).evalf(45)
        errors = []
        for ep in (mp.mpf("1e-5"), mp.mpf("1e-6")):
            V = (sv - 2 * mv) ** 2 - 2 * mv * mv / (1 + ep)
            p = 4 * V / q
            aa = -7 * sv / 4 + 4 * mv + 2 * mv * mv / (sv * (1 + ep))
            bb = -q * q / (4 * sv)
            weight = lambda x, ep=ep: (1 - x * x) ** ep
            norm = mp.quad(weight, [0, 1])
            pint = 1 / (2 * ep) + mp.quad(
                lambda x, weight=weight: weight(x) / (1 + x), [0, 1]
            )
            regular = mp.quad(
                lambda x, weight=weight, p=p, aa=aa, bb=bb: (
                    weight(x)
                    * (
                        -p * kh * r * r / ((1 - r * r) * (1 - r * r * x * x))
                        + H(x) * (aa + bb * x * x)
                    )
                ),
                [0, 1],
            )
            phase = (
                beta
                / (8 * mp.pi)
                * (4 * mp.pi * nu) ** (-ep)
                * q**ep
                * mp.gamma(1 + ep)
                / mp.gamma(2 + 2 * ep)
            )
            literal = phase * (p * H(1) * pint + regular) / (2 * norm)
            library = cuts.light_cut(
                *(s.Rational(str(x)) for x in (sv, mv, nv, ep)),
                1,
                1,
                s.Rational(str(nu)),
            ).evalf(45)
            assert abs(literal - mp.mpf(str(library))) < mp.mpf("1e-38")
            finite = literal - Air * (E2 / nu) ** ep / ep
            errors.append(abs(finite - mp.mpf(str(target))))
        assert errors[1] < errors[0] / 9


@pytest.mark.parametrize("mv,nv", ((1, 9), (2, 17), (1, 100), (3, 13)))
def test_independent_entire_physical_weighted_log(mv, nv):
    with mp.workdps(55):
        mass, heavy = mp.mpf(mv), mp.mpf(nv)
        beta = mp.sqrt(1 - 4 * mass / heavy)

        def integrand(x):
            return (
                -(1 - x * x)
                / 4
                * (
                    mp.log(heavy / 4)
                    + mp.log(abs(x * x - beta * beta))
                    - (1j * mp.pi if x < beta else 0)
                )
            )

        literal = mp.quad(integrand, [0, beta, 1])
        expression = cuts.weighted_bubble(nv, mv)
        expected = mp.mpc(
            str(s.re(expression).evalf(45)), str(s.im(expression).evalf(45))
        )
        assert abs(literal - expected) < mp.mpf("1e-40")


@pytest.mark.parametrize(
    "mv,nv,sv", ((1, 9, 40), (1, 9, 50), (2, 17, 100), (3, 13, 80))
)
def test_independent_entire_HH_interference(mv, nv, sv):
    with mp.workdps(55):
        mass, heavy, energy = map(mp.mpf, (mv, nv, sv))
        A = energy / 2 - heavy
        r = mp.sqrt((energy - 4 * mass) * (energy - 4 * heavy)) / (2 * A)
        alpha = energy / 4 + 2 * mass * heavy / energy
        b = -(energy - 4 * mass) * (energy - 4 * heavy) / (4 * energy)
        integral = mp.quad(
            lambda x: 2 * (alpha + b * x * x) / (A * (1 - r * r * x * x)), [0, 1]
        )
        phase = mp.sqrt(1 - 4 * heavy / energy) / (8 * mp.pi)
        expected = cuts.heavy_cut(sv, mv, nv, 0, 1, 1, 1).evalf(45)
        assert abs(phase * integral / 2 - mp.mpf(str(expected))) < mp.mpf("1e-40")


@pytest.mark.parametrize(
    "args",
    (
        (2, 17, 1),
        (1, 7, 1),
        (1, 9, 0),
        (1, 9, 2),
        (1, 9, 1.0),
        (1.0, 9, 1),
        (True, 9, 1),
        (1, False, 1),
        (1, 9, "1"),
        (1, 9, None),
        (1, 9, s.Symbol("E")),
        (1, 9, s.Rational(1, 8)),
    ),
)
def test_reject_outside_original_bound_domain(args):
    with pytest.raises((TypeError, ValueError)):
        forward.require_bound_domain(*args)


@pytest.mark.parametrize(
    "args",
    (
        ("light", 4, 1, 9),
        ("light", 9, 1, 9),
        ("heavy", 36, 1, 9),
        ("heavy_graviton", 9, 1, 9),
        ("light", 6, 1, 4),
        ("unknown", 6, 1, 9),
        (True, 6, 1, 9),
        ("light", 6.0, 1, 9),
        ("light", 6, False, 9),
        ("light", 6, 1, "9"),
    ),
)
def test_reject_unsupported_cut_domain(args):
    with pytest.raises((TypeError, ValueError)):
        cuts.require_cut_domain(*args)


@pytest.mark.parametrize(
    "species,energy", (("light", 6), ("heavy", 40), ("heavy_graviton", 10))
)
def test_open_cut_domains(species, energy):
    assert cuts.require_cut_domain(species, energy, 1, 9) == (energy, 1, 9)


@pytest.mark.parametrize(
    "name", tuple(forward.data()["whole_exact_majorant_safety_margins"])
)
def test_each_analytic_majorant_has_exact_nonnegative_margin(name):
    value = s.sympify(forward.data()["whole_exact_majorant_safety_margins"][name])
    assert value >= 0 and not value.atoms(s.Float)


@pytest.mark.parametrize("aa,bb", ((2, 0), (2, 1)))
def test_actual_subtracted_box_negative_finite_K_coefficient(aa, bb):
    ep = s.Rational(1, 16)
    expression = graphs.subtracted_box(aa, bb, 1, 9, ep, 1)
    kval = graphs.finite_box_K(aa, bb, 1, 9, ep)
    atom = s.Symbol("independent_box_K")
    replaced = expression.xreplace({kval: atom})
    coefficient = s.diff(replaced, atom)
    expected = -s.gamma(1 - ep) * (4 * s.pi) ** (-ep) * (1 - ep) * aa / (9 - aa)
    assert s.simplify(coefficient - expected) == 0


@pytest.mark.parametrize("dimension", (4, 5, 6))
def test_two_metric_cubic_contact_whole_projector_trace(dimension):
    eta = s.diag(1, *([-1] * (dimension - 1)))
    contraction = s.trace(eta * eta * eta * eta) - s.trace(eta * eta) ** 2 / s.Integer(
        dimension - 2
    )
    assert contraction == -s.Rational(2 * dimension, dimension - 2)


@pytest.mark.parametrize("curvature", (1, -1))
def test_RH_anchor_is_retained_not_folded_into_known_bound(curvature):
    whole = cuts.mixing_transition(1, 9, 1, 1, curvature)
    known = cuts.mixing_transition(1, 9, 1, 1, 0)
    assert s.simplify(whole - known - 11 * curvature) == 0


def test_all_twelve_labelled_box_multiplicities():
    rows = graphs.data()["whole_labelled_mixed_box_counts"]
    assert len(rows) == 6 and sum(row[2] for row in rows) == 12
    assert {(a, b) for a, b, _ in rows} == {
        (a, b) for a in "stu" for b in "stu" if a != b
    }


def test_cached_parent_checks_not_mutated():
    parent = source.previous.data()
    before = dict(parent["checks"])
    source.data()
    assert parent["checks"] == before and len(before) == 114
    assert source.data()["checks"] is not parent["checks"]


@pytest.mark.parametrize("coefficient", (2001, 2002))
def test_exact_original_known_bound_not_numeric_quadrature(coefficient):
    value = coefficient * source.CUBIC**2 / (source.KAPPA * source.HEAVY_MASS2)
    assert value < s.Rational(1, 10**1001)
    assert not value.atoms(s.Float)
