"""Independent complex, integral, fixed-mu and forest checks for the mixed row."""

import copy

import mpmath as mp
import pytest
from p8_vacuum_fermion_mixed_quartic import (
    audit,
    bounds,
    catalog,
    joint,
)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_exact_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_inputs_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("x", ("0", ".001", ".3", "1", "100", "1000000"))
@pytest.mark.parametrize("theta", ("0", ".4", "1.5", "2.7", "4.3", "6"))
def test_complex_pair_denominator_gap(x, theta):
    with mp.workdps(40):
        x, theta = mp.mpf(x), mp.mpf(theta)
        z = 2 + mp.exp(1j * theta)
        for fraction in (0, mp.mpf(".25"), 1):
            P = (1 + x - z / 4) ** 2 + z * x * fraction
            assert abs(P) >= (1 + x) ** 2 / 16
            diff = 1 / P - 1 / (1 + x) ** 2
            assert abs(diff) <= 33 / (1 + x) ** 3


@pytest.mark.parametrize("a,b", joint.data()["exponents"])
def test_sunset_beta_integrals_independently(a, b):
    with mp.workdps(50):
        a, b = mp.mpf(str(a.p)) / int(a.q), mp.mpf(str(b.p)) / int(b.q)
        g = mp.mpf("1.75")
        # Three separately integrated Schwinger/Beta factors, not Gamma-form substitution.
        first = mp.quad(lambda v: v ** (g - 1) / (1 + v) ** 2, [0, 1, mp.inf])
        second = mp.quad(lambda t: t ** (a + b + g - 5) * mp.exp(-t), [0, 1, mp.inf])
        third = mp.quad(
            lambda x: x ** (a + g - 3) * (1 - x) ** (b + g - 3), [0, 0.5, 1]
        )
        numerical = first * second * third / (mp.gamma(a) * mp.gamma(b) * mp.gamma(g))
        exact = (
            mp.gamma(2 - g)
            * mp.gamma(a + b + g - 4)
            * mp.gamma(a + g - 2)
            * mp.gamma(b + g - 2)
            / (mp.gamma(a) * mp.gamma(b) * mp.gamma(a + b + 2 * g - 4))
        )
        assert abs(numerical / exact - 1) < mp.mpf("1e-12")
        assert 0 < exact < 80000


@pytest.mark.parametrize("m,x", ((2, 0), (2, 1), (2, 100), (5, 0), (5, 10), (20, 1000)))
def test_zero_soft_spectral_kernel_against_fixed_mu_mass_derivative(m, x):
    with mp.workdps(45):
        m, x = mp.mpf(m), mp.mpf(x)
        mu = m

        def f(mass):
            J = mp.quad(lambda a: mp.log((mass**2 + a * (1 - a) * x) / mu**2), [0, 1])
            return -(4 * mass**2 + x) * J + 2 * mass**2 * (1 - mp.log(mass**2 / mu**2))

        direct = mp.diff(f, m, 2)
        T = 4 * m * m
        # Units CY. z=1-u^2 removes the endpoint square-root singularity.
        spectral = -32 + 24 * (x / T) * mp.quad(
            lambda u: (1 - 2 * u * u) / (1 + (x / T) * (1 - u * u)), [0, 1]
        )
        assert abs(direct - spectral) < mp.mpf("1e-35")
        assert abs(direct) <= 68 + 24 * mp.log(1 + x)


@pytest.mark.parametrize("e", (".15", ".25", ".35"))
def test_raw_box_spectral_zero_reference_at_nonzero_regulator(e):
    with mp.workdps(65):
        e = mp.mpf(e)
        m = mp.mpf(3)
        T = 4 * m * m
        p = mp.mpf("1.5") - e
        pref = (
            mp.exp(mp.euler * e)
            * m ** (2 * e)
            * 4**e
            * mp.sqrt(mp.pi)
            / (2 * mp.gamma(p))
            * 8
            * p
        )
        # Direct convergent moment over the threshold variable.
        integral = mp.quad(
            lambda z: (
                z ** (e - 1) * (1 - z) ** (-mp.mpf(".5") - e) * ((2 - 2 * e) * z - 1)
            ),
            [0, 0.5, 1],
        )
        raw = -pref * T ** (-e) * integral
        expected = mp.exp(mp.euler * e) * mp.gamma(e) * (12 - 32 * e + 16 * e * e)
        assert abs(raw / expected - 1) < mp.mpf("1e-8")


def leading_parts(e, ell):
    soft = mp.exp(ell * e) * (
        mp.exp(2 * mp.euler * e) * (12 - 32 * e + 16 * e * e) * mp.gamma(e) ** 2
        - 12 * mp.exp(mp.euler * e) * mp.gamma(e) / e
    )
    hard = (
        -8
        * mp.sqrt(mp.pi)
        * mp.exp(2 * mp.euler * e)
        * 4 ** (-e)
        * (mp.mpf("1.5") - e)
        * (1 - 4 * e)
        * mp.gamma(e)
        * mp.gamma(2 * e)
        / ((1 - e) * (1 + 2 * e) * mp.gamma(mp.mpf(".5") + e))
    )
    return soft, hard


@pytest.mark.parametrize("ell", ("0", "3.2", "900"))
def test_independent_complex_laurent_finite_extraction(ell):
    with mp.workdps(65):
        ell = mp.mpf(ell)
        values = []
        for j in range(32):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 32)
            a, b = leading_parts(e, ell)
            values.append(a + b + 6 / e**2 - 2 / e)
        finite = sum(values) / 32
        assert abs(finite + 46 + 32 * ell) < mp.mpf("1e-45")


@pytest.mark.parametrize("e", (".1", ".2", ".3"))
def test_beta_regions_match_full_regulated_leading_reference(e):
    with mp.workdps(45):
        e = mp.mpf(e)
        m = mp.mpf(5)
        T = 4 * m * m
        ell = mp.log(m * m)
        H = (
            mp.exp(2 * mp.euler * e)
            * 4 ** (-e)
            * mp.sqrt(mp.pi)
            / (2 * mp.gamma(mp.mpf("1.5") - e))
        )
        I = lambda n: (
            mp.beta(n * e, mp.mpf(".5") - e)
            * ((2 - 2 * e) * n * e / ((n - 1) * e + mp.mpf(".5")) - 1)
        )
        raw = (
            -8 * (mp.mpf("1.5") - e) * mp.gamma(e) * H * (T**e * I(1) - I(2) / (1 - e))
        )
        ct = -12 * mp.exp((mp.euler + ell) * e) * mp.gamma(e) / e
        a, b = leading_parts(e, ell)
        assert abs(raw + ct - a - b) < mp.mpf("1e-38")


def correction_4d(r):
    def integrand(u):
        z = 1 - u * u
        if z == 0:
            return mp.mpf(0)  # An endpoint value on a measure-zero set only.
        k = r * z
        h = k * (2 - k) / (1 - k) ** 2
        return -24 * (2 * z - 1) / z * (h * (-mp.log(k) - 1) + k / (1 - k) ** 2)

    return mp.quad(integrand, [0, 0.5, 1])


@pytest.mark.parametrize("r", (".0625", ".01", ".000001", "1e-100"))
def test_exact_finite_mass_ratio_correction_is_bounded_and_not_zero(r):
    with mp.workdps(55):
        r = mp.mpf(r)
        v = correction_4d(r)
        assert v != 0
        assert abs(v) < r * (300 - 100 * mp.log(r))


@pytest.mark.parametrize("v", ("16", "100", "1000000"))
def test_triangle_integral_and_leading_difference(v):
    with mp.workdps(45):
        v = mp.mpf(v)
        direct = mp.quad(lambda x: x / ((x + v) * (x + 1) ** 2), [0, 1, v, mp.inf])
        exact = (v * mp.log(v) - v + 1) / (v - 1) ** 2
        assert abs(direct - exact) < mp.mpf("1e-38")
        assert exact > (mp.log(v) - 1) / v


@pytest.mark.parametrize("M", (10000, 10**8, 10**40))
def test_heavy_triangle_radial_majorants(M):
    with mp.workdps(45):
        M = mp.mpf(M)
        ell = mp.log(2 * M)
        # Logarithmic coordinate integrates both sides of the moving region split.
        f = lambda u: mp.exp(2 * u) / ((1 + mp.exp(u)) ** 2 * (M + mp.exp(u)))
        J = mp.quad(f, [-mp.inf, 0, mp.log(M), mp.inf])
        Jlog = mp.quad(
            lambda u: f(u) * mp.log(1 + mp.exp(u)), [-mp.inf, 0, mp.log(M), mp.inf]
        )
        assert 0 < J < (ell + 1) / M
        assert 0 < Jlog < (ell * ell / 2 + ell + 1) / M


@pytest.mark.parametrize("n", (1, 2, 3, 7, 16, 1000, 10**400))
def test_dyadic_bound_is_least(n):
    k = bounds.dyadic(n)
    assert n <= 2**k
    assert k == 0 or 2 ** (k - 1) < n


@pytest.mark.parametrize("word", catalog.data()["cyclic_box_words"])
def test_word_routes_and_actual_complex_fermion_norm(word):
    # Check all external routes using 4x4 Euclidean Clifford matrices.
    with mp.workdps(35):
        sigma = (
            mp.matrix([[0, 1], [1, 0]]),
            mp.matrix([[0, -1j], [1j, 0]]),
            mp.matrix([[1, 0], [0, -1]]),
        )
        gam = []
        for s in sigma:
            gam.append(
                mp.matrix(
                    [
                        [0, 0, -1j * s[0, 0], -1j * s[0, 1]],
                        [0, 0, -1j * s[1, 0], -1j * s[1, 1]],
                        [1j * s[0, 0], 1j * s[0, 1], 0, 0],
                        [1j * s[1, 0], 1j * s[1, 1], 0, 0],
                    ]
                )
            )
        gam.append(mp.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]))
        m = mp.mpf(720)
        k = mp.matrix([10, 20, 30, 40])
        q = mp.matrix([50, -20, 10, 0])
        l = k + q
        S = lambda v: m * m + sum(t * t for t in v)
        R = min(mp.sqrt(S(k)), mp.sqrt(S(l))) / 360
        z = 2 + mp.exp(0.7j)
        E = mp.sqrt(z) / 2
        P = mp.sqrt(z / 4 - 1)
        p1 = mp.matrix([-1j * E, 0, 0, -P])
        p2 = mp.matrix([-1j * E, 0, 0, P])
        K = p1 + p2
        vectors = {"A": q - K / 2, "B": -q - K / 2, "1": p1, "2": p2}
        soft = {"A": -K / 2, "B": -K / 2, "1": p1, "2": p2}
        assert mp.norm(sum(vectors.values(), mp.zeros(4, 1))) < mp.mpf("1e-30")
        base = k.copy()
        shift = mp.zeros(4, 1)
        for label in word:
            if label == "A":
                base += q
            if label == "B":
                base -= q
            shift += soft[label]
            assert sum(abs(a) for a in shift) < 18
            momentum = base + R * mp.exp(0.3j) * shift
            D = m * mp.eye(4) + 1j * sum(
                (gam[i] * momentum[i] for i in range(4)), mp.zeros(4)
            )
            prop = D ** (-1)
            herm = 4 / S(base) * mp.eye(4) - prop.H * prop
            vals = mp.eighe(herm, eigvals_only=True)
            assert min(vals) > 0


@pytest.mark.parametrize("i", range(9))
def test_each_frontier_row_cannot_disappear(i):
    rows = audit.frontier()
    del rows[i]
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


def test_scope_and_zero_coupling():
    r = bounds.enclosure(720, 0, 1, 1, 10000, 144)
    assert r["complete_paired_mixed_b2_upper"] == 0
    assert len([r for r in audit.frontier() if r["status"] == "UNEVALUATED"]) == 4
    assert "other" in bounds.data()["enclosure"]["scope"]
    altered = copy.deepcopy(audit.frontier())
    altered[0]["status"] = "CLOSED"
    with pytest.raises(ValueError):
        audit.validate_frontier(altered)


@pytest.mark.parametrize("r", (".0625", ".01"))
def test_full_dimensional_mass_ratio_difference_finite_extraction(r):
    with mp.workdps(55):
        r = mp.mpf(r)
        m = mp.sqrt(1 / (4 * r))

        def delta_dimensional(e):
            p = mp.mpf("1.5") - e
            pref = (
                8
                * p
                * r
                * mp.exp(2 * mp.euler * e)
                * m ** (2 * e)
                * mp.sqrt(mp.pi)
                * mp.gamma(1 + e)
                / (mp.gamma(p) * (1 - e))
            )

            def f(u):
                z = 1 - u * u
                if z == 0 or u == 0:
                    return mp.mpf(0)
                k = r * z
                # Exact regulated triangle difference, rewritten without cancellation.
                numerator = (1 - k) + (2 - k) * mp.expm1(e * mp.log(k)) / e
                return (
                    u ** (-2 * e)
                    * z**e
                    * ((2 - 2 * e) * z - 1)
                    * numerator
                    / (1 - k) ** 2
                )

            return pref * mp.quad(f, [0, 0.5, 1])

        average = (
            sum(
                delta_dimensional(mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 8))
                for j in range(8)
            )
            / 8
        )
        assert abs(average - correction_4d(r)) < mp.mpf("1e-27")
