"""Independent matrix, dimensional spectral, finite-part and ownership checks."""

import copy

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_insertion_ms import audit as previous
from p8_vacuum_fermion_ms_slopes import tensors as t
from p8_vacuum_fermion_vacuum import audit, calibration, forest, scalar


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_exact_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_invalid_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def gamma_matrices():
    sigma = (
        mp.matrix([[0, 1], [1, 0]]),
        mp.matrix([[0, -1j], [1j, 0]]),
        mp.matrix([[1, 0], [0, -1]]),
    )
    g = []
    for x in sigma:
        g.append(
            mp.matrix(
                [
                    [0, 0, -1j * x[0, 0], -1j * x[0, 1]],
                    [0, 0, -1j * x[1, 0], -1j * x[1, 1]],
                    [1j * x[0, 0], 1j * x[0, 1], 0, 0],
                    [1j * x[1, 0], 1j * x[1, 1], 0, 0],
                ]
            )
        )
    g.append(mp.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]))
    g.append(g[0] * g[1] * g[2] * g[3])
    return g


@pytest.mark.parametrize("dim", (3, 4, 5))
@pytest.mark.parametrize("config", (0, 1))
@pytest.mark.parametrize("sector", ("scalar", "gauge"))
def test_explicit_matrix_vacuum_numerator_and_dimension(dim, config, sector):
    with mp.workdps(50):
        gam = gamma_matrices()[:dim]
        kv = [sp.Rational(j + 1 + config, 5) for j in range(dim)]
        lv = [sp.Rational((-1) ** j * (j + 2), 7 + config) for j in range(dim)]
        km = [mp.mpf(int(v.p)) / int(v.q) for v in kv]
        lm = [mp.mpf(int(v.p)) / int(v.q) for v in lv]

        def S(v):
            return (
                mp.eye(4) - 1j * sum((g * x for g, x in zip(gam, v)), mp.zeros(4))
            ) / (1 + sum(x * x for x in v))

        sk, sl = S(km), S(lm)
        trace = lambda a: sum(a[j, j] for j in range(4))
        if sector == "scalar":
            matrix = trace(sk * sl) / 2
            expr = t.trace(t.Sk, t.Sl) / 2
        else:
            matrix = -sum(trace(g * sk * g * sl) for g in gam) / 2
            expr = -t.trace(t.mu, t.Sk, t.mu, t.Sl) / 2
        x, y, z = (
            sum(v * v for v in kv),
            sum(v * v for v in lv),
            sum(a * b for a, b in zip(kv, lv)),
        )
        exact = expr.subs({t.x: x, t.y: y, t.z: z, t.A: 1 + x, t.B: 1 + y, t.d: dim})
        assert abs(matrix - mp.mpf(str(sp.N(exact, 50)))) < mp.mpf("1e-42")


@pytest.mark.parametrize("dim", (3, 4, 5))
@pytest.mark.parametrize("config", (0, 1))
def test_vector_Ward_resolvent_by_explicit_noncommuting_matrices(dim, config):
    with mp.workdps(50):
        gam = gamma_matrices()[:dim]
        k = [mp.mpf(j + 1) / 5 for j in range(dim)]
        q = [mp.mpf((-1) ** j * (j + 1 + config)) / 7 for j in range(dim)]
        D0 = mp.eye(4) + 1j * sum((g * x for g, x in zip(gam, k)), mp.zeros(4))
        vertex = 1j * sum((g * x for g, x in zip(gam, q)), mp.zeros(4))
        D1 = D0 + vertex
        residual = D0 ** (-1) * vertex * D1 ** (-1) - (D0 ** (-1) - D1 ** (-1))
        assert max(abs(x) for x in residual) < mp.mpf("1e-44")


def vacuum_normal(e, sector, L=0):
    A = mp.exp(mp.euler * e) * mp.gamma(1 + e)
    if sector == "scalar":
        V = 4 / ((1 - e) * (2 * e - 1)) + 1 / (1 - e) ** 2
        c = 6 / (1 - e)
    else:
        V = -4 / ((1 - e) * (2 * e - 1)) + 2 / (1 - e)
        c = -12 / (1 - e)
    return A * A * V * mp.exp(-2 * e * L) + A * c * mp.exp(-e * L)


@pytest.mark.parametrize("sector", ("scalar", "gauge"))
@pytest.mark.parametrize("L", ("0", "2", "-1"))
def test_independent_complex_finite_part_at_fixed_scale(sector, L):
    with mp.workdps(70):
        L = mp.mpf(L)
        vals = []
        for j in range(32):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 32)
            # The Fourier average selects the constant Laurent coefficient
            # without importing any reported pole or finite coefficients.
            vals.append(vacuum_normal(e, sector, L) / e**2)
        expected = (
            -3 * L * L + 14 * L - 19 if sector == "scalar" else 6 * L * L - 16 * L + 18
        )
        assert abs(sum(vals) / 32 - expected) < mp.mpf("1e-48")


@pytest.mark.parametrize("sector", ("scalar", "gauge"))
def test_independent_finite_vacuum_fixed_mu_mass_derivative(sector):
    with mp.workdps(55):
        mu = mp.mpf("1.3")

        def V(m):
            L = mp.log(m * m / (mu * mu))
            polynomial = (
                -3 * L * L + 14 * L - 19
                if sector == "scalar"
                else 6 * L * L - 16 * L + 18
            )
            return m**4 * polynomial

        expected = -56 if sector == "scalar" else 40
        assert abs(mp.diff(V, mu, 2) / mu**2 - expected) < mp.mpf("1e-42")


def H(e):
    return (
        mp.exp(2 * mp.euler * e)
        * 4 ** (-e)
        * mp.sqrt(mp.pi)
        / (2 * mp.gamma(mp.mpf("1.5") - e))
    )


def F0(e):
    return (
        mp.exp(2 * mp.euler * e)
        * 4 ** (-e)
        * mp.sqrt(mp.pi)
        / 2
        * (mp.mpf("1.5") - e)
        / (1 - e)
        * mp.gamma(e)
        * mp.gamma(2 * e)
        / mp.gamma(mp.mpf("2.5") + e)
    )


@pytest.mark.parametrize("e", (".23", "0.6", "1.1", "1.3"))
@pytest.mark.parametrize("s", ("-3", "-10"))
def test_dimensional_inner_spectral_identity_in_vacuum_convergence_region(e, s):
    with mp.workdps(60):
        e, s = map(mp.mpf, (e, s))
        m = mp.mpf(2)
        T = 4 * m * m
        pref = mp.exp(mp.euler * e) * m ** (2 * e) * mp.gamma(e)

        def B(s):
            return pref * mp.quad(
                lambda x: (m * m - x * (1 - x) * s) ** (-e), [0, mp.mpf(".5"), 1]
            )

        B1 = B(1)
        first_derivative = -B1 + (T - 1) * pref * e * mp.quad(
            lambda x: x * (1 - x) * (m * m - x * (1 - x)) ** (-e - 1),
            [0, mp.mpf(".5"), 1],
        )
        direct = (T - s) * B(s) - (T - 1) * B1 - (s - 1) * first_derivative
        rho_pref = (
            mp.exp(mp.euler * e)
            * m ** (2 * e)
            * 4**e
            * mp.sqrt(mp.pi)
            / (2 * mp.gamma(mp.mpf("1.5") - e))
        )
        compact = mp.quad(
            lambda z: (
                z**e
                * (1 - z) ** (mp.mpf("1.5") - e)
                / ((1 - z / T) ** 2 * (1 - s * z / T))
            ),
            [0, mp.mpf(".5"), 1],
        )
        spectral = -((s - 1) ** 2) * rho_pref * T ** (-1 - e) * compact
        assert abs(direct / spectral - 1) < mp.mpf("1e-40")


@pytest.mark.parametrize("e", ("1.1", "1.25", "1.4"))
@pytest.mark.parametrize("term", (0, 1, 2))
def test_all_three_beta_anchors_by_convergent_integrals(e, term):
    with mp.workdps(70):
        e = mp.mpf(e)
        integral = mp.quad(
            lambda z: z ** (2 * e - 3 + term) * (1 - z) ** (mp.mpf("1.5") - e),
            [0, mp.mpf(".5"), 1],
        )
        direct = mp.gamma(e - 1) * H(e) * integral
        expected = (
            -F0(e)
            * (e + mp.mpf(".5"))
            * (e + mp.mpf("1.5"))
            / ((2 * e - 2) * (2 * e - 1)),
            F0(e) * (mp.mpf("1.5") + e) / (1 - 2 * e),
            -F0(e),
        )[term]
        assert abs(direct / expected - 1) < mp.mpf("1e-13")


@pytest.mark.parametrize("T", ("16", "100"))
@pytest.mark.parametrize("e", ("1.1", "1.3"))
def test_actual_scalar_vacuum_spectral_and_compact_integrals(T, e):
    with mp.workdps(70):
        T, e = map(mp.mpf, (T, e))
        m = mp.sqrt(T) / 2
        pref = (
            mp.exp(2 * mp.euler * e)
            * m ** (4 * e)
            * 4**e
            * mp.sqrt(mp.pi)
            * mp.gamma(e - 1)
            / (2 * mp.gamma(mp.mpf("1.5") - e))
        )
        spectral = pref * mp.quad(
            lambda v: v ** (2 - 2 * e) * (1 - T / v) ** (mp.mpf("1.5") - e) / (v - 1),
            [T, 2 * T, mp.inf],
        )
        compact = (
            T
            * T
            * H(e)
            * mp.gamma(e - 1)
            * mp.quad(
                lambda z: (
                    z ** (2 * e - 3) * (1 - z) ** (mp.mpf("1.5") - e) / (1 - z / T)
                ),
                [0, mp.mpf(".5"), 1],
            )
        )
        assert abs(spectral / compact - 1) < mp.mpf("1e-12")


def vacuum_remainder_D(T, e):
    return mp.quad(
        lambda z: z ** (2 * e) * (1 - z) ** (mp.mpf("1.5") - e) / (T**3 * (1 - z / T)),
        [0, mp.mpf(".5"), 1],
    )


def vacuum_remainder_finite(T):
    return mp.quad(
        lambda z: (
            (1 - z) ** mp.mpf("1.5")
            / (T * (1 - z / T))
            * (4 * mp.log(2) - 3 - 2 * mp.log(z) + mp.log(1 - z))
        ),
        [0, mp.mpf(".5"), 1],
    )


@pytest.mark.parametrize("T", ("16", "100", "1e6"))
def test_whole_vacuum_remainder_finite_pole_product(T):
    with mp.workdps(55):
        T = mp.mpf(T)
        D0 = vacuum_remainder_D(T, 0)
        vals = []
        for j in range(16):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 16)
            vals.append(
                T * T * (mp.gamma(e - 1) * H(e) * vacuum_remainder_D(T, e) + D0 / e)
            )
        assert abs(sum(vals) / 16 - vacuum_remainder_finite(T)) < mp.mpf("1e-35")
        assert D0 > 0


@pytest.mark.parametrize("T", ("16", "100", "1e6", "1e40", "1e400"))
def test_complete_scalar_vacuum_remainder_majorant(T):
    with mp.workdps(70):
        T = mp.mpf(T)
        assert abs(vacuum_remainder_finite(T)) < 12 / T


@pytest.mark.parametrize("k", ("0", ".0001", ".01", ".0625"))
def test_three_term_scalar_remainder_identity_and_majorant(k):
    with mp.workdps(60):
        k = mp.mpf(k)
        h3 = k**3 / (1 - k)
        assert abs(1 / (1 - k) - 1 - k - k * k - h3) < mp.mpf("1e-55")
        assert 0 <= h3 <= 2 * k**3


@pytest.mark.parametrize(
    "m,Y,Q",
    ((2, 0, 144), (3, sp.Rational(1, 7), 144), (10**5, sp.Rational(1, 10**9), 100)),
)
def test_exact_scalar_vacuum_enclosure_keeps_mass_correction(m, Y, Q):
    d = scalar.enclosure(m, Y, Q)
    twice = scalar.enclosure(m, 2 * Y, Q)
    assert (
        twice["complete_scalar_vacuum_absolute_upper"]
        == 2 * d["complete_scalar_vacuum_absolute_upper"]
    )
    assert (
        d["complete_scalar_vacuum_absolute_upper"]
        == d["massless_scalar_vacuum_absolute_upper"]
        + d["actual_scalar_vacuum_mass_correction_absolute_upper"]
    )
    if Y > 0:
        assert d["finite_scalar_vacuum_remainder_absolute_upper"] > 0


@pytest.mark.parametrize(
    "m,a,Q",
    ((2, 0, 144), (3, sp.Rational(1, 7), 144), (10**5, sp.Rational(1, 10**9), 100)),
)
def test_gauge_vacuum_all_flavor_normalization(m, a, Q):
    result = calibration.gauge_enclosure(m, a, Q)
    assert result == 7 * (6 * a * sp.Rational(4, 3) * 18 * m**4 / Q**2)
    assert calibration.gauge_enclosure(m, 2 * a, Q) == 2 * result


def test_three_overlapping_multigraph_subcycles():
    edges = ("fermion_1", "fermion_2", "boson")
    actual = {frozenset(v) for v in forest.data()["proper_subcycles"]}
    expected = {
        frozenset((edges[i], edges[j])) for i in range(3) for j in range(i + 1, 3)
    }
    assert actual == expected
    assert all(a & b for a in actual for b in actual if a != b)


def test_flavor_ledger_counts_active_and_inert_separately():
    rows = forest.data()["flavor_ledger"]
    assert len(rows) == 42
    assert len({(r["flavor"], r["color"]) for r in rows}) == 42
    assert sum(r["active_Yukawa"] for r in rows) == 6
    assert sum(not r["active_Yukawa"] for r in rows) == 36
    assert {r["flavor"] for r in rows} == set(range(14))


def test_actual_absolute_and_relative_reference_scales_both_retained():
    d = calibration.data()
    assert all(v is True for v in d["bounds"].values())
    assert 10**590 < d["both_vacuum_primitive_absolute_upper"] < 10**595
    assert d["relative_to_one_loop_vacuum_upper"] < sp.Rational(1, 10**203)
    assert (
        d["both_vacuum_primitive_absolute_upper"]
        > d["all_flavor_gauge_vacuum_absolute_upper"]
        > 0
    )
    assert (
        d["scalar_vacuum_enclosure"][
            "actual_scalar_vacuum_mass_correction_absolute_upper"
        ]
        > 0
    )


def test_exact_two_vacuum_row_advance_without_full_matching_claim():
    old = copy.deepcopy(previous.frontier())
    new = audit.frontier()
    assert [a["id"] for a, b in zip(old, new) if a != b] == list(audit.TARGETS)
    assert previous.frontier() == old
    assert all(r["status"].startswith("BOUNDED") for r in new)
    assert all(r["status"] != "COMPLETE" for r in new)
    assert sum(a == b for a, b in zip(old, new)) == 7


def test_controls_keep_complete_vacuum_and_other_ledgers_open():
    controls = audit.controls()
    assert len(controls) == 9
    assert controls["vacuum_energy_reference_not_yet_full_model_sum"] is True
    assert controls["first_source_square_and_other_vacuum_terms_open"] is True
    assert controls["other_matching_canonical_and_truncation_terms_open"] is True
    assert controls["original_P8_not_closed"] is True
