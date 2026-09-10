"""Independent regulated integrals, finite-part extraction and ownership checks."""

import copy

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_insertion_ms import audit, bubble, calibration, tadpole
from p8_vacuum_fermion_ms_mass import audit as previous


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


def H(e):
    return (
        mp.exp(2 * mp.euler * e)
        * 4 ** (-e)
        * mp.sqrt(mp.pi)
        / (2 * mp.gamma(mp.mpf("1.5") - e))
    )


def leading_tad(e):
    return (
        mp.exp(2 * mp.euler * e)
        * 4 ** (-e)
        * mp.sqrt(mp.pi)
        / 2
        * (mp.mpf("1.5") - e)
        * mp.gamma(e - 1)
        * mp.gamma(2 * e - 1)
        / mp.gamma(mp.mpf("1.5") + e)
    )


def subleading_tad(e):
    return (
        mp.exp(2 * mp.euler * e)
        * 4 ** (-e)
        * mp.sqrt(mp.pi)
        * (mp.mpf("1.5") - e)
        * mp.gamma(e - 1)
        * mp.gamma(2 * e)
        / mp.gamma(mp.mpf("2.5") + e)
    )


@pytest.mark.parametrize("which", ("leading", "subleading"))
def test_independent_complex_Laurent_anchors(which):
    with mp.workdps(70):
        if which == "leading":
            f, p2, p1, finite = (
                leading_tad,
                mp.mpf(".75"),
                mp.mpf(".25"),
                mp.mpf("3.25") + mp.pi**2 / 8,
            )
        else:
            f, p2, p1, finite = (
                subleading_tad,
                -1,
                mp.mpf(7) / 3,
                -mp.mpf(47) / 9 - mp.pi**2 / 6,
            )
        vals = []
        for j in range(32):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 32)
            vals.append(f(e) - p2 / e**2 - p1 / e)
        assert abs(sum(vals) / 32 - finite) < mp.mpf("1e-48")


@pytest.mark.parametrize("e", (".6", ".75", ".85"))
@pytest.mark.parametrize("which", ("leading", "subleading"))
def test_independent_convergent_beta_anchors(e, which):
    with mp.workdps(65):
        e = mp.mpf(e)
        power = 2 * e - 2 if which == "leading" else 2 * e - 1
        compact = mp.quad(
            lambda z: z**power * (1 - z) ** (mp.mpf("1.5") - e), [0, mp.mpf(".5"), 1]
        )
        value = H(e) * mp.gamma(e - 1) * compact
        expected = leading_tad(e) if which == "leading" else subleading_tad(e) / 2
        assert abs(value / expected - 1) < mp.mpf("1e-12")


@pytest.mark.parametrize("T", ("16", "100"))
@pytest.mark.parametrize("e", (".6", ".8"))
def test_complete_regulated_tadpole_spectral_change_of_variables(T, e):
    with mp.workdps(65):
        T, e = map(mp.mpf, (T, e))
        m = mp.sqrt(T) / 2
        spectral_pref = (
            mp.exp(2 * mp.euler * e)
            * m ** (4 * e)
            * 4**e
            * mp.sqrt(mp.pi)
            * mp.gamma(e - 1)
            / (2 * mp.gamma(mp.mpf("1.5") - e))
        )
        spectral = (
            mp.quad(
                lambda v: (
                    v ** (2 - 2 * e) * (1 - T / v) ** (mp.mpf("1.5") - e) / (v - 1) ** 2
                ),
                [T, 2 * T, mp.inf],
            )
            * spectral_pref
        )
        compact = (
            T
            * H(e)
            * mp.gamma(e - 1)
            * mp.quad(
                lambda z: (
                    z ** (2 * e - 2) * (1 - z) ** (mp.mpf("1.5") - e) / (1 - z / T) ** 2
                ),
                [0, mp.mpf(".5"), 1],
            )
        )
        assert abs(spectral / compact - 1) < mp.mpf("1e-11")


def h2(k):
    return k * k * (3 - 2 * k) / (1 - k) ** 2


def tad_compact_remainder(T, e):
    return mp.quad(
        lambda z: z ** (2 * e - 2) * (1 - z) ** (mp.mpf("1.5") - e) * h2(z / T),
        [0, mp.mpf(".5"), 1],
    )


def tad_finite_remainder(T):
    return T * mp.quad(
        lambda z: (
            (1 - z) ** mp.mpf("1.5")
            / z**2
            * h2(z / T)
            * (4 * mp.log(2) - 3 - 2 * mp.log(z) + mp.log(1 - z))
        ),
        [0, mp.mpf(".5"), 1],
    )


@pytest.mark.parametrize("T", ("16", "100", "1e6"))
def test_complete_tadpole_remainder_finite_pole_product(T):
    with mp.workdps(55):
        T = mp.mpf(T)
        D0 = tad_compact_remainder(T, 0)
        vals = []
        for j in range(16):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 16)
            vals.append(
                T * (mp.gamma(e - 1) * H(e) * tad_compact_remainder(T, e) + D0 / e)
            )
        expected = tad_finite_remainder(T)
        assert abs(sum(vals) / 16 - expected) < mp.mpf("1e-36")
        assert D0 > 0


@pytest.mark.parametrize("T", ("16", "100", "1e6", "1e40", "1e400"))
def test_finite_tadpole_remainder_full_majorant(T):
    with mp.workdps(70):
        T = mp.mpf(T)
        assert abs(tad_finite_remainder(T)) < 24 / T


@pytest.mark.parametrize("k", ("0", ".0001", ".01", ".0625"))
def test_tadpole_two_power_remainder_pointwise(k):
    with mp.workdps(60):
        k = mp.mpf(k)
        assert abs((1 - k) ** (-2) - 1 - 2 * k - h2(k)) < mp.mpf("1e-55")
        assert 0 <= h2(k) <= 4 * k * k


def h1(k):
    return k * (2 - k) / (1 - k) ** 2


def bubble_D(T, M, e):
    def f(z):
        k = z / T
        j = M * z / T
        A = 1 / (1 - k) ** 2
        # Stable exact subtraction of the leading F0 at arbitrary e.
        delta = h1(k) / (1 - e) - A * j * mp.expm1(-e * mp.log(j)) / ((1 - e) * (1 - j))
        return z ** (2 * e - 1) * (1 - z) ** (mp.mpf("1.5") - e) * delta

    return mp.quad(f, [0, mp.mpf(".5"), 1])


@pytest.mark.parametrize("T,M", (("16", "1"), ("100", "4"), ("1e8", "10")))
def test_full_general_mass_bubble_finite_part_includes_pole_product(T, M):
    with mp.workdps(55):
        T, M = map(mp.mpf, (T, M))
        D0 = bubble_D(T, M, 0)
        vals = []
        for j in range(16):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 16)
            vals.append(mp.gamma(e) * H(e) * bubble_D(T, M, e) - D0 / e)
        assert abs(sum(vals) / 16 - bubble_finite_difference(T, M)) < mp.mpf("1e-35")
        assert D0 > 0
        assert abs(D0 - bubble_D(T, 1, 0)) < mp.mpf("1e-48")


@pytest.mark.parametrize(
    "T,M",
    (("16", "1"), ("100", "1"), ("100", "6.25"), ("1e8", "10"), ("1e400", "1e197")),
)
def test_full_two_ratio_bubble_finite_remainder_majorant(T, M):
    with mp.workdps(70):
        T, M = map(mp.mpf, (T, M))
        bound = 18 / T + 2 * M / T * (mp.log(T / M) + 1)
        assert abs(bubble_finite_difference(T, M)) < bound


@pytest.mark.parametrize("j,e", (("0.01", ".3"), (".0625", ".7"), (".0001", "-.2")))
def test_general_mass_J_by_independent_Feynman_parameter(j, e):
    with mp.workdps(60):
        j, e = map(mp.mpf, (j, e))
        numeric = mp.quad(lambda x: (x + (1 - x) * j) ** (-e), [0, mp.mpf(".1"), 1])
        exact = (1 - j ** (1 - e)) / ((1 - j) * (1 - e))
        assert abs(numeric / exact - 1) < mp.mpf("1e-45")


@pytest.mark.parametrize("T,M", (("16", "1"), ("100", "5")))
@pytest.mark.parametrize("s", ("0", "1"))
def test_uniform_derivative_bound_by_independent_spectral_integral(T, M, s):
    with mp.workdps(30):
        T, M, s = map(mp.mpf, (T, M, s))

        def outer(z):
            v = T / z
            inner = mp.quad(
                lambda x: x * (1 - x) / (x * v + (1 - x) * M - x * (1 - x) * s),
                [0, mp.mpf(".1"), 1],
            )
            return (1 - z) ** mp.mpf("1.5") / (z * (1 - z / T) ** 2) * inner

        value = mp.quad(outer, [0, mp.mpf(".5"), 1])
        assert 0 < value < 4 / T


def test_complete_stationary_field_source_cancellation_independently():
    H, Phi, G, M, J, L, U = sp.symbols("H Phi G M J L U")
    V = M * H**2 / 2 + G * H * Phi**2 / 2 + J * H
    solution = sp.solve(sp.diff(V, H), H)[0]
    eliminated = sp.expand(V.subs(H, solution))
    source_mass = sp.diff(eliminated, Phi, 2).subs(Phi, 0)
    assert sp.factor(source_mass + G * J / M) == 0
    fixed = -G * U / 2
    local = (L - G**2 / M) * U / 2 + source_mass.subs(J, fixed)
    assert sp.factor(local - L * U / 2) == 0
    assert eliminated.subs(Phi, 0) == -(J**2) / (2 * M)


@pytest.mark.parametrize("order", (2, 3, 4))
def test_source_vacuum_cross_terms_have_actual_loop_orders(order):
    h, J1, J2, M = sp.symbols("h J1 J2 M")
    vacuum = sp.expand(-((h * J1 + h * h * J2) ** 2) / (2 * M))
    expected = {2: -J1 * J1 / (2 * M), 3: -J1 * J2 / M, 4: -J2 * J2 / (2 * M)}
    assert sp.factor(vacuum.coeff(h, order) - expected[order]) == 0
    assert vacuum.coeff(h, order) != 0


@pytest.mark.parametrize(
    "m,Y,Q",
    ((2, 0, 144), (3, sp.Rational(1, 7), 144), (10**5, sp.Rational(1, 10**9), 100)),
)
def test_exact_tadpole_bound_scaling(m, Y, Q):
    d = tadpole.enclosure(m, Y, Q)
    doubled = tadpole.enclosure(m, 2 * Y, Q)
    assert (
        doubled["complete_tadpole_absolute_upper"]
        == 2 * d["complete_tadpole_absolute_upper"]
    )
    assert d["complete_tadpole_absolute_upper"] == sum(
        d[k]
        for k in (
            "leading_tadpole_absolute_upper",
            "subleading_tadpole_absolute_upper",
            "finite_tadpole_remainder_absolute_upper",
        )
    )


@pytest.mark.parametrize(
    "m,Y,M,Q",
    (
        (2, 0, 1, 144),
        (4, sp.Rational(1, 7), 4, 144),
        (10**5, sp.Rational(1, 10**9), 10**6, 100),
    ),
)
def test_exact_heavy_bubble_bound_scaling(m, Y, M, Q):
    d = bubble.enclosure(m, Y, M, Q)
    doubled = bubble.enclosure(m, 2 * Y, M, Q)
    assert (
        doubled["complete_finite_bubble_absolute_upper"]
        == 2 * d["complete_finite_bubble_absolute_upper"]
    )
    n = d["dyadic_T_over_M"]
    assert sp.Integer(2) ** (n - 1) < sp.Rational(4 * m * m, M) <= sp.Integer(2) ** n


def test_actual_reference_bounds_retain_every_local_piece():
    d = calibration.data()
    assert all(v is True for v in d["bounds"].values())
    assert d["complete_row_on_shell_mass_reference_absolute_upper"] < sp.Rational(
        1, 10**10
    )
    assert (
        d["complete_row_zero_mass_reference_absolute_upper"]
        > d["local_tadpole_mass_absolute_upper"]
        > 0
    )
    assert (
        d["complete_row_on_shell_mass_reference_absolute_upper"]
        > d["complete_row_zero_mass_reference_absolute_upper"]
    )
    assert d["assigned_H_source_absolute_upper"] > 10**180
    assert d["assigned_H_shift_absolute_upper"] < sp.Rational(1, 10**6)


def test_exact_one_row_frontier_change_and_parent_immutability():
    old = copy.deepcopy(previous.frontier())
    new = audit.frontier()
    assert [a["id"] for a, b in zip(old, new) if a != b] == [audit.TARGET]
    assert previous.frontier() == old
    assert sum(r["status"] == "UNEVALUATED" for r in new) == 2
    assert all(r["status"] != "COMPLETE" for r in new)
    for a, b in zip(old, new):
        if a["id"] != audit.TARGET:
            assert a == b


def test_controls_keep_remaining_scientific_obligations_open():
    controls = audit.controls()
    assert len(controls) == 9
    assert controls["vacuum_source_cross_terms_remain_open"] is True
    assert controls["other_matching_and_canonical_terms_remain_open"] is True
    assert controls["original_P8_not_closed"] is True


def bubble_finite_difference(T, M):
    def f(z):
        k = z / T
        j = M * z / T
        A = 1 / (1 - k) ** 2
        L = -j * mp.log(j) / (1 - j)
        return (
            (1 - z) ** mp.mpf("1.5")
            / z
            * ((3 - 4 * mp.log(2) + 2 * mp.log(z) - mp.log(1 - z)) * h1(k) - A * L)
        )

    return mp.quad(f, [0, mp.mpf(".5"), 1])
