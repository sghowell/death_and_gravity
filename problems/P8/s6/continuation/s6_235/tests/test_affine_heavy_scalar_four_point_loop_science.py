"""Independent complete one-loop graphs, parameter domains and finite matching checks."""

import itertools

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_scalar_four_point_loop import amplitude as amp
from p8_vacuum_affine_heavy_scalar_four_point_loop import audit
from p8_vacuum_affine_heavy_scalar_four_point_loop import symmetric as sym
from p8_vacuum_affine_heavy_scalar_one_loop import audit as previous
from p8_vacuum_affine_heavy_scalar_tree_matching import model


def mpq(value):
    value = s.Rational(value)
    return mp.mpf(int(value.p)) / int(value.q)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_exact_scalar_entry(name, value):
    assert s.cancel(value) == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_every_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "external", [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 1, 1), (1, 1, 1, 1)]
)
def test_complete_trace_log_fourth_tensor_against_three_bubble_pairings(external):
    p = s.symbols("p0 p1")
    G = s.Matrix([[2, s.Rational(1, 3)], [s.Rational(1, 3), 3]])
    GH = s.Matrix(
        [[s.Rational(2, 9), s.Rational(1, 30)], [s.Rational(1, 30), s.Rational(1, 7)]]
    )
    f = s.Matrix([v * v for v in p])
    action = (
        s.Rational(3, 20) * (p[0] ** 4 + p[1] ** 4)
        - s.Rational(1, 13) * (f.T * GH * f)[0]
    )
    W = s.hessian(action, p)
    jet = -s.trace(G * W * G * W) / 4
    actual = s.diff(jet, *(p[i] for i in external))
    vertex = lambda i, j, a, b: s.diff(action, p[i], p[j], p[a], p[b])
    i, j, k, l = external
    predicted = 0
    for left, right in (((i, j), (k, l)), ((i, k), (j, l)), ((i, l), (j, k))):
        for a, b, c, d in itertools.product(range(2), repeat=4):
            predicted -= (
                vertex(*left, a, b) * G[a, c] * G[b, d] * vertex(*right, c, d) / 2
            )
    assert s.cancel(actual - predicted) == 0


def test_labelled_derivative_partition_count_and_all_box_assignments():
    partitions = {}
    for route in itertools.permutations(range(4)):
        pair = tuple(sorted([tuple(sorted(route[:2])), tuple(sorted(route[2:]))]))
        partitions[pair] = partitions.get(pair, 0) + 1
    assert len(partitions) == 3 and set(partitions.values()) == {8}
    assert s.Rational(8, 16) == s.Rational(1, 2)
    ordered = {(a, b) for a in ("s", "t", "u") for b in ("s", "t", "u") if a != b}
    assert len(ordered) == 6
    A, g, H = s.symbols("A g H")
    product = (A + 2 * g * g * H) ** 2 / 2
    assert s.expand(product - A * A / 2 - 2 * A * g * g * H - 2 * g**4 * H * H) == 0


def test_all_UV_counterterms_and_missing_piece_controls():
    C, g, n, delta = s.symbols("C g n delta")
    channels = (s.Integer(1), s.Integer(2), s.Integer(3))
    pole = delta / (32 * s.pi**2) * sum((C + g * g / (n - v)) ** 2 for v in channels)
    dc = -3 * C * C * delta / (32 * s.pi**2)
    dg = -C * g * delta / (32 * s.pi**2)
    dm = g * g * delta / (32 * s.pi**2)
    gpart = sum(2 * g * dg / (n - v) for v in channels)
    mpart = sum(-g * g * dm / (n - v) ** 2 for v in channels)
    assert s.cancel(pole + dc + gpart + mpart) == 0
    fixture = {C: -2, g: s.Rational(1, 3), n: 10, delta: 1}
    for missing in (dc, gpart, mpart):
        assert s.cancel((pole + dc + gpart + mpart - missing).subs(fixture)) != 0


def test_independent_parameter_polynomials_and_ordered_Jacobians():
    z, x, y, n, e, S, T = s.symbols("z x y n e S T")
    l1, l2, h1, h2 = (1 - z) * x, (1 - z) * (1 - x), z * y, z * (1 - y)
    direct = (
        l1 + l2 + n * (h1 + h2) - e * (l1 + l2) * (h1 + h2) - S * l1 * l2 - T * h1 * h2
    )
    expected = (
        1
        - z
        + n * z
        - e * z * (1 - z)
        - S * (1 - z) ** 2 * x * (1 - x)
        - T * z * z * y * (1 - y)
    )
    assert s.expand(direct - expected) == 0
    jac = s.Matrix([l1, h1, l2]).jacobian([z, x, y]).det()
    assert s.expand(jac + z * (1 - z)) == 0
    assert s.expand(expected.subs({e: 0, S: 0, T: 0}) - (1 - z + n * z)) == 0
    assert (
        s.expand(expected.subs({e: 1, S: 0, T: 0}) - (1 - z + n * z) + z * (1 - z)) == 0
    )


@pytest.mark.parametrize("n", [3, 9, 100])
def test_box_order_is_not_symmetric_even_below_cuts(n):
    z = s.Symbol("z")
    Q = (1 - z) ** 2 + n * z
    assert s.expand(Q.subs(z, 1 - z) - Q - (n - 1) * (1 - 2 * z)) == 0
    with mp.workdps(60):
        ds = (
            mp.quad(lambda y: y * (1 - y) ** 3 / (((1 - y) ** 2 + n * y) ** 3), [0, 1])
            / 3
        )
        dt = (
            mp.quad(lambda y: y**3 * (1 - y) / (((1 - y) ** 2 + n * y) ** 3), [0, 1])
            / 3
        )
        assert ds > dt > 0


@pytest.mark.parametrize("n,S", [(9, 5), (100, 16), (1000, 25)])
def test_triangle_delta_function_cut_and_both_forward_box_cuts(n, S):
    with mp.workdps(80):
        nn, ss = mp.mpf(n), mp.mpf(S)
        beta = mp.sqrt(1 - 4 / ss)
        k = (ss - 4) / 2
        a = nn + k
        low, high = (1 - beta) / 2, (1 + beta) / 2

        def density(x):
            h = ss * x * (1 - x) - 1
            z = 2 * h / (2 * h + nn + mp.sqrt(nn * nn + 4 * nn * h))
            return (1 - z) / (nn + 2 * h * (1 - z))

        independent = mp.pi * mp.quad(density, [low, (low + high) / 2, high])
        L = mp.log((a + k) / (a - k))
        assert abs(independent - mp.pi * beta * L / (2 * k)) < mp.mpf("1e-70")
        box0 = mp.pi * beta / 2 * mp.quad(lambda y: 1 / (a - k * y) ** 2, [-1, 1])
        boxu = (
            mp.pi * beta / 2 * mp.quad(lambda y: 1 / (a * a - k * k * y * y), [-1, 1])
        )
        assert abs(box0 - mp.pi * beta / (a * a - k * k)) < mp.mpf("1e-70")
        assert abs(boxu - mp.pi * beta * L / (2 * a * k)) < mp.mpf("1e-70")


@pytest.mark.parametrize("n,S", [(9, 5), (100, 16), (1000, 25)])
def test_full_loop_cut_against_full_rational_tree_angular_quadrature(n, S):
    with mp.workdps(80):
        nn, ss = mp.mpf(n), mp.mpf(S)
        gg = mp.mpf(1) / 49
        dd = nn - 2
        contact = -gg * (3 / dd - 2 / dd**2)
        fixed = contact + gg / (nn - ss)
        beta = mp.sqrt(1 - 4 / ss)
        k = (ss - 4) / 2
        a = nn + k
        logarithm = mp.log((a + k) / (a - k))
        cuts = (
            mp.pi * beta,
            mp.pi * beta * logarithm / (2 * k),
            mp.pi * beta / (a * a - k * k),
            mp.pi * beta * logarithm / (2 * a * k),
        )
        loop = (
            fixed * fixed * cuts[0] / 2
            + 2 * fixed * gg * cuts[1]
            + gg * gg * (cuts[2] + cuts[3])
        ) / (16 * mp.pi**2)
        tree = lambda y: fixed + gg / (a - k * y) + gg / (a + k * y)
        optical = beta / (64 * mp.pi) * mp.quad(lambda y: tree(y) ** 2, [-1, 1])
        assert abs(loop - optical) < mp.mpf("1e-70")
        assert optical > 0
        wrong = loop - gg * gg * cuts[3] / (16 * mp.pi**2)
        assert abs(wrong - optical) > mp.mpf("1e-20")


@pytest.mark.parametrize("n", [3, 9, 100, 10000])
def test_complete_zero_external_integrals_and_mass_derivative(n):
    with mp.workdps(80):
        nn = mp.mpf(n)
        first = mp.quad(lambda y: (1 - y) / (1 + (nn - 1) * y), [0, 1])
        second = mp.quad(lambda y: y * (1 - y) / (1 + (nn - 1) * y) ** 2, [0, 1])
        expected1 = (nn * mp.log(nn) - nn + 1) / (nn - 1) ** 2
        expected2 = ((nn + 1) * mp.log(nn) - 2 * (nn - 1)) / (nn - 1) ** 3
        assert abs(first - expected1) < mp.mpf("1e-70")
        assert abs(second - expected2) < mp.mpf("1e-70")
        assert 0 < second < first / nn
        assert first < mp.log(nn) / (nn - 1)
        deriv = -mp.diff(lambda v: (v * mp.log(v) - v + 1) / (v - 1) ** 2, nn)
        assert abs(second - deriv) < mp.mpf("1e-70")


def test_independent_full_constant_background_two_field_hessian():
    phi, H, k, n, g, C = s.symbols("phi H k n g C")
    V = phi * phi / 2 + n * H * H / 2 - g * H * phi * phi / 2 - C * phi**4 / 24
    full = s.hessian(V, (phi, H)) + k * s.eye(2)
    full = full.subs(H, g * phi * phi / (2 * n))
    reduced = s.cancel(full.det() / (k + n))
    fourth = s.diff(s.log(reduced / (k + 1)) / 2, phi, 4).subs(phi, 0)
    w = -(C + g * g / n) / 2 - g * g / (k + n)
    assert s.cancel(fourth + 6 * w * w / (k + 1) ** 2) == 0


@pytest.mark.parametrize("A,B", [(s.Rational(2, 3), s.Rational(4, 3)), (2, 1), (10, 3)])
def test_complete_angular_integrals_against_independent_quadrature(A, B):
    with mp.workdps(80):
        aa, bb = mpq(A), mpq(B)
        j1 = 2 * mp.atan(mp.sqrt(bb / aa) / 2) / mp.sqrt(aa * bb)
        j2 = 1 / (2 * aa * (aa + bb / 4)) + mp.atan(mp.sqrt(bb / aa) / 2) / (
            aa * mp.sqrt(aa * bb)
        )
        direct1 = mp.quad(lambda x: 1 / (aa + bb * (x - mp.mpf("0.5")) ** 2), [0, 1])
        direct2 = mp.quad(
            lambda x: 1 / (aa + bb * (x - mp.mpf("0.5")) ** 2) ** 2, [0, 1]
        )
        assert abs(j1 - direct1) < mp.mpf("1e-70")
        assert abs(j2 - direct2) < mp.mpf("1e-70")


def test_full_symmetric_denominator_minima_and_common_linear_lower():
    D, z, v, w = s.symbols("D z v w")
    S = s.Rational(4, 3)
    x, y = (1 + v) / 2, (1 + w) / 2
    qc = (1 - z) ** 2 + (D + 2) * z - S * (1 - z) ** 2 * x * (1 - x)
    qd = qc - S * z * z * y * (1 - y)
    lower = s.Rational(2, 3) + D * z
    assert (
        s.expand(qc - lower - (1 - z) ** 2 * v * v / 3 - s.Rational(2, 3) * (z + z * z))
        == 0
    )
    assert (
        s.expand(
            qd
            - lower
            - (1 - z) ** 2 * v * v / 3
            - z * z * w * w / 3
            - s.Rational(2, 3) * z
            - z * z / 3
        )
        == 0
    )
    assert s.Rational(9, 4) + 924 - 1425 == -s.Rational(1995, 4)
    assert s.Rational(39300, 101) > 375


def symmetric_integrals_log_coordinate(nn):
    ss = mp.mpf(4) / 3

    def all_terms(q):
        z = mp.expm1(q) / (nn - 1)
        jac = mp.exp(q) / (nn - 1)
        a = (1 - z) ** 2 + nn * z
        b = ss * (1 - z) ** 2
        A = a - b / 4
        if not b:
            return z, jac, 1 / A, 1 / A**2
        theta = mp.atan(mp.sqrt(b / A) / 2)
        J1 = 2 * theta / mp.sqrt(A * b)
        J2 = 1 / (2 * A * a) + theta / (A * mp.sqrt(A * b))
        return z, jac, J1, J2

    def triangle(q):
        z, jac, J1, _ = all_terms(q)
        return jac * (1 - z) * J1

    def box0(q):
        z, jac, _, J2 = all_terms(q)
        return jac * z * (1 - z) * J2

    path = [0, 1, mp.log(nn) / 2, mp.log(nn)]
    c = mp.quad(triangle, path)
    d0 = mp.quad(box0, path)
    return c, d0, d0 / (1 - ss / (4 * nn)) ** 2


@pytest.mark.parametrize("n", [9, 100])
def test_symmetric_complete_integral_bounds_at_independent_masses(n):
    with mp.workdps(80):
        nn = mp.mpf(n)
        D = nn - 2
        c, d0, dup = symmetric_integrals_log_coordinate(nn)
        assert 0 < c < mp.log(1 + 3 * D / 2) / D
        assert 0 < d0 < dup
        Q = lambda z, x: (
            (1 - z) ** 2 + nn * z - mp.mpf(4) / 3 * (1 - z) ** 2 * x * (1 - x)
        )
        # Independently integrate both original parameter variables.
        original = mp.quad(
            lambda z: (1 - z) * mp.quad(lambda x: 1 / Q(z, x), [0, 1]), [0, 1]
        )
        assert abs(original - c) < mp.mpf("1e-65")
        assert d0 < mp.log(1 + 3 * D / 2) / D**2


def test_actual_symmetric_loop_value_contact_and_tree_hierarchy():
    with mp.workdps(260):
        n = mpq(model.MASS2)
        D = mpq(model.GAP)
        gg = mpq(model.G2)
        C = mpq(model.CONTACT)
        ss = mp.mpf(4) / 3
        cs, d0, dup = symmetric_integrals_log_coordinate(n)
        b = -mp.quad(lambda x: mp.log(1 - ss * x * (1 - x)), [0, 1])
        A = C + gg / (n - ss)
        low = 3 * (A * A * b / 2 + 2 * A * gg * cs + 2 * gg * gg * d0) / (16 * mp.pi**2)
        high = (
            3 * (A * A * b / 2 + 2 * A * gg * cs + 2 * gg * gg * dup) / (16 * mp.pi**2)
        )
        tree = 4 * gg / (D * D * (3 * D + 2))
        valley = gg * (D - 1) / (6 * D * D * (D + 2))
        assert low < high < 0
        assert -high / tree > mp.mpf("1e190")
        assert -low / (24 * valley) < mp.mpf("1e-6")
        assert -low < mpq(sym.ABS_UPPER)
        assert 0 < b < mp.mpf("0.5")
        assert cs > 375 / D and abs(A) > mp.mpf(19) * gg / (10 * D)
        assert tree > 0 and high + (-high) == 0
        assert high - (-high) < 0


def test_contact_only_classical_comparison_retains_entire_valley_margin():
    phi, H, D, g, eps = s.symbols("phi H D g eps")
    C = -g * g * (3 / D - 2 / D**2)
    q = g * g * (D - 1) / (6 * D * D * (D + 2))
    contact_shift = 24 * eps * q
    potential = (
        phi * phi / 2
        + (D + 2) * H * H / 2
        - g * H * phi * phi / 2
        - (C + contact_shift) * phi**4 / 24
    )
    square = (
        phi * phi / 2
        + (D + 2) * (H - g * phi * phi / (2 * (D + 2))) ** 2 / 2
        + (1 - eps) * q * phi**4
    )
    assert s.cancel(potential - square) == 0
    assert 0 < sym.VALLEY_RATIO_UPPER < s.Rational(1, 10**6) < 1
    # A separate finite heavy source counterterm can move the bare classical minimum.
    j = -g / (32 * s.pi**2)
    assert s.diff(potential + j * H, H).subs({phi: 0, H: 0}) == j
    assert j != 0


@pytest.mark.parametrize("S", [5, 16, model.MASS2 - 1])
def test_exact_first_cut_domain(S):
    assert audit.require_cut(S) == S


@pytest.mark.parametrize(
    "stage",
    [
        "complete_first_loop_integral_representation",
        "complete_UV_counterterm_and_first_cut_check",
        "off_shell_zero_jet_control",
        "separate_symmetric_first_loop_value_matching",
        "contact_only_classical_quartic_margin",
    ],
)
def test_positive_named_scopes(stage):
    assert audit.require_stage(stage) == stage


def test_exact_matching_point_and_all_frozen_frontiers():
    expected = ("V2S-T1-OS4", s.Rational(4, 3), s.Rational(4, 3), s.Rational(4, 3), 1)
    assert audit.require_matching(*expected) == expected
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 91
    assert len(audit.residuals()) == 55 and audit.scalar_entry_count() == 55
    assert len(audit.gates()) == 25 and len(audit.controls()) == 9
    assert audit.rejected_inputs() == 270
    assert "SIX ordered boxes" in amp.data()["ordered_boxes"]
    assert "not real2-to2 scattering" in sym.data()["explicit_new_finite_condition"]
    assert (
        "not a proof of strong coupling" in sym.data()["large_unadjusted_loop_to_tree"]
    )
    assert (
        "CONTACT-ONLY COMPARISON"
        in sym.data()["finite_contact_and_classical_quartic_margin"]
    )
    assert "not a bound" in amp.data()["boundary"]
