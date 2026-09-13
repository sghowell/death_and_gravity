"""Independent complete nonscalar curvature, Ward, force-domain and inverse diagnostics."""

from functools import cache
from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_two_mass_nonscalar_inverse import audit, geometry, local


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_complete_exact_identity(name):
    value = audit.residuals()[name]
    assert all(x == 0 for x in value) if isinstance(value, s.MatrixBase) else value == 0


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_every_complete_written_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_scope_is_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_complete_current_parameters_and_unpromoted_frontier():
    g = geometry.reference
    assert audit.require_parameters(g.PROCA_MASS2, g.HEAVY_MASS2, g.KAPPA) == (
        g.PROCA_MASS2,
        g.HEAVY_MASS2,
        g.KAPPA,
    )
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 105
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching()) is True
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        audit.rejected_inputs(),
    ) == (108, 522, 29, 244)
    assert all(v is True for k, v in audit.controls().items() if k != "rejected_inputs")


@cache
def independent_weyl(entries, frequency, momentum):
    """Build Christoffels, then Riemann, Ricci and the full Weyl tensor."""
    h = s.zeros(4)
    h[1:, 1:] = -s.Matrix(3, 3, entries)
    metric = s.diag(1, -1, -1, -1)
    derivative = (frequency, 0, 0, s.I * momentum)
    connection = {
        (a, b, c): sum(
            metric[a, d]
            * (
                derivative[b] * h[d, c]
                + derivative[c] * h[d, b]
                - derivative[d] * h[b, c]
            )
            / 2
            for d in range(4)
        )
        for a, b, c in product(range(4), repeat=3)
    }
    curvature = {
        (a, b, c, d): sum(
            metric[a, e]
            * (
                derivative[c] * connection[e, d, b]
                - derivative[d] * connection[e, c, b]
            )
            for e in range(4)
        )
        for a, b, c, d in product(range(4), repeat=4)
    }
    ricci = s.Matrix(
        4,
        4,
        lambda b, d: sum(
            metric[a, c] * curvature[a, b, c, d] for a, c in product(range(4), repeat=2)
        ),
    )
    scalar = s.trace(metric * ricci)
    return {
        (a, b, c, d): s.expand(
            curvature[a, b, c, d]
            - (
                metric[a, c] * ricci[b, d]
                - metric[a, d] * ricci[b, c]
                - metric[b, c] * ricci[a, d]
                + metric[b, d] * ricci[a, c]
            )
            / 2
            + scalar * (metric[a, c] * metric[b, d] - metric[a, d] * metric[b, c]) / 6
        )
        for a, b, c, d in product(range(4), repeat=4)
    }


@pytest.mark.parametrize(
    "frequency,momentum",
    [
        (s.Rational(2, 3), s.Rational(5, 7)),
        (s.Rational(-4, 5), s.Rational(3, 2)),
        (s.I * s.Rational(7, 5), s.Rational(8, 3)),
    ],
)
@pytest.mark.parametrize("i,j", list(product(range(2, 6), repeat=2)))
def test_independent_full_christoffel_weyl_all_nonscalar_pairs(
    frequency, momentum, i, j
):
    D, G = geometry.BASIS[i], geometry.BASIS[j]
    CD = independent_weyl(tuple(D), frequency, momentum)
    CG = independent_weyl(tuple(G), frequency, momentum)
    signs = (1, -1, -1, -1)
    literal = s.expand(
        2
        * sum(
            signs[a] * signs[b] * signs[c] * signs[d] * CD[a, b, c, d] * CG[a, b, c, d]
            for a, b, c, d in product(range(4), repeat=4)
        )
    )
    produced = local.weyl_mixed(D, G).subs({local.lam: frequency, local.P: momentum})
    target = (
        (frequency**2 + momentum**2) ** 2
        if i == j and i >= 4
        else frequency**2 * (frequency**2 + momentum**2)
        if i == j
        else 0
    )
    assert s.expand(literal - target) == 0
    assert s.expand(literal - produced) == 0
    if i == j and i < 4:
        assert s.expand(literal - (frequency**2 + momentum**2) ** 2) != 0


@pytest.mark.parametrize("case", range(8))
def test_independent_whole_spatial_ricci_bilinear_with_trace_and_cross_terms(case):
    D = s.Matrix([[case + 1, 2, -1], [2, 3 - case, case + 2], [-1, case + 2, -2]])
    G = s.Matrix([[2, case - 1, 3], [case - 1, -1, 2 - case], [3, 2 - case, case + 4]])
    n = s.Matrix([0, 0, 1])
    p2 = local.P**2
    expected = p2 * (
        -s.trace(D * G) / 2
        + (D * n).dot(G * n)
        - ((n.T * D * n)[0] * s.trace(G) + (n.T * G * n)[0] * s.trace(D)) / 2
        + s.trace(D) * s.trace(G) / 2
    )
    assert s.expand(local.spatial_mixed(D, G) - expected) == 0
    assert s.expand(expected + p2 * s.trace(D * G) / 2) != 0


@pytest.mark.parametrize(
    "axis",
    [
        s.Matrix([0, 0, 1]),
        s.Matrix([1, 2, 2]) / 3,
        s.Matrix([2, 3, 6]) / 7,
        s.Matrix([2, -2, 1]) / 3,
    ],
)
@pytest.mark.parametrize("case", range(6))
def test_independent_coordinate_free_orthogonal_projectors_at_nonaxial_momenta(
    axis, case
):
    Q = s.Matrix([[1 + case, 2, -case], [2, case - 3, 1], [-case, 1, 4 - case]])
    sectors = ("scalar", "vector", "tensor")
    projected = [geometry.project(Q, axis, label) for label in sectors]
    assert sum(projected, s.zeros(3)) == Q
    for i, U in enumerate(projected):
        assert U == U.T
        for j, label in enumerate(sectors):
            assert geometry.project(U, axis, label) == (U if i == j else s.zeros(3))
            if i != j:
                assert s.trace(U * projected[j]) == 0
    vector, tensor = projected[1:]
    assert tensor * axis == s.zeros(3, 1) and s.trace(tensor) == 0
    plane = s.eye(3) - axis * axis.T
    assert plane * vector * plane == s.zeros(3)
    assert (axis.T * vector * axis)[0] == 0
    assert sum(s.trace(U * U) for U in projected) == s.trace(Q * Q)


@pytest.mark.parametrize("case", range(4))
def test_independent_full_O2_operator_commutes_without_scalar_time_symmetry(case):
    rotation = s.Matrix(
        [
            [s.Rational(3, 5), -s.Rational(4, 5), 0],
            [s.Rational(4, 5), s.Rational(3, 5), 0],
            [0, 0, 1],
        ]
    )
    if case % 2:
        rotation = rotation * s.diag(-1, 1, 1)
    action = s.Matrix(
        6,
        6,
        lambda i, j: s.trace(
            geometry.BASIS[i] * rotation * geometry.BASIS[j] * rotation.T
        ),
    )
    whole = s.diag(
        s.Matrix([[2, 3], [5, 7]]), (case + 11) * s.eye(2), (case + 17) * s.eye(2)
    )
    assert (whole * action - action * whole).applyfunc(s.simplify) == s.zeros(6)
    assert whole != whole.T
    if case % 2:
        odd = whole.copy()
        odd[2, 3], odd[3, 2] = 1, -1
        assert odd * action != action * odd


@pytest.mark.parametrize("case", range(6))
@pytest.mark.parametrize("sector", ("tensor", "vector"))
def test_independent_full_weighted_action_variation_and_current_tree(case, sector):
    t, u = local.t, s.Symbol("u", real=True)
    aa = (1 + u * u) ** 2
    HH = s.diff(aa, u) / aa
    qq = s.Rational(case + 1, 3) ** 2 / aa**2
    R = 6 * (s.diff(HH, u) + 2 * HH**2)
    # Diagnostic parameters test a literal formula, not a different admitted physical model.
    mm, nn, ll = s.Rational(case + 2), s.Rational(case + 7), s.Rational(case + 11, 3)
    AA = (5 * mm + nn * (ll - 1)) / 12 - (ll + 2) * R / 72
    bb = (ll + 2) / 60
    f = s.Function("independent_history")(u)
    D0 = s.diff(f, u, 2) + HH * s.diff(f, u)
    density = (
        aa**3 * (AA * (s.diff(f, u) ** 2 - qq * f * f) - bb * (D0 + qq * f) ** 2 / 2)
        if sector == "tensor"
        else aa**3
        * (AA * s.diff(f, u) ** 2 - bb * (D0**2 - qq * s.diff(f, u) ** 2) / 2)
    )
    Euler = (
        sum(
            (-1) ** j * s.diff(s.diff(density, s.diff(f, u, j)), u, j) for j in range(3)
        )
        / aa**3
    )
    classical = -local.C0 * (
        s.diff(f, u, 2) + 3 * HH * s.diff(f, u) + (qq * f if sector == "tensor" else 0)
    )
    produced = local.data()["entire_nonscalar_local_plus_classical_operators"][sector]
    oldf = s.Function("unit_nonscalar_metric")(t)
    produced = (
        produced.subs(
            {
                local.mu2: mm,
                local.n: nn,
                local.ell: ll,
                local.P: s.Rational(case + 1, 3),
            }
        )
        .subs(t, u)
        .subs(oldf.subs(t, u), f)
        .doit()
    )
    fixture = (u + s.Rational(1, 2)) ** 6 * (1 + u + u**3)
    difference = (Euler + classical - produced).subs(f, fixture).doit()
    assert s.cancel(difference) == 0
    derivative_term = -2 * s.diff(AA, u) * s.diff(f, u)
    assert s.cancel(derivative_term.subs(f, fixture).doit()) != 0
    kinetic = (
        local.state.KAPPA
        * aa**3
        * s.diff(2 * f / s.sqrt(local.state.KAPPA), u) ** 2
        / 8
    )
    assert s.cancel(kinetic - aa**3 * s.diff(f, u) ** 2 / 2) == 0


@pytest.mark.parametrize("case", range(4))
def test_independent_unexpanded_noncommuting_metric_volume_and_clock_contraction(case):
    t, x, y, z = geometry.metric.COORDS
    e, g = s.symbols("e g")
    D = geometry.BASIS[2] + (case + 1) * geometry.BASIS[4]
    G = geometry.BASIS[3] + (case + 2) * geometry.BASIS[5]
    assert D * G != G * D
    Q = e * D + g * G
    exponential = s.eye(3) + Q + Q * Q / 2
    determinant = s.expand(exponential.det())
    assert determinant.coeff(e, 1).coeff(g, 1) == 0
    assert determinant.coeff(e, 1).subs(g, 0) == 0
    assert determinant.coeff(g, 1).subs(e, 0) == 0
    mixed_metric = (D * G + G * D) / 2
    assert exponential.diff(e, g).subs({e: 0, g: 0}) == mixed_metric
    assert mixed_metric != s.zeros(3)
    # Full inhomogeneous synchronous metric: spatial Hessian is not set to zero.
    h = s.Matrix(
        [
            [2 + t * t + x * x, t * x, y * z],
            [t * x, 3 + y * y, t * z],
            [y * z, t * z, 4 + z * z + t],
        ]
    )
    metric = s.diag(1, -1, -1, -1)
    metric[1:, 1:] = -(local.a**2) * h
    gradient = s.Matrix([1, 0, 0, 0])
    assert metric * gradient == gradient
    connection0 = s.Matrix(
        4,
        4,
        lambda mu, nu: (
            (
                s.diff(metric[0, nu], geometry.metric.COORDS[mu])
                + s.diff(metric[0, mu], geometry.metric.COORDS[nu])
                - s.diff(metric[mu, nu], t)
            )
            / 2
        ),
    )
    assert connection0[0, :] == s.zeros(1, 4)
    assert connection0[1:, 1:] != s.zeros(3)
    assert gradient.dot(gradient) == 1


@pytest.mark.parametrize("direction", range(2))
@pytest.mark.parametrize("case", range(3))
def test_independent_unexpanded_both_ordered_vector_Ward_contacts_with_nonzero_mean(
    direction, case
):
    m = geometry.metric
    t, z = m.t, m.z
    p = s.Rational(case + 1, 3)
    e, g = s.symbols("detector_parameter source_parameter")
    phase = s.exp(s.I * p * z)
    chi = 1 + t + t**3
    xi = s.zeros(4, 1)
    xi[direction + 1] = chi * phase
    bg = s.zeros(3, 1)
    bg[direction] = s.diff(chi, t) * phase
    qg = s.sqrt(2) * s.I * p * chi * phase * geometry.BASIS[direction + 2]
    bd = s.zeros(3, 1)
    bd[direction] = (1 - t + t * t) / phase
    qd = (1 + t) ** 2 * geometry.BASIS[direction + 2] / phase
    Q = e * qd + g * qg
    spatial = m.a**2 * (s.eye(3) + Q + Q * Q / 2)
    shift = e * bd + g * bg
    entire = s.zeros(4)
    entire[0, 0] = 1 - (shift.T * spatial * shift)[0]
    entire[0, 1:] = -(spatial * shift).T
    entire[1:, 0] = -spatial * shift
    entire[1:, 1:] = -spatial
    mixed = entire.diff(e, g).subs({e: 0, g: 0})
    hD = entire.diff(e).subs({e: 0, g: 0})
    rho, pressure = 2 + t * t, 3 + t + t * t
    E = s.diag(
        -(m.a**3) * rho / 2,
        -m.a * pressure / 2,
        -m.a * pressure / 2,
        -m.a * pressure / 2,
    )
    J = xi.jacobian(m.COORDS)
    transport = E.applyfunc(
        lambda f: sum(xi[j] * s.diff(f, c) for j, c in enumerate(m.COORDS))
    )
    lieE = transport - J * E - E * J.T + s.trace(J) * E
    first = s.trace(hD.T * lieE)
    contact = s.trace(E.T * mixed)
    assert s.cancel(first + contact) == 0
    assert (mixed - m.second((0, bd, qd), (0, bg, qg))).applyfunc(s.cancel) == s.zeros(
        4
    )
    missing_lapse_contact = -E[0, 0] * mixed[0, 0]
    missing_spatial_contact = -s.trace(E[1:, 1:] * mixed[1:, 1:])
    assert s.cancel(missing_lapse_contact) != 0
    assert s.cancel(missing_spatial_contact) != 0
    # Interchange source and detector roles while keeping the opposite Fourier phases.
    xiD = xi / phase**2
    hG = m.first(
        (0, s.zeros(3, 1), (1 + t + t * t) * phase * geometry.BASIS[direction + 2])
    )
    JD = xiD.jacobian(m.COORDS)
    lieh = (
        hG.applyfunc(
            lambda f: sum(xiD[j] * s.diff(f, c) for j, c in enumerate(m.COORDS))
        )
        + JD.T * hG
        + hG * JD
    )
    detector_first = -s.trace(E.T * lieh)
    detector_contact = s.trace(
        E.T
        * m.second(
            m.gauge(xiD),
            (0, s.zeros(3, 1), (1 + t + t * t) * phase * geometry.BASIS[direction + 2]),
        )
    )
    assert s.cancel(detector_first + detector_contact) == 0
    assert s.cancel(detector_first) != 0


@pytest.mark.parametrize("case", range(5))
def test_independent_distinct_advanced_detector_and_retarded_weighted_transpose(case):
    t, u = local.t, s.Symbol("independent_time", real=True)
    lower, upper = -s.Rational(1, 2), s.Rational(1, 2)
    f = (t - lower) ** (case + 2) * (1 + t)
    detector = (upper - t) ** (case + 2) * (2 + t)
    Iplus = s.integrate(detector.subs(t, u), (u, upper, t))
    Iminus = s.integrate(detector.subs(t, u), (u, lower, t))
    IW = s.integrate((local.a**3 * f).subs(t, u), (u, lower, t)) / local.a**3
    lhs = s.integrate(local.a**3 * f * Iplus, (t, lower, upper))
    rhs = -s.integrate(local.a**3 * detector * IW, (t, lower, upper))
    wrong = s.integrate(local.a**3 * f * Iminus, (t, lower, upper))
    assert lhs == rhs and lhs != wrong
    assert Iplus.subs(t, upper) == 0 and Iplus.subs(t, lower) != 0
    assert IW.subs(t, lower) == 0
    assert s.diff(Iplus, t) == detector.expand()
    assert s.cancel(s.diff(IW, t) + 3 * local.H * IW - f) == 0


@pytest.mark.parametrize("case", range(6))
def test_independent_full_compatible_force_and_original_shift_constraint_recovery(case):
    t, u = local.t, s.Symbol("force_time", real=True)
    lower = -s.Rational(1, 2)
    p = s.Rational(case + 1, 7)
    kappa = s.Rational(case + 2) ** 2
    f = (t - lower) ** 5 * (1 + t + t * t)
    ev = s.sqrt(kappa) * f / 2
    eb = (
        -s.sqrt(2)
        * s.I
        * p
        * s.integrate((local.a**3 * ev).subs(t, u), (u, lower, t))
        / local.a**3
    )
    Ev = -ev
    Ebeta = (
        s.sqrt(2)
        * s.I
        * p
        * s.integrate((local.a**3 * ev).subs(t, u), (u, lower, t))
        / local.a**3
    )
    W = lambda q: s.diff(q, t) + 3 * local.H * q
    assert s.cancel(W(eb) + s.sqrt(2) * s.I * p * ev) == 0
    assert s.cancel(W(Ebeta) + s.sqrt(2) * s.I * p * Ev) == 0
    assert s.cancel(Ebeta + eb) == 0
    unprepared = Ebeta + (case + 1) / local.a**3
    assert s.cancel(W(unprepared) + s.sqrt(2) * s.I * p * Ev) == 0
    assert (unprepared + eb).subs(t, lower) != 0
    assert s.cancel(2 * ev / s.sqrt(kappa) - f) == 0


@pytest.mark.parametrize(
    "epsilon",
    (s.Rational(1, 2), s.Rational(1, 8), s.Rational(1, 32), s.Rational(1, 128)),
)
def test_independent_infrared_shift_force_annulus_counterexample(epsilon):
    p = s.Symbol("radial_momentum", positive=True)
    # Angular factors cancel: raw shift amplitude is supported on [eps,2eps].
    raw_squared = s.integrate(p * p, (p, epsilon, 2 * epsilon))
    quotient_squared = s.integrate(1, (p, epsilon, 2 * epsilon))
    assert quotient_squared / raw_squared == s.Rational(3, 7) / epsilon**2
    assert quotient_squared / raw_squared > 1 / (4 * epsilon**2)


@pytest.mark.parametrize("size", (3, 4, 5, 6))
def test_independent_ordered_causal_inverse_and_wrong_order_defect(size):
    F = s.Matrix(size, size, lambda i, j: s.Rational(i + j + 1, i + 2) if i >= j else 0)
    V = s.Matrix(
        size,
        size,
        lambda i, j: s.Rational(2 * i - j + 1, 7 * i + j + 3) if i > j else 0,
    )
    B = F.inv()
    assert F * V != V * F
    E = (s.eye(size) + B * V).inv() * B
    assert (F + V) * E == s.eye(size)
    assert E * (F + V) == s.eye(size)
    assert E == B * (s.eye(size) + V * B).inv()
    wrong = B * (s.eye(size) + B * V).inv()
    assert (F + V) * wrong != s.eye(size)


@pytest.mark.parametrize("scale", (0, 1, 3, 10, 100))
def test_independent_complete_weak_log_weight_and_strict_contraction(scale):
    with mp.workdps(55):
        beta = mp.mpf(scale)
        weight = (4 * beta + 1) ** 2
        complete = (
            mp.quad(lambda x: mp.exp(-weight * x) * (1 - mp.log(x)), [0, 1 / weight, 1])
            if weight > 1
            else mp.quad(lambda x: mp.exp(-x) * (1 - mp.log(x)), [0, 1])
        )
        envelope = 2 / mp.sqrt(weight)
        assert complete < envelope
        assert beta * envelope < mp.mpf("0.5")
        assert abs(mp.mpf("0.5") - beta * envelope - 1 / (2 * (4 * beta + 1))) < mp.mpf(
            "1e-50"
        )
