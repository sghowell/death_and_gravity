"""Independent full-source, CTP, Gaussian, graph-ancestry and KG diagnostics."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_source_filtration import audit, graphs, source


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_full_exact_identity(name):
    value = audit.residuals()[name]
    assert all(e == 0 for e in value) if isinstance(value, s.MatrixBase) else value == 0


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_every_full_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_scope_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_all_current_parameters_and_unpromoted_frontier():
    args = (source.KAPPA, source.MASS2, source.G, source.A)
    assert audit.require_parameters(*args) == args
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 106
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching()) is True
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        audit.rejected_inputs(),
    ) == (152, 201, 26, 254)
    assert all(v is True for k, v in audit.controls().items() if k != "rejected_inputs")


@pytest.mark.parametrize("case", range(8))
def test_independent_full_inverse_metric_canonical_source_leading_coefficient(case):
    eps = s.Symbol("full_variation", real=True)
    # Diagnostic canonical scales exercise the exact general map; actual constants are pinned separately.
    kappa = s.Integer(case + 3) ** 2
    root = s.sqrt(kappa)
    time = s.Rational(case + 1, 19)
    a = (1 + time * time) ** 2
    background = s.diag(1, -a * a, -a * a, -a * a)
    h = s.Matrix(
        [
            [2 + case, 1, -1, 2],
            [1, case + 3, 2, 1],
            [-1, 2, -case, 3],
            [2, 1, 3, case + 5],
        ]
    )
    pi = s.Matrix([case + 4, 2, -3, 1]) / 7
    entire = background + 2 * eps * h / root
    gradient = s.Matrix([1, 0, 0, 0]) + eps * pi / root
    X = s.cancel((gradient.T * entire.inv() * gradient)[0])
    quotient = s.cancel((X - 1) / eps)
    Xfirst = quotient.subs(eps, 0)
    assert X.subs(eps, 0) == 1
    assert Xfirst == 2 * (pi[0] - h[0, 0]) / root
    assert Xfirst != 0
    assert s.diff(quotient, eps).subs(eps, 0) != 0
    alpha = s.Integer(case + 2)
    # Cancel ONLY the exact order-eight factor; keep the entire inverse and localizer.
    J_over_eps8 = (
        source.G
        * kappa
        * (time + eps / root) ** 2
        * quotient**8
        * s.exp(-alpha * X * X)
        / 2
    )
    leading = s.simplify(J_over_eps8.subs(eps, 0))
    expected = (
        128 * source.G * time * time * s.exp(-alpha) * (pi[0] - h[0, 0]) ** 8 / kappa**3
    )
    assert s.simplify(leading - expected) == 0
    wrong = (
        128
        * source.G
        * time
        * time
        * s.exp(-alpha)
        * (pi[0] - h[0, 0] / 2) ** 8
        / kappa**3
    )
    assert s.simplify(leading - wrong) != 0


@pytest.mark.parametrize("case", range(6))
def test_independent_entire_source_finite_extension_and_generic_clock_boundaries(case):
    d = s.Symbol("clock_displacement", real=True)
    time = s.Rational(case + 1, 17)
    alpha = s.Integer(case + 3)
    X = 1 + d
    whole = source.physical_source(time, X, alpha)
    h = d**8 * s.exp(-alpha * (1 + d) ** 2)
    V = d**1024 / ((1 + d) ** 1024 + d**1024)
    expected = source.G * source.KAPPA * time * time * h / 2 - source.J1 * V
    assert s.cancel(whole - expected) == 0
    assert s.limit(s.cancel(V / d**1024), d, 0) == 1
    reduced = source.G * source.KAPPA * time * time * s.exp(
        -alpha * (1 + d) ** 2
    ) / 2 - source.J1 * d**1016 / ((1 + d) ** 1024 + d**1024)
    assert s.cancel(d**8 * reduced - whole) == 0
    assert (
        reduced.subs(d, 0) == source.G * source.KAPPA * time * time * s.exp(-alpha) / 2
    )
    assert whole.subs(d, 0) == 0
    vacuum = source.physical_source(time, 0, alpha)
    assert not vacuum.has(s.Float)
    assert vacuum == source.physical_source(time, s.S.Zero, alpha)
    assert vacuum == source.G * source.KAPPA * time * time / 2 - source.J1
    # A fixed coherent heavy mean permits a source current seven degrees earlier.
    mean = s.Rational(case + 2, 5)
    current = s.diff(mean * source.G * source.KAPPA * time * time * h / 2, d)
    assert s.limit(
        current / d**7, d, 0
    ) == 4 * mean * source.G * source.KAPPA * time * time * s.exp(-alpha)


@pytest.mark.parametrize("case", range(6))
def test_independent_complete_CTP_oscillator_ordering_noise_and_causal_leg(case):
    frequency = s.Rational(case + 2, 3)
    angles = [0, s.pi / 6, s.pi / 2, 5 * s.pi / 6]
    W = s.Matrix(
        4,
        4,
        lambda i, j: s.expand_complex(
            s.exp(-s.I * (angles[i] - angles[j])) / (2 * frequency)
        ),
    )
    T = s.Matrix(4, 4, lambda i, j: W[i, j] if i >= j else W[j, i])
    AT = s.Matrix(4, 4, lambda i, j: W[j, i] if i >= j else W[i, j])
    jp = s.Matrix(s.symbols("jp0:4", real=True))
    jm = s.Matrix(s.symbols("jm0:4", real=True))
    logF = -(jp.T * T * jp + jm.T * AT * jm - 2 * jm.T * W * jp)[0] / 2
    N = s.Matrix(4, 4, lambda i, j: s.cos(angles[i] - angles[j]) / (2 * frequency))
    R = s.Matrix(
        4, 4, lambda i, j: s.sin(angles[i] - angles[j]) / frequency if i > j else 0
    )
    jd, jc = jp - jm, (jp + jm) / 2
    target = (s.I * jd.T * R * jc - jd.T * N * jd / 2)[0]
    assert s.expand(logF - target) == 0
    cosine = s.Matrix([s.cos(q) for q in angles])
    sine = s.Matrix([s.sin(q) for q in angles])
    assert N == (cosine * cosine.T + sine * sine.T) / (2 * frequency)
    assert s.expand(logF.subs(dict(zip(jm, jp, strict=True)))) == 0
    c = s.Matrix(s.symbols("average0:4", real=True))
    d = s.Matrix(s.symbols("difference0:4", real=True))
    action = (d.T * R * c + s.I * d.T * N * d / 2)[0]
    actual = s.Matrix([s.diff(action, v).subs(dict.fromkeys(d, 0)) for v in d])
    wrong = s.Matrix([s.diff((c.T * R * c)[0] / 2, v) for v in c])
    assert actual == R * c and wrong != actual
    assert s.diff(wrong[0], c[1]) != 0 and s.diff(actual[0], c[1]) == 0
    assert s.expand(target - (s.I * jd.T * R * jc - jd.T * N * jd)[0]) != 0


@pytest.mark.parametrize("case", range(8))
def test_independent_two_parameter_full_density_operator_inverse_and_source_variations(
    case,
):
    x, y = s.symbols("metric_x metric_y", real=True)
    K0 = s.Matrix([[case + 4, 1], [1, case + 5]])
    Kx = s.Matrix([[2, 1], [1, -1]])
    Ky = s.Matrix([[1, -2], [-2, 3]])
    Kxy = s.Matrix([[case + 1, 2], [2, -1]])
    K = K0 + x * Kx + y * Ky + x * y * Kxy
    J0 = s.Matrix([case + 1, 2])
    Jx = s.Matrix([1, case - 2])
    Jy = s.Matrix([-2, case + 3])
    Jxy = s.Matrix([3, -1])
    # Unequal density factors on the two components are varied, not commuted.
    density = s.diag(1 + x + 2 * y, 2 - x + y)
    bare = J0 + x * Jx + y * Jy + x * y * Jxy
    j = density * bare
    zero = {x: 0, y: 0}
    G = K.inv()
    G0 = K0.inv()
    Gxy = G0 * Kx * G0 * Ky * G0 + G0 * Ky * G0 * Kx * G0 - G0 * Kxy * G0
    assert (G.diff(x, y).subs(zero) - Gxy).applyfunc(s.cancel) == s.zeros(2)
    assert Kx * G0 * Ky != Ky * G0 * Kx
    j0, jx, jy, jxy = [f.subs(zero) for f in (j, j.diff(x), j.diff(y), j.diff(x, y))]
    effect = (j.T * G * j)[0] / 2
    expected = (
        jxy.T * G0 * j0
        + jx.T * G0 * jy
        - jx.T * G0 * Ky * G0 * j0
        - jy.T * G0 * Kx * G0 * j0
        + j0.T * Gxy * j0 / 2
    )[0]
    assert s.cancel(effect.diff(x, y).subs(zero) - expected) == 0
    frozen_density = s.diag(1, 2)
    wrong = (bare.T * frozen_density * G * frozen_density * bare)[0] / 2
    assert s.cancel(effect.diff(x, y).subs(zero) - wrong.diff(x, y).subs(zero)) != 0


def independent_graph(degrees, edges, grades=None):
    """Actually check vertex incidence and connectivity, not just a degree formula."""
    incidence = [0] * len(degrees)
    neighbors = [set() for _ in degrees]
    for left, right in edges:
        assert 0 <= left < len(degrees) and 0 <= right < len(degrees)
        incidence[left] += 1
        incidence[right] += 1
        neighbors[left].add(right)
        neighbors[right].add(left)
    reached = {0}
    while True:
        new = reached | set().union(*(neighbors[v] for v in reached))
        if new == reached:
            break
        reached = new
    assert len(reached) == len(degrees)
    assert all(used <= degree for used, degree in zip(incidence, degrees, strict=True))
    external = sum(
        degree - used for used, degree in zip(incidence, degrees, strict=True)
    )
    topo = len(edges) - len(degrees) + 1
    grade = sum(grades or [0] * len(degrees)) + topo
    return external, topo, grade


@pytest.mark.parametrize("cross_light_edges", range(9))
def test_independent_all_two_source_graphs_and_countergraph_contractions(
    cross_light_edges,
):
    edges = [(0, 1)] * (1 + cross_light_edges)
    external, topo, grade = independent_graph((9, 9), edges)
    assert external == 16 - 2 * cross_light_edges and topo == grade == cross_light_edges
    full = graphs.complete_graph((9, 9), len(edges))
    assert (full["external"], full["total_loop_grade"]) == (external, grade)
    assert external + 2 * grade == 16
    assert graphs.weighted_vertex(external, grade) == 14
    if grade > 0:
        assert graphs.weighted_vertex(external, 0) < 14


@pytest.mark.parametrize("self_edges", range(5))
def test_independent_actual_single_source_tadpole_edges_and_weight(self_edges):
    external, topo, grade = independent_graph((9,), [(0, 0)] * self_edges)
    assert external == 9 - 2 * self_edges and grade == topo == self_edges
    assert graphs.weighted_vertex(external, grade) == 7
    if self_edges == 4:
        assert external == 1 and grade == 4


def test_independent_generic_mean_and_undeleted_minimal_gravitational_graphs():
    external, topo, grade = independent_graph((10, 9), [(0, 1)] * 9)
    assert (external, topo, grade) == (1, 8, 8)
    # These three edges are two H lines and a graviton line, with NO source marker.
    external, topo, grade = independent_graph((3, 3), [(0, 1)] * 3)
    assert (external, topo, grade) == (0, 2, 2)
    assert grade < graphs.minimum_loops(0, "pure_nonheavy_source")


@pytest.mark.parametrize("counter_grade", range(1, 8))
def test_independent_countergraded_whole_graph_identity_and_contraction(counter_grade):
    # One source descendant degree9-2r (r<=4); attach an ordinary cubic by one edge.
    r = min(counter_grade, 4)
    degree = 9 - 2 * r
    degrees = (degree, 3)
    edges = [(0, 1)]
    external, topo, total = independent_graph(degrees, edges, (r, counter_grade - r))
    expected = sum(
        d - 2 + 2 * v for d, v in zip(degrees, (r, counter_grade - r), strict=True)
    )
    assert expected == external + 2 * total - 2
    assert graphs.weighted_vertex(external, total) == expected
    assert total == counter_grade and topo == 0


@pytest.mark.parametrize("external", range(17))
def test_independent_integer_loop_search_with_actual_higher_source_vertices(external):
    candidates = []
    for d1, d2 in product(range(9, 28), repeat=2):
        legs = d1 + d2 - external
        if legs >= 2 and legs % 2 == 0:
            internal = legs // 2
            # At least the heavy edge is present; other edges are light contractions.
            L = internal - 1
            if L >= 0 and internal <= min(d1, d2):
                candidates.append(L)
    assert candidates
    assert min(candidates) == graphs.minimum_loops(external, "pure_nonheavy_source")


@pytest.mark.parametrize("external", (0, 2, 4, 6, 8, 10, 12, 14, 16))
def test_independent_shifted_Gaussian_source_coefficient_has_required_hbar_grade(
    external,
):
    background, q, hbar, variance = s.symbols(
        "background q hbar variance", positive=True
    )
    polynomial = s.diff((background + q) ** 16, background, external).subs(
        background, 0
    )
    order = 16 - external
    coefficient = s.expand(polynomial).coeff(q, order)
    moment = (
        s.factorial2(order - 1) * (hbar * variance) ** (order // 2)
        if order
        else s.S.One
    )
    result = coefficient * moment
    assert result != 0
    assert s.degree(result, hbar) == graphs.minimum_loops(
        external, "pure_nonheavy_source"
    )
    if external == 2:
        assert result == 240 * s.factorial2(13) * (hbar * variance) ** 7


def test_independent_heavy_mean_wick_contraction_and_wrong_early_claim():
    x, hbar, variance = s.symbols("x hbar variance", positive=True)
    moment = s.integrate(
        x**8 * s.exp(-x * x / (2 * hbar * variance)), (x, -s.oo, s.oo)
    ) / s.sqrt(2 * s.pi * hbar * variance)
    assert s.simplify(moment - 105 * (hbar * variance) ** 4) == 0
    assert s.degree(s.expand(moment), hbar) == 4


@pytest.mark.parametrize("case", range(6))
def test_independent_finite_Gaussian_density_source_integral(case):
    with mp.workdps(55):
        K = mp.matrix([[case + 3, mp.mpf(1) / 3], [mp.mpf(1) / 3, case + 4]])
        j = mp.matrix([mp.mpf(case + 1) / 9, mp.mpf(2 - case) / 11])
        # Triangular change of variables diagonalizes the FULL positive quadratic form.
        L00 = mp.sqrt(K[0, 0])
        L10 = K[1, 0] / L00
        L11 = mp.sqrt(K[1, 1] - L10 * L10)
        L = mp.matrix([[L00, 0], [L10, L11]])
        transformed = L**-1 * j
        literal = mp.mpf(1)
        for coefficient in transformed:
            literal *= mp.quad(
                lambda z, coefficient=coefficient: mp.exp(-z * z / 2 + coefficient * z),
                [-mp.inf, 0, mp.inf],
            ) / mp.sqrt(2 * mp.pi)
        expected = mp.exp((j.T * (K**-1) * j)[0] / 2)
        assert abs(literal / expected - 1) < mp.mpf("1e-48")
        determinant = K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0]
        assert determinant > 0
        wrong = mp.exp((j.T * K * j)[0] / 2)
        assert abs(literal - wrong) > mp.mpf("1e-4")


@pytest.mark.parametrize("mass,momentum", list(product((1, 3, 11), (0, 2, 7))))
def test_independent_complete_forced_KG_energy_and_retarded_solution(mass, momentum):
    import numpy as np
    from scipy.integrate import solve_ivp

    t = s.Symbol("actual_time", real=True)
    a = (1 + t * t) ** 2
    H = s.diff(a, t) / a
    history = (t + s.Rational(1, 2)) ** 4 * (1 + t)
    velocity = s.diff(history, t)
    omega2 = mass * mass + momentum * momentum / a**2
    J = s.diff(history, t, 2) + 3 * H * velocity + omega2 * history
    E = (velocity * velocity + omega2 * history * history) / 2
    identity = (
        s.diff(E, t)
        + 3 * H * velocity**2
        + H * momentum**2 * history**2 / a**2
        - J * velocity
    )
    assert s.cancel(identity) == 0
    assert (
        history.subs(t, -s.Rational(1, 2)) == 0
        and velocity.subs(t, -s.Rational(1, 2)) == 0
    )
    forcing = s.lambdify(t, J, "numpy")
    exact = s.lambdify(t, history, "numpy")
    exact_velocity = s.lambdify(t, velocity, "numpy")

    def rhs(time, state):
        aa = (1 + time * time) ** 2
        hh = 4 * time / (1 + time * time)
        return [
            state[1],
            float(forcing(time))
            - 3 * hh * state[1]
            - (mass * mass + momentum * momentum / aa**2) * state[0],
        ]

    points = np.linspace(-0.5, 0.5, 81)
    solved = solve_ivp(
        rhs, (-0.5, 0.5), (0.0, 0.0), t_eval=points, rtol=2e-11, atol=2e-13
    )
    assert solved.success
    assert np.max(np.abs(solved.y[0] - exact(points))) < 2e-9
    assert np.max(np.abs(solved.y[1] - exact_velocity(points))) < 2e-8
    with mp.workdps(50):
        force = s.lambdify(t, J, "mpmath")
        h = s.lambdify(t, history, "mpmath")
        hd = s.lambdify(t, velocity, "mpmath")
        for time in (
            mp.mpf("-0.4"),
            mp.mpf("-0.2"),
            mp.mpf(0),
            mp.mpf("0.3"),
            mp.mpf("0.5"),
        ):
            integral = mp.quad(lambda u: abs(force(u)), [-mp.mpf(".5"), time])
            assert abs(h(time)) <= 15 * integral / mass
            assert abs(hd(time)) <= 15 * integral
            assert abs(momentum * h(time) / (1 + time * time) ** 2) <= 15 * integral
    # Dividing only the source by kappa is not the same physical forced equation.
    normalized_only = J / source.KAPPA
    assert s.cancel(J - normalized_only) != 0


@pytest.mark.parametrize(
    "value", (True, False, -1, 1.5, s.Rational(1, 2), "2", None, s.oo)
)
def test_independent_invalid_exact_graph_counts_rejected(value):
    for args in ((value, "pure_nonheavy_source"), (value, "one_heavy_source")):
        with pytest.raises((ValueError, TypeError)):
            graphs.minimum_loops(*args)


@pytest.mark.parametrize("case", range(5))
def test_independent_source_free_determinant_and_spatial_source_zero(case):
    q = s.Symbol("spatial_metric_variation", real=True)
    K = s.Matrix([[3 + q + case, 1], [1, 4 + 2 * q + case]])
    assert s.diff(s.log(K.det()), q).subs(q, 0) != 0
    # Full synchronous inverse time row is fixed despite arbitrary off-diagonal spatial change.
    spatial = s.Matrix([[2 + q, 1, q], [1, 3 + 2 * q, 1], [q, 1, 4 - q]])
    metric = s.zeros(4)
    metric[0, 0] = 1
    metric[1:, 1:] = -spatial
    gradient = s.Matrix([1, 0, 0, 0])
    assert (gradient.T * metric.inv() * gradient)[0] == 1
    assert source.physical_source(s.Rational(case + 1, 9), 1) == 0
    assert K.diff(q) != s.zeros(2)
