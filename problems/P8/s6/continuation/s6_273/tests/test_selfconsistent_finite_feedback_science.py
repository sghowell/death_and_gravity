"""Independent finite-hybrid source, field, canonical and operator diagnostics."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase as original_scalar
from p8_vacuum_affine_selfconsistent_finite_feedback import (
    audit,
    dynamics,
    geometry,
    homogeneous,
    quantum,
    source,
    symmetry,
)

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")
PACKETS = (
    source.data,
    homogeneous.data,
    symmetry.data,
    geometry.data,
    quantum.data,
    dynamics.data,
)


@pytest.mark.parametrize("packet", PACKETS)
def test_entire_exact_packet_and_all_gates(packet):
    data = packet()
    for name, value in data["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(entry) == 0 for entry in entries), name
    assert all(bool(value) for value in data["gates"].values())


@pytest.mark.parametrize("index", range(500))
def test_every_whole_source_derivative_envelope(index):
    name, expression = list(source.expressions().items())[index]
    bound = source.magnitude(source.q.eliminate_N_primitive(expression))
    assert bound == source.derivative_bounds()[name]
    assert bound.is_Rational is True and 0 <= bound < source.MAJORANT


@pytest.mark.parametrize("index", range(63))
def test_every_componentwise_jet_and_its_actual_clock_binding(index):
    key, bound = list(source.jet_data()["bounds"].items())[index]
    function, i, k = key
    assert i + k <= 5
    if function == source.q.j:
        assert bound == source.OFF_CLOCK
    else:
        expression, clock_bound = source.jet_data()["clocks"][key]
        assert source.absolute_interval(expression) == clock_bound
        assert (
            clock_bound + source.OFF_CLOCK
            < bound
            <= 10 * (clock_bound + source.OFF_CLOCK)
        )
    assert bound.is_Rational is True and bound > 0


@pytest.mark.parametrize("order", range(196))
def test_every_full_phase_derivative_growth_ratio(order):
    for amplitude in (quantum.H_AMPLITUDE, quantum.F_AMPLITUDE):
        for parameter_order in range(3):
            ratio = quantum.jet(order + 1, amplitude, parameter_order) / quantum.jet(
                order, amplitude, parameter_order
            )
            assert ratio == 2056 * (order + 1) ** 3 / quantum.RADIUS
            assert ratio < s.Rational(1, 10**9)


@pytest.mark.parametrize("z", source.q.COORDS)
def test_full_constraint_directional_derivative_before_clock_restriction(z):
    shift = s.Symbol("independent_full_invariant_shift", real=True)
    shifted = source.q.CONSTRAINT.subs(z, z + shift)
    assert (
        s.expand(s.diff(shifted, shift).subs(shift, 0) - s.diff(source.q.CONSTRAINT, z))
        == 0
    )
    # Degree-two in invariants is not a Taylor truncation in phase fields.
    assert s.diff(shifted, shift, 3) == 0


@pytest.mark.parametrize("i,k", [(i, k) for i in range(6) for k in range(7 - i)])
def test_extra_normal_derivative_uses_at_most_five_real_time_derivatives(i, k):
    assert i <= 5 and i + k <= 6
    assert all(value < source.SIX_JET for value in source.jet_data()["six"][(i, k)])
    assert (6, 0) not in source.jet_data()["six"]


def reconstruct_fields(coefficients, point):
    """Literal full shell reconstruction, independent of the signed-action code."""
    scalars = [mp.mpf(0) for _ in range(3)]
    vector = mp.matrix(3, 1)
    tensor = mp.matrix(3)
    for axis in range(3):
        e, f = (axis + 1) % 3, (axis + 2) % 3
        plus, cross = mp.matrix(3), mp.matrix(3)
        plus[e, e], plus[f, f] = 1 / mp.sqrt(2), -1 / mp.sqrt(2)
        cross[e, f] = cross[f, e] = 1 / mp.sqrt(2)
        for wave in range(2):
            basis = mp.cos(point[axis]) if wave == 0 else mp.sin(point[axis])
            start = (2 * axis + wave) * 8
            for target, channel in enumerate((0, 1, 4)):
                scalars[target] += coefficients[start + channel] * basis
            tensor += basis * (
                coefficients[start + 2] * plus + coefficients[start + 3] * cross
            )
            for component in range(3):
                vector[component] += coefficients[start + 5 + component] * basis
    return scalars, vector, tensor


@pytest.mark.parametrize("index", range(24))
def test_actual_signed_actions_on_independently_reconstructed_fields(index):
    with mp.workdps(65):
        g = mp.matrix(symmetry.rotations()[index].tolist())
        U = symmetry.representations()[index]
        coefficients = [mp.mpf((7 * i + 3) % 19 - 9) / 17 for i in range(48)]
        transformed = [
            mp.fsum(int(U[i, j]) * coefficients[j] for j in range(48))
            for i in range(48)
        ]
        x = mp.matrix([mp.mpf("0.13"), mp.mpf("-0.27"), mp.mpf("0.39")])
        left = reconstruct_fields(transformed, x)
        right = reconstruct_fields(coefficients, g.T * x)
        assert max(abs(a - b) for a, b in zip(left[0], right[0], strict=True)) < mp.mpf(
            "1e-58"
        )
        assert mp.norm(left[1] - g * right[1]) < mp.mpf("1e-58")
        assert mp.norm(left[2] - g * right[2] * g.T) < mp.mpf("1e-58")
        # A nonlinear equivariance fixture, not a replacement shape solver.
        lhs = mp.exp(2 * left[0][0]) * mp.expm(-left[2])
        rhs = mp.exp(2 * right[0][0]) * g * mp.expm(-right[2]) * g.T
        assert mp.norm(lhs - rhs) < mp.mpf("1e-56")


@pytest.mark.parametrize("index", range(24))
def test_whitening_covariance_implies_orthogonal_canonical_action(index):
    # An independent nontrivial squeeze commuting with this scalar-only action.
    U = symmetry.representations()[index]
    factors = [
        s.Rational(2 + channel, 3)
        for _ in range(6)
        for channel in (0, 1, 2, 2, 4, 5, 5, 5)
    ]
    L = s.diag(*factors)
    S0 = s.diag(L, L.inv())
    physical = s.diag(U, U)
    V0 = S0 * S0.T / 2
    assert physical * V0 * physical.T == V0
    white = S0.inv() * physical * S0
    assert white * white.T == s.eye(96)
    J = s.zeros(48).row_join(s.eye(48)).col_join((-s.eye(48)).row_join(s.zeros(48)))
    assert white * J * white.T == J


@pytest.mark.parametrize("normal", ((1, 0, 0), (0, -1, 0), (0, 0, 1)))
def test_original_scalar_forms_are_radial_not_new_state_choices(normal):
    radial, scale = s.symbols("independent_radial positive_scale", positive=True)
    k = s.Matrix(normal) * radial
    q = k.dot(k) / scale**2
    assert q == radial**2 / scale**2
    for form in (
        original_scalar.canonical_generator(),
        original_scalar.selection_matrix(),
        original_scalar.central_map(),
        original_scalar.tensor_selection_matrix(),
    ):
        directional = form.subs(
            {original_scalar.charts.q: q, original_scalar.a: scale}, simultaneous=True
        )
        radial_form = form.subs(
            {original_scalar.charts.q: radial**2 / scale**2, original_scalar.a: scale},
            simultaneous=True,
        )
        assert directional == radial_form


def test_anisotropic_seed_is_not_made_radial_by_translation_invariance():
    covariance = s.eye(48)
    covariance[0, 0] = covariance[8, 8] = 2  # Equal cosine/sine variance on one axis.
    assert any(U * covariance * U.T != covariance for U in symmetry.representations())


@pytest.mark.parametrize("case", range(3))
def test_independent_full_mixed_energy_canonical_density_force(case):
    forces = homogeneous.quantum_force_chain()["forces"]
    functions = set().union(*(value.atoms(s.Function) for value in forces.values()))
    K = next(
        value
        for value in functions
        if value.func.__name__ == "whole_centered_quantum_energy"
    )
    alpha, p, pm, eta, ph = K.args
    polynomial = (
        (case + 1) * alpha * p**2
        + pm * ph
        + eta * ph**2
        + alpha * eta
        + p * pm * eta
        + eta**2
    )
    a = s.Symbol("independent_a", positive=True)
    PV, PM, PH = s.symbols("independent_PV independent_PM independent_PH", real=True)
    kappa, volume = s.symbols("kappa torus_volume", real=True)
    density_rule = {alpha: s.log(a), p: PV / (3 * a**3), pm: PM / a**3, ph: PH / a**3}
    E = polynomial.subs(density_rule, simultaneous=True)
    adot = a * s.diff(E, PV) / (kappa * volume)
    pv_dot = -a * s.diff(E, a) / (kappa * volume)
    ph_dot = -(10**100) * s.diff(E, eta) / (kappa * volume)
    direct = {
        "alpha": adot / a,
        "p": pv_dot / (3 * a**3) - PV * adot / a**4,
        "pm": -3 * PM * adot / a**4,
        "eta": 10**100 * s.diff(E, PH) / (kappa * volume),
        "ph": ph_dot / a**3 - 3 * PH * adot / a**4,
        "M1": s.diff(E, PM) / (kappa * volume),
    }
    for name, force in forces.items():
        evaluated = force.xreplace({K: polynomial}).doit()
        evaluated = evaluated.subs(s.Symbol("absolute_alpha", real=True), s.log(a))
        evaluated = evaluated.subs(density_rule, simultaneous=True)
        assert s.simplify(direct[name] - evaluated) == 0, name
    assert s.factor(direct["pm"] + 3 * (PM / a**3) * direct["alpha"]) == 0


def test_scalar_center_changes_force_even_when_only_a_quantum_phase():
    q = s.Symbol("independent_homogeneous_coordinate", real=True)
    phase_energy = 7 * q**2 + 3 * q
    full_force = -s.diff(phase_energy, q)
    assert full_force.subs(q, 1) == -17
    assert s.diff(phase_energy - phase_energy, q) == 0
    assert full_force != 0  # Factoring a phase cannot erase the classical force.


def test_weighted_maximum_metric_not_an_unproved_sum_metric():
    weight = quantum.STATE_WEIGHT
    B = s.Matrix([[s.Rational(3, 5), 0], [s.Rational(3, 5), 0]])
    W = s.diag(1, weight)
    A = W.inv() * B * W
    delta = s.Matrix([1, 0])
    transformed = W * A * delta
    assert max(sum(abs(B[i, j]) for j in range(2)) for i in range(2)) == s.Rational(
        3, 5
    )
    assert max(abs(value) for value in transformed) == s.Rational(3, 5)
    assert sum(abs(value) for value in transformed) == s.Rational(6, 5) > 1
    assert "max{" in dynamics.data()["whole_path_space_and_map"]


@pytest.mark.parametrize("time", ("0.01", "0.7", "3"))
@pytest.mark.parametrize("delta", ("0.02", "-0.13"))
def test_noncommuting_unitary_Duhamel_bound_without_exponential_growth(time, delta):
    with mp.workdps(65):
        X, Z = mp.matrix([[0, 1], [1, 0]]), mp.matrix([[1, 0], [0, -1]])
        t, d = mp.mpf(time), mp.mpf(delta)
        A, B = 100 * X + mp.mpf("0.3") * Z, 100 * X + (mp.mpf("0.3") + d) * Z
        U, V = mp.expm(-mp.j * t * A), mp.expm(-mp.j * t * B)
        difference = U - V
        eigenvalues = mp.eighe(difference.H * difference, eigvals_only=True)
        opnorm = mp.sqrt(max(eigenvalues))
        assert mp.norm(U.H * U - mp.eye(2)) < mp.mpf("1e-58")
        assert opnorm <= abs(t * d) + mp.mpf("1e-58")


@pytest.mark.parametrize("epsilon", ("0.01", "-0.01", "0.005"))
def test_full_implicit_average_not_root_at_average_input(epsilon):
    with mp.workdps(60):
        amplitude = mp.mpf(epsilon)

        def delta_root(x):
            z = amplitude * mp.cos(x)
            return (
                -1 + mp.sqrt(1 - mp.mpf("0.4") * (z + mp.mpf("0.2") * z * z))
            ) / mp.mpf("0.2")

        average = mp.quad(delta_root, [0, mp.pi, 2 * mp.pi]) / (2 * mp.pi)
        assert abs(average) < amplitude**2  # NZ<=2, NZZ<=2 on this fixture.
        assert average < -(amplitude**2) / 10  # avg z=0 does not give avg N=1.
        for x in (mp.mpf("0.2"), mp.mpf("1.3")):
            z = amplitude * mp.cos(x)
            n = delta_root(x)
            assert abs(n + mp.mpf("0.1") * n * n + z + mp.mpf("0.2") * z * z) < mp.mpf(
                "1e-55"
            )


@pytest.mark.parametrize("epsilon", ("0.01", "0.003", "-0.007"))
def test_full_density_exponential_keeps_generated_spatial_mean(epsilon):
    with mp.workdps(60):
        e = mp.mpf(epsilon)
        mean = mp.quad(lambda x: mp.exp(3 * e * mp.cos(x)), [0, mp.pi, 2 * mp.pi]) / (
            2 * mp.pi
        )
        assert 1 < mean < 1 + 10 * e * e
        assert abs(mean - 1) > e * e


@pytest.mark.parametrize("phase", (0, s.pi / 2, s.pi, 3 * s.pi / 2))
def test_both_operator_error_and_turnaround_allow_large_state_motion(phase):
    X = s.Matrix([[0, 1], [1, 0]])
    U = s.cos(phase) * s.eye(2) - s.I * s.sin(phase) * X
    psi = U * s.Matrix([1, 0])
    epsilon = dynamics.VOLUME_ERROR
    readout = s.eye(2) + epsilon * s.diag(1, -1)
    mean = s.simplify((psi.H * readout * psi)[0])
    assert 1 - epsilon <= mean <= 1 + epsilon
    assert (1 + source.TIME**2) ** 6 * mean > 1 + epsilon


@pytest.mark.parametrize(
    "time,error,gap",
    (
        (source.TIME, dynamics.VOLUME_ERROR, True),
        (source.TIME, s.Rational(1, 10**300), False),
        (s.Rational(1, 10**2000), dynamics.VOLUME_ERROR, False),
    ),
)
def test_endpoint_gap_requires_volume_error_smaller_than_time_squared(time, error, gap):
    assert bool((1 + time * time) ** 6 * (1 - error) > 1 + error) is gap


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_reject_every_unsupported_input_and_scope_claim(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_full_joint_root_self_map_keeps_nonzero_phase_image():
    b = geometry.bounds()
    assert b["full_center_residual"] > b["homogeneous_center_residual"]
    assert (
        b["full_center_residual"] / 3 + b["contraction"] * source.LAPSE_RADIUS
        < source.LAPSE_RADIUS
    )
    assert b["deltaN_mean"] > geometry.NZ * b["mean_sum"]


def test_exact_admitted_guards_and_unchanged_original_frontiers():
    for time in (-source.TIME, 0, source.TIME):
        assert audit.require_time(time) == time
    assert audit.require_homogeneous_deviation([0] * 5) == (0,) * 5
    assert (
        audit.require_homogeneous_deviation([source.REAL_RADIUS] * 5)
        == (source.REAL_RADIUS,) * 5
    )
    assert audit.require_model(audit.MODEL) == audit.MODEL
    for ordering in ("calibrated_coherent", "Weyl"):
        assert audit.require_ordering(ordering) == ordering
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 129
    assert "homogeneous" in audit.observable()["not_established"]
    assert "OPEN" in audit.observable()["not_established"]


def test_full_reference_not_small_state_displacement_claim():
    assert source.TIME * quantum.bounds()["centered_operator"] > 2
    assert "No long-time core leakage" in dynamics.data()["whole_open_boundary"]
