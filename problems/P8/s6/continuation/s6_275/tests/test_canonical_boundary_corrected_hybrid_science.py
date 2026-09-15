"""Independent boundary, canonical transport, complete-symbol and hybrid tests."""

from functools import cache

import pytest
import sympy as s
from p8_vacuum_affine_canonical_boundary_corrected_hybrid import (
    audit,
    boundary,
    dynamics,
    geometry,
    quantum,
    source,
)
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_hybrid_core_regulator_comparison import quadratic as raw_quadratic
from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar
from p8_vacuum_affine_selfconsistent_finite_feedback import homogeneous, symmetry
from p8_vacuum_affine_weyl_operator_comparison import (
    derivatives as original_derivatives,
)

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")
PACKETS = (boundary.data, source.data, geometry.data, quantum.data, dynamics.data)


@pytest.mark.parametrize("packet", PACKETS)
def test_all_new_corrected_exact_residuals_and_proof_gates(packet):
    data = packet()
    for name, value in data["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(entry) == 0 for entry in entries), name
    assert all(bool(value) for value in data["gates"].values())


@pytest.mark.parametrize("coefficient", (0, 1, -7, s.Rational(3, 5)))
def test_original_raw_and_prepared_velocity_derivatives_keep_arbitrary_retuning(
    coefficient,
):
    packet = boundary.raw_momenta()
    v, sigma, vd, sd, n = scalar.v, scalar.sigma, scalar.vd, scalar.sd, scalar.n
    potential = coefficient * (3 * n**2 - 2 * n * v + 5 * v**2 + n * sigma)
    raw = packet["whole_original_raw_quadratic_ADM"] + potential
    prepared = packet["whole_original_prepared_quadratic_action"] + potential
    delta_v = s.factor(s.diff(raw, vd) - s.diff(prepared, vd))
    delta_m = s.factor(s.diff(raw, sd) - s.diff(prepared, sd))
    assert s.factor(delta_v + 18 * scalar.tree.H * v - 3 * scalar.tree.ell * sigma) == 0
    assert s.factor(delta_m - 3 * scalar.tree.ell * v) == 0
    assert delta_v.subs(scalar.tree.u, 0) == 3 * sigma / 10
    assert delta_m.subs(scalar.tree.u, 0) == 3 * v / 10


def test_tree_binding_diagnostic_not_generalized_to_actual_profile_root():
    tree = boundary.tree_hessian()
    assert tree["whole_nonzero_unshifted_tree_H2_defect"] != 0
    assert tree["checks"]["complete_tree_boundary_H2_transport"] == 0
    # The ACTUAL nonzero fixed profiles remain visible in the independent
    # source binding; the tree diagnostic never substitutes them there.
    q = source.q
    actual = source.previous.previous.previous.clock_at(q.CONSTRAINT)
    assert actual.has(q.physical.RHO)
    assert s.factor(actual) != 0
    assert source.PROFILE > 0
    assert geometry.remainder()["whole_retained_independent_reference_amplitude"] > 0


@pytest.mark.parametrize("time", (0, s.Rational(1, 5), -s.Rational(1, 7)))
def test_full_boundary_weight_and_time_generator_by_direct_differentiation(time):
    u = boundary.u
    x, y = s.symbols("unit_CCR_configuration_x unit_CCR_configuration_y", real=True)
    a = (1 + u**2) ** 2
    H, ell = 4 * u / (1 + u**2), 1 / (10 * a**3)
    direct = a**3 * (-9 * H * x**2 + 3 * ell * x * y)
    assert s.factor(direct - boundary.canonical_generator()) == 0
    expected = -9 * a**3 * (s.diff(H, u) + 3 * H**2) * x**2
    assert s.factor(s.diff(direct, u) - expected) == 0
    assert s.diff(direct, x, y).subs(u, time) == s.Rational(3, 10)
    # Discarding the derivative because H(0)=0 would lose this term.
    assert s.diff(direct, u).subs(u, 0) == -36 * x**2


@pytest.mark.parametrize(
    "power_x,power_y", tuple((i, j) for i in range(4) for j in range(4 - i))
)
@pytest.mark.parametrize("time", (0, s.Rational(1, 5)))
def test_independent_entire_differential_operator_boundary_conjugation(
    power_x, power_y, time
):
    x, y = s.symbols("independent_x independent_y", real=True)
    u = s.Rational(time)
    a = (1 + u**2) ** 2
    F = -36 * u * (1 + u**2) ** 5 * x**2 + s.Rational(3, 10) * x * y
    assert a > 0
    wave = 1 + x**3 + 2 * x * y + y**2
    phase_factor = s.exp(s.I * F)
    literal = phase_factor * wave
    expected = wave
    for variable, power in ((x, power_x), (y, power_y)):
        for _ in range(power):
            literal = -s.I * s.diff(literal, variable)
            expected = s.expand(
                -s.I * s.diff(expected, variable) + s.diff(F, variable) * expected
            )
    residual = s.expand(s.expand(literal) / phase_factor - expected)
    assert s.simplify(residual) == 0
    if power_x + power_y == 1:
        unshifted = (-s.I) ** (power_x + power_y) * s.diff(wave, x, power_x, y, power_y)
        assert s.expand(expected - unshifted) != 0


@pytest.mark.parametrize("seed_id", range(8))
def test_generic_positive_full_covariance_requires_both_boundary_momenta(seed_id):
    k = s.Integer(seed_id + 1)
    A = s.Matrix([[1, k, 2, -1], [0, 2, -k, 3], [2, 1, 3, 1], [-1, 0, k, 2]])
    V = A * A.T + s.eye(4)
    S = s.eye(4)
    S[2, 1] = S[3, 0] = s.Rational(3, 10)
    corrected = S * V * S.T
    assert V[0, 0] > 0 and corrected[0, 3] - V[0, 3] == 3 * V[0, 0] / 10
    assert S * phase.OMEGA * S.T == phase.OMEGA
    assert corrected != V and s.det(corrected) == s.det(V)
    # A scalar global phase would preserve EVERY covariance entry.
    assert corrected[0, 3] != V[0, 3]


@pytest.mark.parametrize("seed_id", range(5))
def test_independent_time_dependent_full_canonical_vector_field_and_wrong_sign(seed_id):
    u = s.Symbol("independent_reference_time", real=True)
    C = s.Matrix(
        [[-72 * u * (1 + u**2) ** 5, s.Rational(3, 10)], [s.Rational(3, 10), 0]]
    )
    S = s.eye(4)
    S[2:, :2] = C
    G = s.Matrix(
        [[2, 1, -3, seed_id], [1, 4, 2, 1], [-3, 2, 5, -2], [seed_id, 1, -2, 7]]
    )
    extra = s.zeros(4)
    extra[:2, :2] = C.diff(u)
    standard = S.T * G * S + extra
    raw = S.diff(u) * S.inv() + S * phase.OMEGA * standard * S.inv()
    assert (raw - phase.OMEGA * G).applyfunc(s.factor) == s.zeros(4)
    missing = S.diff(u) * S.inv() + S * phase.OMEGA * (standard - extra) * S.inv()
    wrong = S.diff(u) * S.inv() + S * phase.OMEGA * (standard - 2 * extra) * S.inv()
    assert (missing - phase.OMEGA * G).applyfunc(s.factor) != s.zeros(4)
    assert (wrong - phase.OMEGA * G).applyfunc(s.factor) != s.zeros(4)


@pytest.mark.parametrize("group_id", range(24))
def test_each_actual_cubic_action_preserves_complete_boundary_and_corrected_covariance(
    group_id,
):
    U = symmetry.representations()[group_id]
    C = boundary.group_contacts()["whole_same_shear_on_all_six_scalar_real_waves"]
    assert U * C * U.T == C
    # The generic original PREPARED covariance includes all q-p/internal rows.
    V = symmetry.radial_covariance()
    diagonal = s.diag(U, U)
    assert diagonal * V * diagonal.T == V
    # Commutation of the shear and group, rather than a new averaged state,
    # implies invariance of Sb V Sb^T. Check the nontrivial block directly.
    transformed_qp = V[:48, :48] * C.T + V[:48, 48:]
    assert U * transformed_qp * U.T == transformed_qp


@pytest.mark.parametrize("time", (0, source.TIME, -source.TIME))
def test_actual_small_time_boundary_rows_and_generator_fit_corrected_bounds(time):
    old = geometry.field.field_bounds()
    v = 4 * old["v_A2"] / (1 + geometry.P) ** 2
    sigma = old["M1_A1"] / (1 + geometry.P)
    C = boundary.matrix().subs(boundary.u, time)
    extra_v = abs(C[0, 0]) * v + abs(C[0, 1]) * sigma
    extra_m = abs(C[1, 0]) * v
    assert extra_v < old["Pi_v_A0"]
    assert extra_m < old["delta_Pi_M_A0"]
    assert old["Pi_v_A0"] + extra_v < geometry.field_rows()["Pi_v_A0"]
    assert old["delta_Pi_M_A0"] + extra_m < geometry.field_rows()["delta_Pi_M_A0"]
    assert (
        abs(s.diff(boundary.matrix()[0, 0], boundary.u).subs(boundary.u, time)) / 2
        < 1000
    )


@pytest.mark.parametrize("radius", (source.RADIUS, source.ANALYSIS_RADIUS))
def test_both_complete_domains_retain_all_changed_trace_and_matter_contacts(radius):
    b = geometry.bounds(radius)
    before = geometry.field.field_bounds()
    assert b["fields"]["Pi_v_A0"] == 16 * radius * before["Pi_v_A0"]
    assert b["fields"]["delta_Pi_M_A0"] == 16 * radius * before["delta_Pi_M_A0"]
    assert (
        b["whole_boundary_time_generator_amplitude"]
        == 10**6 * geometry.KAPPA * b["v"] ** 2
    )
    assert b["D0"] == sum(
        b["full_D" + key] for key in ("trace", "shape", "vector", "matter")
    )
    assert len(b["images"]) == len(b["averages"]) == 12
    assert b["whole_center_residual"] > geometry.NZ * b["sum_images"]
    assert b["contraction"] < s.Rational(1, 20)
    assert (
        b["whole_center_residual"] / 3 + b["contraction"] * source.LAPSE_RADIUS
        < source.LAPSE_RADIUS
    )
    assert b["fields"]["eta_A0"] > 0 and b["full_Dvector"] > 0
    assert b["curvature_remainder"] > 0 and b["shear"] > 0


@pytest.mark.parametrize(
    "name",
    ("trace_shear", "matter", "potential", "spatial", "vector", "full_implicit_lapse"),
)
def test_every_complete_raw_quadratic_row_is_retained_and_rebounded(name):
    old = raw_quadratic.bounds()["whole_quadratic_row_budgets"]
    corrected = geometry.remainder()["whole_corrected_positive_full_quadratic_rows"]
    assert corrected[name] == 4 * old[name] > 0


def test_entire_remainder_and_time_generator_not_erased_by_free_Hessian_identification():
    r = geometry.remainder()
    ratio = source.RADIUS / source.ANALYSIS_RADIUS
    large = geometry.bounds(source.ANALYSIS_RADIUS)
    assert r["whole_complete_degree_at_least_three_remainder"] == large[
        "whole_centered_amplitude"
    ] * ratio**3 / (1 - ratio)
    assert r["whole_corrected_centered_interaction_amplitude"] == sum(
        r[key]
        for key in (
            "whole_raw_constrained_quadratic_amplitude",
            "whole_retained_independent_reference_amplitude",
            "whole_retained_boundary_time_generator_amplitude",
            "whole_complete_degree_at_least_three_remainder",
        )
    )
    assert 0 < r["whole_corrected_centered_interaction_amplitude"] < source.H_AMPLITUDE
    assert 0 < r["whole_corrected_centered_volume_amplitude"] < source.F_AMPLITUDE


@cache
def independent_radial_partitions():
    bell, sums = [[s.Integer(1)]], [s.Integer(1)]
    for n in range(1, 207):
        row = [s.Integer(0)] * (n + 1)
        for k in range(1, n + 1):
            first = bell[n - 1][k - 1]
            second = bell[n - 2][k - 1] if n >= 2 and k - 1 < len(bell[n - 2]) else 0
            row[k] = 4 * first + 2 * (n - 1) * second
        bell.append(row)
        sums.append(sum(row[k] * 256**k * s.factorial(k) ** 2 for k in range(n + 1)))
    return sums


@pytest.mark.parametrize("order", range(207))
def test_all207_radial_bounds_by_independent_Bell_recurrence(order):
    actual = independent_radial_partitions()[order]
    assert actual == quantum.radial_partition(order)
    assert actual <= 2048**order * s.factorial(order) ** 3


@pytest.mark.parametrize("order", range(206))
def test_all206_full_phase_ratios_for_both_amplitudes_and_live_derivatives(order):
    for amplitude in (source.H_AMPLITUDE, source.F_AMPLITUDE):
        for j in range(3):
            ratio = quantum.jet(order + 1, amplitude, j) / quantum.jet(
                order, amplitude, j
            )
            assert ratio == 2056 * (order + 1) ** 3 / source.RADIUS
            assert ratio < s.Rational(1, 10**9)


@pytest.mark.parametrize("order", range(197))
def test_all197_earlier_universal_phase_coefficients_retain_their_exact_values(order):
    assert quantum.coefficient(order) == original_derivatives.coefficient(order)


@pytest.mark.parametrize("degree", (0, 2, 4, 6, 8, 10, 12, 14))
def test_exact_sixth_antiheat_polynomial_including_nonzero_remainder(degree):
    x = s.Symbol("independent_heat_coordinate", real=True)
    f = x**degree
    powers = [f]
    for _ in range(7):
        powers.append(s.diff(powers[-1], x, 2) / 4)
    finite = sum((-1) ** j * powers[j] / s.factorial(j) for j in range(7))
    heat, term = finite, finite
    for j in range(1, degree // 2 + 1):
        term = s.diff(term, x, 2) / 4
        heat += term / s.factorial(j)
    remainder = -powers[7] / s.factorial(7)
    assert s.expand(f - heat - remainder) == 0
    if degree == 14:
        assert remainder != 0 and s.expand(f - heat) != 0


@pytest.mark.parametrize("amplitude", (source.H_AMPLITUDE, source.F_AMPLITUDE))
@pytest.mark.parametrize("order", range(3))
def test_both_complete_operator_and_force_packets_keep_entire_heat_remainder(
    amplitude, order
):
    b = quantum.operator_packet(amplitude, order)
    terms = [
        96**j
        * 2
        * amplitude
        * 2056 ** (2 * j)
        * s.factorial(2 * j) ** 3
        / source.RADIUS ** (2 * j)
        * s.factorial(order)
        / source.HOMOGENEOUS_CAUCHY**order
        / (4**j * s.factorial(j))
        for j in range(8)
    ]
    assert b["whole_all_eight_heat_terms"] == terms
    assert (
        b["whole_both_complete_operator_bound"] == sum(terms[:7]) + 10**112 * terms[7]
    )
    assert (
        b["whole_Weyl_minus_calibrated_operator_difference"]
        == sum(terms[2:7]) + 10**112 * terms[7]
    )
    assert b["whole_retained_Weyl_heat_remainder"] > 0


def test_full_homogeneous_canonical_force_chain_is_retained_for_new_energy():
    data = homogeneous.quantum_force_chain()
    assert set(data["forces"]) == {"alpha", "p", "pm", "eta", "ph", "M1"}
    assert all(s.simplify(value) == 0 for value in data["checks"].values())
    assert homogeneous.density().has(source.q.j, source.q.eta, source.q.ph)
    assert source.q.primitive in homogeneous.density().atoms(s.Function)
    assert homogeneous.clock_residuals()["p"] == -source.q.physical.PRESSURE


def test_weighted_maximum_coupled_map_uses_corrected_not_archived_generator():
    b, q = dynamics.existence_bounds(), quantum.bounds()
    T, weight = source.TIME, quantum.STATE_WEIGHT
    matrix = s.Matrix(
        [
            [
                T * (homogeneous.LIPSCHITZ + q["quantum_force_Lipschitz"]),
                2 * T * q["quantum_force"] / weight,
            ],
            [weight * 5 * T * q["each_first_homogeneous_derivative"], 0],
        ]
    )
    assert b["classical_contraction_row"] == sum(matrix.row(0))
    assert b["quantum_contraction_row"] == sum(matrix.row(1))
    assert max(sum(matrix.row(i)) for i in range(2)) < s.Rational(1, 10**50)
    assert b["self_map_deviation"] < source.REAL_RADIUS / 2
    assert b["physical_volume_error"] < dynamics.VOLUME_ERROR


@pytest.mark.parametrize("comparison", ("cutoff", "ordering"))
def test_independent_two_by_two_coupled_stability_solution(comparison):
    b = dynamics.comparison_bounds()
    a, cross, state = b["whole_coupled_a_b_c"]
    M = s.Matrix([[1 - a, -cross], [-state, 1]])
    force = b["whole_direct_" + comparison + "_force"]
    error_name = (
        "whole_direct_cutoff_unitary_comparison"
        if comparison == "cutoff"
        else "whole_direct_ordering_unitary_comparison"
    )
    rhs = s.Matrix([source.TIME * force, b[error_name]])
    exact = M.inv() * rhs
    assert exact[0] == b["whole_coupled_" + comparison + "_five_coordinate_difference"]
    assert (
        exact[1]
        == b["whole_coupled_" + comparison + "_phase_factored_state_difference"]
    )
    assert M.det() > s.Rational(99, 100)
    assert cross > 0 and state > 0 and force > 0


def test_physical_turnaround_core_probability_and_regulator_scopes_remain_distinct():
    e, c = dynamics.existence_bounds(), dynamics.comparison_bounds()
    assert e["endpoint_gap"] > 5 * source.TIME**2
    assert e["minimum_location_margin"] > 0 and dynamics.MINIMUM_RADIUS < source.TIME
    assert 0 < c["whole_evolved_positive_coherent_POVM_leakage"] < dynamics.LEAKAGE
    assert c["whole_coupled_cutoff_physical_volume_difference"] < dynamics.CUTOFF_VOLUME
    assert (
        c["whole_coupled_ordering_physical_volume_difference"] < dynamics.ORDER_VOLUME
    )
    assert audit.STATE != audit.previous.STATE and audit.MODEL != audit.previous.MODEL
    assert audit.require_prepared_state(audit.PREPARED_STATE) == audit.previous.STATE


def test_historical_statuses_preserved_but_all_six_physical_claims_qualified():
    old = audit.previous.matching()
    assert audit.matching()[:-1] == old
    assert audit.frontier() == audit.previous.frontier()
    assert len(audit.frontier()) == 9
    qualified = audit.qualifications()
    assert len(qualified) == 6
    for row in qualified:
        original = next(item for item in old if item["id"] == row["id"])
        assert original["status"] == row["archived_status_preserved_not_reendorsed"]
        assert "refuted" in row["current_physical_qualification"]
    assert len(audit.matching()) == 131
    assert "OPEN" in audit.observable()["not_established"]


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_invalid_or_overclaimed_corrected_inputs_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_warmed_analysis_and_operator_caches_do_not_accept_bool_or_float_aliases():
    for radius in (source.RADIUS, source.ANALYSIS_RADIUS):
        geometry.bounds(radius)
    for amplitude in (source.H_AMPLITUDE, source.F_AMPLITUDE):
        for order in range(3):
            quantum.operator_packet(amplitude, order)
    for order in range(207):
        quantum.radial_partition(order)
    for wrong in (True, False, 0.0, 1.0, s.Float(1), s.Float(0)):
        with pytest.raises((TypeError, ValueError)):
            quantum.radial_partition(wrong)
        with pytest.raises((TypeError, ValueError)):
            quantum.operator_packet(source.H_AMPLITUDE, wrong)
        with pytest.raises((TypeError, ValueError)):
            geometry.bounds(wrong)
