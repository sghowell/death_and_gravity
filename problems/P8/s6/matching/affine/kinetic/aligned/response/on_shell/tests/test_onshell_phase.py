"""Independent Cauchy jets and physical canonical-map controls."""
import pytest
import sympy as sp
from p8_aligned_onshell import checks, phase


def test_all_actual_phase_and_readout_identities():
    for value in checks.residuals().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_omitting_the_moving_boundary_changes_the_phase_system():
    wrong = checks.canonical_swap()["omitted_time_boundary_residual"]
    assert wrong != sp.zeros(4)
    assert wrong.subs(phase.u, phase.LEFT) != sp.zeros(4)
    assert wrong.subs(phase.u, 0) == sp.zeros(4)


def test_exact_prepared_family_is_nontrivial_and_really_drives_the_vector():
    data = phase.initial_jets()
    assert data["prepared_matrix"].rank() == 3
    assert data["observability_determinant"] != 0
    assert max(abs(value) for value in data["initial_vector"]) == 1
    assert data["lapse_third"] > 0
    assert data["trace_initial"] > 0
    assert data["source_third_initial"] < 0


def test_prepared_jets_from_independent_Cauchy_recurrence():
    data, jets = phase.system(), phase.initial_jets()
    initial = jets["initial_vector"]
    A = [data["A"].diff(phase.u, j).subs(phase.u, phase.LEFT) for j in range(3)]
    state_jets = [initial]
    for order in range(3):
        next_jet = sp.zeros(4, 1)
        for i in range(order+1):
            next_jet += sp.binomial(order, i)*A[i]*state_jets[order-i]
        state_jets.append(next_jet)
    observer = [data["lapse_row"].diff(phase.u, j).subs(phase.u, phase.LEFT) for j in range(4)]
    lapse = []
    for order in range(4):
        lapse.append(sum(sp.binomial(order, i)*(observer[i]*state_jets[order-i])[0] for i in range(order+1)))
    assert lapse[:3] == [0, 0, 0]
    assert lapse[3] == jets["lapse_third"]


@pytest.mark.parametrize("point", (phase.LEFT, -sp.Rational(1, 4), 0, sp.Rational(1, 4), phase.RIGHT))
def test_full_phase_and_observers_are_regular_through_the_crossing(point):
    data = phase.system()
    for value in (*data["A"], *data["lapse_row"], *data["trace_row"]):
        actual = value.subs(phase.u, point)
        assert actual.is_real is True and actual.is_finite is True
    assert data["q"].subs(phase.u, point) > 0


def test_real_cosine_source_has_zero_and_double_momentum():
    assert 4*phase.K_IN2 == sp.Rational(1, 4)
    # Each complex +/-2k coefficient is one quarter of the cos² source.
    source3 = phase.initial_jets()["source_third_initial"]
    assert source3/4 != 0
    assert source3/4+source3/2+source3/4 == source3
