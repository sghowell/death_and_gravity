"""Native finite inputs to the pinned global two-cone Hadamard proof."""
import pytest
import sympy as sp
from p8_proca_global_hadamard import (
    audit,
    canonical,
    covariance,
    domains,
    green,
    jets,
    recursion,
    riccati,
    transitions,
)


@pytest.mark.parametrize("name",list(audit.residuals()))
def test_exact_identity(name):
    value=audit.residuals()[name]
    entries=list(value) if isinstance(value,sp.MatrixBase) else [value]
    assert all(item==0 for item in entries),name


@pytest.mark.parametrize("name",list(audit.gates()))
def test_explicit_proof_gate(name):
    assert audit.gates()[name] is True,name


def test_audited_bad_inputs():
    assert audit.controls()["rejected_inputs"]==64


@pytest.mark.parametrize("weight",(0,sp.Rational(1,4),sp.Rational(1,2),1))
def test_accepted_convex_weights(weight):
    assert covariance.cutoff_weight(weight)==weight


@pytest.mark.parametrize("radius",(sp.Rational(1,4),1,2,100))
def test_finite_strip_bounds_stay_strictly_positive(radius):
    assert all(value>0 for value in domains.at_radius(radius).values())


def test_actual_second_graph_jet_keeps_matter_entry():
    assert jets.data()["center_Riccati_second_coefficient"][1,1]==401*sp.I/200


def test_second_order_graph_is_not_an_all_frequency_positive_prescription():
    d=jets.data()
    B=-sp.im(d["center_leading_graph"]+d["center_Riccati_second_coefficient"])
    assert B[1,1]==-sp.Rational(201,200)


def test_naive_leading_graph_has_actual_nonzero_defect():
    assert any(value!=0 for value in jets.data()["center_naive_leading_graph_second_defect"])


def test_losing_Riccati_time_derivative_changes_answer():
    assert any(value!=0 for value in recursion.data()["missing_derivative_control_difference"])


def test_wrong_Green_jump_sign_is_minus_identity():
    Omega=canonical.data()["constant_symplectic_form"]
    assert Omega*Omega==-sp.eye(4)
    assert green.data()["advanced_minus_retarded_kernel"]=="E_Q(t,s)=U_density(t,s)*Omega"


def test_unnormalized_packet_map_changes_CCR():
    d=canonical.data()
    value=(d["density_to_packet"]*d["constant_symplectic_form"]*d["density_to_packet"].T).applyfunc(sp.factor)
    assert value.subs(canonical.k,2)==2*d["constant_symplectic_form"]
    assert value.subs(canonical.k,2)!=d["constant_symplectic_form"]


def test_covariance_proof_does_not_assume_A_B_commute():
    d=covariance.data()
    comm=d["real_symmetric_graph_part"]*d["positive_imaginary_graph_part"]-d["positive_imaginary_graph_part"]*d["real_symmetric_graph_part"]
    assert any(sp.factor(v)!=0 for v in comm)


def test_opposite_graph_sign_reverses_CCR():
    d=covariance.data()
    F=sp.conjugate(d["positive_Cauchy_covariance_Gram_factor"])
    C=F*F.conjugate().T
    wrong=(C-C.T+sp.I*d["positive_hbar"]*d["constant_symplectic_form"]/d["positive_action_normalization"]).applyfunc(sp.factor)
    assert wrong==sp.zeros(4)


def test_exact_global_Laurent_lower_orders_are_zero():
    d=canonical.data()["complete_Laurent_coefficients"]
    assert d[-2]==sp.zeros(4)
    assert d[-3]==sp.zeros(4)
    assert d[-1]!=sp.zeros(4)


def test_overlap_leading_maps_are_not_reset_to_identity():
    d=transitions.data()
    assert d["gamma_to_unitary"]["leading_map"]!=sp.eye(4)
    assert d["canonical_packet_to_gamma_velocity"]["leading_map"]!=sp.eye(4)


def test_Riccati_recursion_requires_all_previous_coefficients():
    with pytest.raises(ValueError,match="preceding"):
        riccati.next_coefficient(2,[sp.eye(2)],[sp.eye(4)],lambda value:value,sp.eye(2),sp.eye(2))


@pytest.mark.parametrize("hamiltonian",([], [sp.eye(2)]))
def test_Riccati_recursion_rejects_incomplete_Hamiltonian_domain(hamiltonian):
    with pytest.raises(ValueError,match="rank-four"):
        riccati.next_coefficient(1,[sp.eye(2)],hamiltonian,lambda value:value,sp.eye(2),sp.eye(2))


def test_Riccati_recursion_rejects_wrong_graph_dimension():
    with pytest.raises(ValueError,match="two-field"):
        riccati.next_coefficient(1,[sp.eye(4)],[sp.eye(4)],lambda value:value,sp.eye(2),sp.eye(2))


def test_Sylvester_solver_rejects_wrong_matrix_dimension():
    with pytest.raises(ValueError,match="two-field"):
        riccati.solve_sylvester(sp.eye(4),sp.eye(2),sp.eye(2))


def test_full_four_mode_transport_keeps_diagonal_amplitudes():
    d=recursion.mode_data()
    assert d["retained_diagonal_transport"]!=sp.zeros(4)
    assert all(d["first_four_mode_change"][j,j]==0 for j in range(4))
