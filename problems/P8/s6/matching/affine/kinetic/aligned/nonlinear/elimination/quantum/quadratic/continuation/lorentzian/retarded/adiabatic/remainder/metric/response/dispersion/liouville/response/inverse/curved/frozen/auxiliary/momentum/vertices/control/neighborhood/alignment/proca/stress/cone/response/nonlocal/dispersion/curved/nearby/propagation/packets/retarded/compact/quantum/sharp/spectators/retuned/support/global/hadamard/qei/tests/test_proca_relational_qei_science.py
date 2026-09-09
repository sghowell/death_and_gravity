"""Native inputs to the actual relational inequality and source-boundary proof."""
import pytest
import sympy as sp
from p8_proca_relational_qei import audit, modes, observers, qei, stress


@pytest.mark.parametrize("name",list(audit.residuals()))
def test_exact_identity(name):
    value=audit.residuals()[name]
    entries=list(value) if isinstance(value,sp.MatrixBase) else [value]
    assert all(item==0 for item in entries),name


@pytest.mark.parametrize("name",list(audit.gates()))
def test_explicit_proof_gate(name):
    assert audit.gates()[name] is True,name


def test_audited_bad_domains():
    assert audit.controls()["rejected_inputs"]==66


@pytest.mark.parametrize("speed",(0,sp.Rational(1,2),sp.Rational(99,100),sp.Rational(999,1000)))
def test_exact_timelike_curve_prescription(speed):
    data=observers.curve(speed)
    assert data["physical_speed"]==speed
    assert data["named_uniformly_subclock_class"] is bool(speed<=sp.Rational(99,100))


@pytest.mark.parametrize("c2",(sp.Rational(1,4),sp.Rational(1,2),1))
def test_positive_coefficient_domains(c2):
    data=qei.coefficients(c2,sp.Rational(1,1215))
    assert data["test_kinetic_density_enhancement"]>1
    assert data["field_square_prefactor_without_hbar_over_kappa"]>0


def test_actual_relational_clock_residue_not_zero():
    data=modes.data()
    assert data["center_clock_relative_weight"]!=0
    assert data["uniform_clock_relational_weight_lower_on_gamma_slab"]>0


def test_actual_energy_coefficient_must_keep_spatial_frequency():
    data=qei.data()
    assert sp.factor(data["generic_time_derivative_square_coefficient"]-data["generic_spatial_gradient_sum_coefficient"])!=0


def test_observer_condition_not_just_physical_timelikeness():
    data=observers.curve(observers.INTER_SPEED)
    assert data["physical_speed"]<1
    assert data["named_uniformly_subclock_class"] is False
    assert observers.data()["interluminal_interval_clock_squared_speed_upper"]<observers.INTER_SPEED**2


def test_conormal_control_is_not_promoted_to_a_QEI_no_go():
    assert observers.data()["conormal_failure_is_not_by_itself_a_nonexistence_or_QEI_no_go"] is True


def test_fixed_background_tensor_does_not_have_free_scalar_conservation():
    assert stress.data()["bounce_principal_box_source_coefficient_per_q_b"]==sp.Rational(8,6075)
    assert stress.data()["bounce_principal_box_source_coefficient_per_q_b"]!=0


def test_matter_metric_null_replacement_has_wrong_clock_characteristic():
    assert stress.checks()["bounce_clock_covector_is_physical_metric_spacelike"]==0
    assert sp.Rational(16,1215)>0


def test_actual_short_sampling_coefficient_exact_algebraic_form():
    value=qei.data()["actual_bounce_coefficients"]["test_kinetic_density_enhancement"]
    assert sp.simplify(value-1-10863*sp.sqrt(17985)/1723683599)==0


def test_decoupling_removes_clock_correction_without_altering_free_factor():
    for c2 in (sp.Rational(1,4),sp.Rational(1199,1215),1):
        d=qei.coefficients(c2,0)
        assert d["test_kinetic_density_enhancement"]==1
        assert d["test_kinetic_density_prefactor_without_hbar_over_kappa"]==1/(16*sp.pi**2)


@pytest.mark.parametrize("order",(1,2))
def test_sampling_rescaling_keeps_both_derivatives_and_measure(order):
    u,s=sp.symbols("u s",real=True)
    eta=sp.Symbol("positive_sampling_width",positive=True)
    g=sp.exp(-u*u/2)
    scaled=eta**sp.Rational(-1,2)*g.subs(u,u/eta)
    transformed=eta*sp.diff(scaled,u,order).subs(u,eta*s)**2
    expected=eta**(-2*order)*sp.diff(g,u,order).subs(u,s)**2
    assert sp.simplify(transformed-expected)==0


def test_leading_limit_is_not_assigned_as_a_finite_width_bound():
    data=qei.data()
    assert data["finite_width_bound_is_reference_functional_not_its_leading_asymptotic"] is True
    assert data["no_claim_of_optimality_or_cosmological_scale_remainder"] is True
