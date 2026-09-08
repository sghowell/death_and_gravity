"""Independent flat-frequency, contour-boundary and regular-chart controls."""
from fractions import Fraction

import sympy as sp
from p8_vector_metric_dispersion import chart, principal, proofs, spectral, verify


def zero(values):
    for value in values.values():
        if isinstance(value, sp.MatrixBase):
            assert all(entry == 0 for entry in value)
        else:
            assert value == 0


def test_physical_pairs_and_static_contacts():
    zero(spectral.checks())


def test_adiabatic_coefficients_equal_independent_frequency_Taylor():
    assert len(spectral.adiabatic_checks()) == 16
    zero(spectral.adiabatic_checks())


def test_original_curved_principal_coefficients_are_replayed():
    zero(principal.local_checks())


def test_real_axis_inequalities_and_finite_radial_limits():
    zero(principal.checks())


def test_nonlinear_chart_chain_and_fixed_tadpole_cancellation():
    zero(chart.checks())


def test_independent_spectral_matrix_at_rational_momenta():
    for z in (Fraction(0), Fraction(1, 3), Fraction(3, 4), Fraction(1)):
        a, b = Fraction(4, 9), Fraction(28, 81)
        transverse = (b*(1-z), 2*(1-z))
        longitudinal = (b+a*z, 2*(1+z))
        expected = [[2*transverse[i]*transverse[j]+longitudinal[i]*longitudinal[j]
                     for j in range(2)] for i in range(2)]
        actual = spectral.actual(spectral.matrix()).subs({spectral.h: 1, spectral.z: sp.Rational(z)})
        assert actual == sp.ImmutableMatrix(expected)


def test_independent_radial_constant_weights():
    odd_sums = [sum((Fraction(1, 2*k-1) for k in range(1, j+1)), Fraction()) for j in (1, 2, 3)]
    assert odd_sums == [Fraction(1), Fraction(4, 3), Fraction(23, 15)]
    expected = sum((spectral.coefficients()[j]*sp.Rational(odd_sums[j])/4 for j in range(3)), sp.zeros(2))
    assert spectral.finite_asymptotic() == expected


def test_exact_moments_differentiate_to_integrand_and_domain_is_not_early_expansion():
    d = sp.Symbol("d", positive=True)
    assert sp.simplify(spectral.moment(2, d)-d*spectral.moment(1, d)+sp.Rational(1, 3)) == 0
    assert sp.simplify(spectral.moment(3, d)-d*spectral.moment(2, d)+sp.Rational(1, 5)) == 0


def test_pole_rank_one_but_finite_complement_is_nonzero():
    assert principal.pole().det() == 0
    assert sp.factor(principal.finite().det()) == -sp.Rational(50992, 59049)/principal.h**2
    Q = principal.change()
    complement = (Q.T*principal.asymptotic()*Q)[0, 0]
    assert sp.factor(complement) == sp.Rational(15616, 98415)/principal.h**2


def test_uniform_transformed_real_axis_inverse_bound_is_rational():
    data = principal.real_axis()
    expected = Fraction(1635582375, 63963136)
    assert data["uniform_interval_inverse_row_sum_upper"] == sp.Rational(expected)
    assert expected < 26
    assert data["determinant_upper"] == -sp.Rational(62464, 98415)/principal.h**2


def test_leading_log_zero_not_in_large_frequency_regime():
    assert -sp.Rational(9137, 58560) < 0
    # Exact real-axis block has a strictly signed determinant, including there.
    assert principal.real_axis()["positive_diagonal_lower"].is_positive
    assert principal.real_axis()["negative_diagonal_upper"] == -4


def test_omitting_contacts_fails_static_control_at_fixed_fixture():
    item = spectral.contact("T")
    value = item["contact"][0, 0].subs({
        spectral.m: 1, spectral.w: 2, spectral.beta: 1, spectral.beta2: 0})
    assert value == -sp.Rational(1, 4)
    assert value != 0


def test_omitting_chart_contact_before_background_cancellation_fails():
    value = 3*chart.data()["second_lapse_jet"].subs(chart.h, 1)
    assert value == -sp.Rational(3, 2)
    assert value != 0


def test_chart_principal_matrices_keep_the_lapse_fourth_derivative():
    data = chart.data()
    assert data["finite_principal"][0, 0] == -sp.Rational(7357, 6561)/chart.h**2
    assert data["pole_principal"][0, 0] == sp.Rational(12769, 13122)/chart.h**2


def test_chart_product_bound_and_one_response_are_not_same_norm_iteration():
    data = chart.norm()
    assert data["chart_C10_product_upper"] == sp.Rational(46090764897, 8)
    assert all(x < sp.Rational(1, 10**25) for x in data["one_retarded_tree_application_C10_to_C0_upper"].values())
    assert "same-norm" in data["scope"]
    assert "laplace" not in data["one_retarded_tree_application_C10_to_C0_upper"]


def test_strict_interfaces_after_cache_population():
    spectral.pair("T")
    spectral.contact("L")
    spectral.moment(1, 2)
    assert verify.controls()["rejected_inputs"] == 52


def test_all_exact_residuals_and_audit_gates():
    zero(proofs.residuals())
    assert len(proofs.residuals()) == 56
    assert len(proofs.checks()) == 14
    assert all(value is True for value in proofs.checks().values())
