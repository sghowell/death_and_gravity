from fractions import Fraction

import pytest
import sympy as sp
from p8_m1_control import bounds, model, oscillator
from p8_m1_tree import majorant as m


def polynomial(series, variable):
    return sum(sp.Rational(value.numerator, value.denominator)*variable**n
               for n, value in enumerate(series.coefficients))


def test_positive_series_exact_polynomial_convolution():
    e = sp.Symbol("e")
    left = m.Series((1, 2, 3, 4, 5))
    right = m.Series((2, Fraction(1, 3), 4, 0, 6))
    actual = sp.expand(polynomial(left, e)*polynomial(right, e))
    for n, value in enumerate((left*right).coefficients):
        assert actual.coeff(e, n) == value
    assert (left+right).coefficients == tuple(a+b for a, b in zip(left.coefficients, right.coefficients))
    assert (left**3).coefficients[4] == sp.expand(polynomial(left, e)**3).coeff(e, 4)


@pytest.mark.parametrize("exponent", (sp.Rational(1, 2), sp.Rational(-1, 2),
                                     sp.Rational(-1), sp.Rational(3, 2)))
def test_binomial_majorant_against_independent_signed_polynomial(exponent):
    e = sp.Symbol("e")
    signed = e-2*e**2+3*e**3-4*e**4
    actual = sp.series((1+signed)**exponent, e, 0, 5).removeO().expand()
    bounded = m.Series((0, 1, 2, 3, 4)).analytic(exponent)
    assert all(abs(actual.coeff(e, n)) <= bound for n, bound in enumerate(bounded.coefficients))


def test_positive_series_bounds_labelled_nilpotent_products():
    left = {0: 2, 1: -3, 2: 5, 3: -7, 8: 2}
    right = {1: 1, 4: -2, 6: 5, 8: 3, 15: -4}
    product = {}
    for i, a in left.items():
        for j, b in right.items():
            if not i & j:
                product[i | j] = product.get(i | j, 0)+a*b
    norm = lambda data: m.Series(tuple(sum(abs(value) for mask, value in data.items()
                                           if mask.bit_count() == n) for n in range(5)))
    actual, bounded = norm(product), norm(left)*norm(right)
    assert all(a <= b for a, b in zip(actual.coefficients, bounded.coefficients))


def test_invalid_positive_series_operations_fail_closed():
    with pytest.raises(ValueError, match="nonnegative"):
        m.Series.at(1, -1)
    with pytest.raises(ValueError, match="five"):
        m.Series((0, 1))
    with pytest.raises(ValueError, match="zero constant"):
        m.Series.at(0, 1).analytic(sp.Rational(1, 2))
    with pytest.raises(ValueError, match="integer"):
        m.Series.at(1, 1)**-1


def test_full_metric_inverse_and_volume_are_dominated():
    e = sp.Symbol("e")
    matrix = sp.Matrix([[2, -1, 1], [-1, -2, 2], [1, 2, 1]])
    metric = sp.eye(3)+e*matrix
    inverse = sum(((-e)**n*matrix**n for n in range(5)), sp.zeros(3))
    volume = sp.series(metric.det()**sp.Rational(1, 2), e, 0, 5).removeO().expand()
    inverse_volume = sp.series(metric.det()**sp.Rational(-1, 2), e, 0, 5).removeO().expand()
    data = m.geometry_majorants(m.Series.at(1, 1), Fraction(4))
    for n in range(5):
        assert all(abs(entry.expand().coeff(e, n)) <= data["inverse"].coefficients[n] for entry in inverse)
        assert abs(volume.coeff(e, n)) <= data["volume"].coefficients[n]
        assert abs(inverse_volume.coeff(e, n)) <= data["inverse_volume"].coefficients[n]
    assert data["volume"].coefficients[0] == data["inverse_volume"].coefficients[0] == 1


@pytest.mark.parametrize("wave", ((1, 2, 3), (2, -5, 1), (0, 0, 7), (-3, -2, 5)))
def test_York_inverse_entrywise_bound(wave):
    k = sp.Matrix(wave)
    q = k.dot(k)
    inverse = (sp.eye(3)-k*k.T/(4*q))/q
    operator = -q*sp.eye(3)-k*k.T/3
    assert operator*inverse == -sp.eye(3)
    assert all(sum(abs(inverse[i, j]) for j in range(3)) <= sp.Rational(2)/q for i in range(3))
    source = sp.Matrix([2, -3, 5])
    W = inverse*source
    LW = sp.Matrix(3, 3, lambda i, j: k[i]*W[j]+k[j]*W[i]
                   -(sp.Rational(2, 3)*k.dot(W) if i == j else 0))
    K = sum(abs(value) for value in k)
    assert all(abs(entry) <= 8*K*5/q for entry in LW)


def test_exact_boundary_bound_and_required_finite_upper_band():
    report = m.canonical_boundary_bounds()
    expected = Fraction(42828424431, 2000000)
    assert Fraction(report["unitary"]["boundary_operator_bound_over_q"]) == expected
    assert 0 < Fraction(report["gamma"]["boundary_operator_bound_over_q"]) < expected
    with pytest.raises(ValueError, match="Unproved"):
        bounds.coefficient_bound(oscillator.derive("unitary")["momentum_boundary11"], "unitary")
    assert model.no_high_frequency_pole(model.z*oscillator.derive("unitary")["momentum_boundary11"])


def test_safe_band_free_column_and_terminal_seed_bounds():
    report = m.domain_checks()
    assert all(Fraction(value) > 0 for value in report["margins"].values())
    assert report["S5_7_q_threshold"] == str(10**20)
    assert report["common_terminal_phase_seed"] == str(10**34)
    seed = Fraction(report["common_terminal_phase_seed"])
    assert all(0 < Fraction(value) <= seed for value in report["terminal_phase_seed_bounds"].values())
    assert report["free_column_coordinate_norm_upper"] == "1"
    assert report["free_column_momentum_norm_upper"] == str(m.MODE_P)


@pytest.mark.parametrize("arguments", ({"lower_k": 10**10}, {"q_max": 10**24},
                                       {"derivative_k": 10**12}, {"mode_p": 1},
                                       {"lower_k": 10**13, "upper_k": 10**12}))
def test_smaller_or_inconsistent_domain_is_not_certified(arguments):
    with pytest.raises(ValueError, match="Unsafe|ordered"):
        m.domain_checks(**arguments)


def test_matter_terms_and_nonlinear_York_orders_are_not_omitted():
    full = m.physical_series(2, 5, Fraction(1, 9))
    vacuum = m.physical_series(2, 5, Fraction(1, 9), include_matter_source=False)
    linear = m.physical_series(2, 5, Fraction(1, 9), york_order=1)
    for degree in (1, 2, 3):
        assert full["momentum"].coefficients[degree] > vacuum["momentum"].coefficients[degree]
    for degree in (2, 3):
        assert full["momentum"].coefficients[degree] > linear["momentum"].coefficients[degree]
    assert full["eta"].coefficients[0] == 0
    for key in ("shear2", "matter_gradient"):
        assert full[key].coefficients[0] == full[key].coefficients[1] == 0
        assert full[key].coefficients[2] > 0
    assert full["matter_density"].coefficients[0] == Fraction(1, 5)
    assert len(full["York_correction_majorants"]) == 3


def test_independent_matter_subalgebra_omission_controls():
    for residual in m.negative_controls().values():
        assert sp.expand(residual) != 0
    e, zeta, pm, matter_l = sp.symbols("e zeta P l", real=True)
    density = matter_l+e*pm+3*matter_l*e*zeta
    direct_eta = sp.series(density*(1+2*e*zeta)**sp.Rational(-3, 2)-matter_l, e, 0, 3).removeO().expand()
    assert direct_eta.coeff(e, 1) == pm
    assert direct_eta.coeff(e, 2) == -3*pm*zeta-sp.Rational(3, 2)*matter_l*zeta**2
    assert direct_eta.coeff(e, 2) == m.negative_controls()["matter_eta_nonlinear_second_order"]
    assert m.EXCHANGE_FACTOR == 24
    assert m.EXCHANGE_FACTOR > 18


def test_homogeneous_scale_restoration():
    assert all(sp.cancel(value) == 0 for value in m.scale_checks().values())


def test_all_invariant_records_and_tree_block_inequalities():
    result = m.build()
    records = result["stationary_coefficient_majorants"]
    assert len(records) == 58
    assert all(-2 <= row["unit_weight"] <= 6 for row in records)
    assert any(int(row["powers"].split(",")[2]) > 0 for row in records)
    assert any(int(row["powers"].split(",")[4]) > 0 for row in records)
    assert all(Fraction(row["fixed_unit_bound"]) == 64*Fraction(row["compact_bound"]) for row in records)
    B3, B4 = (Fraction(result[key]) for key in ("cubic_kernel_bound", "quartic_Hamiltonian_kernel_bound"))
    C3, C4 = (Fraction(result[key]) for key in
              ("cubic_transition_block_bound_numerator", "quartic_connected_tree_block_bound_numerator"))
    schur = Fraction(result["window"]["Schur_measure_majorant"])
    assert C3 == schur*B3 and C4 == schur*(B4+24*B3**2)
    scale = int(result["sufficient_M_tau"])
    assert result["sufficient_M_tau_power10"] == 324
    assert scale == 10**result["sufficient_M_tau_power10"]
    assert C3/scale <= Fraction(1, 1000) and C4/scale**2 <= Fraction(1, 1000)
    assert scale//10 < 1000*C3 or (scale//10)**2 < 1000*C4
    assert result["window"]["propagating_columns"] == 4
    assert result["not_a_necessary_scale_or_all_orders_cutoff"]
