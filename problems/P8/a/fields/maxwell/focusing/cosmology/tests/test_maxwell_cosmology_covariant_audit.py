"""Independent clock, full-history, prescription and source-budget audit."""

from fractions import Fraction as F

import sympy as sp
from p8a_maxwell_cosmology import calibration, controls, dictionary, geometry


def test_reversed_physical_power_law_clock_and_all_four_hubble_jets():
    p, tau = sp.symbols("p tau", positive=True)
    s, x = sp.symbols("s x", real=True)
    age = p*tau/2
    scale = ((age-s)/age)**p
    hubble = sp.simplify(sp.diff(scale, s)/scale)
    assert sp.simplify(hubble+p/(age-s)) == 0
    for j, expected in enumerate(geometry.reference_jets(p, x)):
        actual = (tau**(j+1)*sp.diff(hubble, s, j)).subs(s, tau*x)
        assert sp.simplify(actual-expected) == 0
    assert sp.simplify(-hubble.subs(s, 0)*tau) == 2


def test_all_power_and_history_bounds_have_nonnegative_polynomial_numerators():
    offset, past = sp.symbols("offset past", nonnegative=True)
    p = sp.Rational(1, 2)+offset
    denominator = p/2+past
    caps = (2, 8, 64, 768)
    for j, cap in enumerate(caps):
        numerator = sp.expand(cap*denominator**(j+1)-sp.factorial(j)*p)
        assert all(c >= 0 for c in sp.Poly(numerator, offset, past).coeffs())
    # The history minimum also holds uniformly, not just at sampled p.
    minimum = sp.Rational(25, 13)
    residual = sp.expand(p-minimum*(p/2+sp.Rational(1, 100)))
    assert residual == offset/26
    assert geometry.history_point(sp.Rational(1, 2), -sp.Rational(1, 100))[0] == -minimum


def test_robust_tube_has_strict_history_and_derivative_slack():
    data = geometry.neighborhood()
    actual_caps = [sp.Integer(2)+sp.Rational(1, 100), sp.Rational(17, 2), 67, 784]
    expanded = [sp.Rational(21, 10), 9, 70, 800]
    assert [b-a for a, b in zip(actual_caps, expanded, strict=True)] == data["strict_cap_margins"]
    assert data["strict_contraction_margin"] == sp.Rational(17, 1300)
    assert data["observed_H0_times_tau_range"] == [sp.Rational(199, 100), sp.Rational(201, 100)]
    assert geometry.neighborhood(anchored=True)["observed_H0_times_tau_range"] == [2, 2]
    assert data["future_relative_caps_verified_by_tube"] is False
    assert data["actual_quantum_SEE_solution_asserted"] is False


def test_reference_age_and_normalized_scale_do_not_fix_a_quantum_age():
    data = dictionary.proper_scales(100)
    assert data["actual_short_history_duration"] == 1
    assert data["reference_power_law_age_range"] == [25, sp.Rational(100, 3)]
    assert data["actual_observed_H0_range"] == [sp.Rational(199, 10000), sp.Rational(201, 10000)]
    assert data["reference_power_law_age_is_actual_quantum_age"] is False
    assert geometry.neighborhood()["log_scale_factor_error_upper"] == sp.Rational(1, 10000)


def test_affine_cost_is_reconstructed_from_independent_fraction_integrals():
    c0, c1, c2, c3 = F(21, 10), F(9), F(70), F(800)
    r, q = F(1, 100), F(3, 4)*c0*c0+F(3, 2)*c1
    v0 = 186*(c0**4+2*c0*c0*c1)
    vb = 18*c3+90*c1*c1+90*c0*c2+108*c0*c0*c1
    root_p = F(7, 2)+F(11, 5)*c0*r+F(2, 3)*q*r*r
    f0, f1, f2, f3 = map(F, (4, 128, 16384, 1048576))
    future_q = F(3, 4)*f0*f0+F(3, 2)*f1
    future_v0 = 186*(f0**4+2*f0*f0*f1)
    future_vb = 18*f3+90*f1*f1+90*f0*f2+108*f0*f0*f1
    root_f = F(7, 2)*(1+2*f0)+F(25, 12)*future_q
    past_weight, future_weight = F(13, 35)*r/360, F(13, 3)/360
    independent_c0 = root_p*root_p/r**3+root_f*root_f+past_weight*v0+future_weight*future_v0
    independent_cb = past_weight*vb+future_weight*future_vb
    data = calibration.affine_cost()
    assert independent_c0 == F(str(data["C0"])) == F(1440182116554809, 112500000)
    assert independent_cb == F(str(data["Cbeta"])) == F(3345309847373, 10500000)
    assert data["past_caps"] != data["future_caps"]
    assert independent_c0 < 13000000 and independent_cb < 13000000


def test_joint_finite_prescription_budget_cannot_be_replaced_by_delta_alone():
    for beta in (0, 1, -1, 1000000, -1000000):
        actual = calibration.cost_at_beta(beta)
        assert actual < 13000000*(1+abs(beta))
        weighted_delta = sp.Rational(1, 10**8)
        delta = weighted_delta/(1+abs(beta))
        assert delta*actual < sp.Rational(13, 100)
    # Holding delta fixed while discarding its beta weight fails badly.
    assert calibration.cost_at_beta(1000000)/10**8 > 1


def test_complete_source_and_history_margin_is_independent_fraction_arithmetic():
    h, r = F(19, 10), F(1, 100)
    gain = 3*h+F(39, 35)*h*h*r
    weight = F(13, 35)*(1+r)
    lower = gain-F(18, 5)-F(13, 100)-5*weight
    assert lower == F(47079, 350000) > F(1, 8)
    data = calibration.theorem_gate(sp.Rational(1, 10**8), 5)
    assert F(str(data["strict_focusing_margin_lower"])) == lower
    assert F(str(data["margin_above_one_eighth"])) == F(3329, 350000)
    assert data["future_relative_caps_or_actual_SEE_verified"] is False


def test_signed_sources_and_actual_versus_reference_lambda_fraction():
    for lam, lower, sigma in ((1, 2, 1), (1, -2, 5), (-1, 2, 0), (-1, -2, 4)):
        data = dictionary.source_budget(kappa=2, hbar=0, tau=1, beta_m=-7,
                                        cosmological_constant=lam, other_eed_lower=lower)
        assert data["sigma"] == sigma
        assert data["actual_field_state_source_or_metric_verified"] is False
    reference_fraction = sp.Rational(5, 12)
    actual_fraction = 5/(3*sp.Rational(199, 100)**2)
    assert actual_fraction == sp.Rational(50000, 118803) > reference_fraction
    assert dictionary.calibration_data()["Lambda_positive_fraction_of_actual_3H0_squared_upper_unanchored"] == actual_fraction


def test_planck_seconds_conversion_retains_both_powers_of_c():
    newton, hbar, c, seconds = sp.symbols("G hbar c seconds", positive=True)
    lp2 = newton*hbar/c**3
    tp2 = newton*hbar/c**5
    assert sp.simplify(lp2/(sp.pi*(c*seconds)**2)-tp2/(sp.pi*seconds**2)) == 0
    assert sp.simplify(lp2/(sp.pi*seconds**2)-tp2/(sp.pi*seconds**2)) != 0


def test_observer_jets_alone_do_not_certify_the_full_contraction_history():
    x = sp.Symbol("x", real=True)
    reference = -sp.Rational(1, 2)/(sp.Rational(1, 4)-x)
    bad_history = reference+10**9*x**4
    assert all(sp.diff(bad_history-reference, x, j).subs(x, 0) == 0 for j in range(4))
    assert bad_history.subs(x, -sp.Rational(1, 100)) > 0
    assert reference.subs(x, -sp.Rational(1, 100)) < -sp.Rational(19, 10)


def test_old_tight_future_caps_already_fail_without_the_quantum_inequality():
    x = F(1, 20)
    upper = -F(199, 100)-F(11, 2)*x-F(33, 2)*x*x+F(800)/(1-x)**4*x**3/6
    required = -F(21, 10)/(1-x)
    assert upper < required
    assert required-upper == F(4707907, 62554080)
    data = controls.rejected_tight_future()
    assert F(str(data["strict_incompatibility_gap"])) == required-upper
    assert data["QSEI_used_in_this_exclusion"] is False
    assert data["rejected_as_nontrivial_cosmological_calibration"] is True


def test_compact_geometric_extension_is_c3_and_its_mollification_fits_the_tube():
    x, p = sp.symbols("x p", real=True, positive=True)
    cutoff = 35*(8*x)**4-84*(8*x)**5+70*(8*x)**6-20*(8*x)**7
    past = -p/(p/2-x)
    future = past*(1-cutoff)
    for j in range(4):
        assert sp.simplify(sp.diff(future-past, x, j).subs(x, 0)) == 0
        assert sp.simplify(sp.diff(future, x, j).subs(x, sp.Rational(1, 8))) == 0
    # Direct binomial Leibniz bounds, separate from the primary function.
    base = (4, 32, 512, 12288, 393216)
    cut = (1, sp.Rational(35, 2), 1680, 107520, 216760320)
    raw = [sum(sp.binomial(j, k)*base[j-k]*cut[k] for k in range(j+1)) for j in range(5)]
    assert raw == [4, 102, 8352, 630528, 887218176]
    assert controls.raw_global_jet_caps() == raw
    data = controls.complete_geometry()
    assert data["past_C3_error_upper"] == [value/10**12 for value in raw[1:]]
    assert all(value > 0 for value in data["strict_past_tube_margins"])
    assert all(value >= 0 for value in data["future_cap_margins_from_global_bounds"])
    assert data["geometric_assumptions_alone_force_incompleteness"] is False
    assert data["actual_quantum_SEE_solution_asserted"] is False


def test_complete_extension_has_infinite_future_affine_lengths():
    time, amin, momentum = sp.symbols("time amin momentum", positive=True)
    assert sp.limit(time/sp.sqrt(1+momentum**2/amin**2), time, sp.oo) == sp.oo
    assert sp.limit(amin*time/momentum, time, sp.oo) == sp.oo
    data = controls.complete_geometry()
    assert data["future_static_by_time_over_tau"] == sp.Rational(1, 8)+sp.Rational(1, 10**12)
    assert data["future_log_scale_factor_loss_upper"] == sp.Rational(1, 2)+sp.Rational(4, 10**12)
    assert data["covers_every_reference_power_in_interval"] is True
    assert data["future_timelike_and_null_complete"] is True
