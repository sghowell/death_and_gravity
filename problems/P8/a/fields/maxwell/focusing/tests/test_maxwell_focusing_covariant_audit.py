"""Independent history, shrinking-window, source and completeness controls."""

from fractions import Fraction as F

import sympy as sp
from p8a_maxwell import stress
from p8a_maxwell_focusing import controls, envelope, focusing, history


def test_both_exact_index_primitives_have_the_correct_initial_K_sign():
    t = sp.Symbol("t", real=True)
    h, u, g = (sp.Function(name)(t) for name in ("H", "u", "g"))
    ricci = 3*(sp.diff(h, t)+h*h)
    past = ricci*u*u-3*((sp.diff(u, t)-h*u)**2-sp.diff(u, t)**2)
    future = 3*sp.diff(g, t)**2+ricci*g*g-3*(sp.diff(g, t)-h*g)**2
    assert sp.expand(past-sp.diff(3*h*u*u, t)) == 0
    assert sp.expand(future-sp.diff(3*h*g*g, t)) == 0
    # u rises from0 to1; g falls from1 to0. Their boundary contributions
    # are +K and -K, respectively, not a doubled initial contraction.
    k = sp.Symbol("K", real=True)
    assert (k*1-k*0)+(k*0-k*1) == 0


def test_past_orientation_preserves_the_physical_Maxwell_EED_and_sampler_cost():
    h0, h1, h2, h3, f0, f1, f2 = sp.symbols("H Hdot Hddot Hthird f fdot fddot", real=True)
    beta = sp.Symbol("beta", real=True)
    original = stress.reference_jets(h0, h1, h2, h3, beta)["EED"]
    reversed_eed = stress.reference_jets(-h0, h1, -h2, h3, beta)["EED"]
    assert sp.expand(original-reversed_eed) == 0
    operator = f2-2*h0*f1+(sp.Rational(3, 4)*h0*h0-sp.Rational(3, 2)*h1)*f0
    reverse = operator.subs({h0: -h0, f1: -f1}, simultaneous=True)
    assert sp.expand(operator-reverse) == 0


def test_history_excess_bound_uses_a_nonnegative_remainder_not_a_Ricci_sign():
    u, up, delta_h, magnitude, tau = sp.symbols("u uprime d h tau", nonnegative=True)
    h = -magnitude/tau-delta_h
    remainder = 3*((up-h*u)**2-up**2)-6*magnitude*u*up/tau
    remainder -= 3*magnitude**2*u*u/tau**2
    expected = 6*delta_h*u*up+6*delta_h*magnitude*u*u/tau+3*delta_h**2*u*u
    assert sp.expand(remainder-expected) == 0
    x = sp.Symbol("x", real=True)
    p = 3*x*x-2*x**3
    cross = sp.integrate(p*sp.diff(p, x), (x, 0, 1))
    zeroth = sp.integrate(p*p, (x, 0, 1))
    h, ratio = sp.Rational(3, 2), sp.Rational(1, 100)
    assert 6*h*cross+3*h*h*ratio*zeroth == history.lower(h, ratio)


def test_shrinking_pole_integrals_are_finite_due_to_quadratic_endpoint_zero():
    x = sp.Symbol("x", positive=True)
    p = 3*x*x-2*x**3
    assert sp.cancel(p/x**2) == 3-2*x
    assert sp.integrate((3-2*x)**2, (x, 0, 1)) == sp.Rational(13, 3)
    assert sp.integrate((6-6*x)**2, (x, 0, 1)) == 12
    assert sp.integrate((sp.diff(p, x, 2))**2, (x, 0, 1)) == 12
    # A merely linear endpoint zero cannot be substituted in this estimate.
    assert sp.integrate(1/x**2, (x, 0, 1)) == sp.oo


def test_separate_duration_and_anomaly_scalings_match_independent_Fraction_costs():
    ratio, h = F(1, 100), F(3, 2)
    c0, c1, c2, c3 = map(F, (2, 4, 16, 96))
    q = F(3, 4)*c0*c0+F(3, 2)*c1
    v = 186*(c0**4+2*c0*c0*c1)+18*c3+90*c1*c1+90*c0*c2+108*c0*c0*c1
    past = (F(7, 2)+F(11, 5)*c0*ratio+F(2, 3)*q*ratio**2)**2/ratio**3
    past += F(13, 12600)*v*ratio
    future = (F(7, 2)*(1+2*c0)+F(25, 12)*q)**2+F(13, 1080)*v
    gain = 3*h+F(39, 35)*h*h*ratio
    data = envelope.coefficients((2, 4, 16, 96), (2, 4, 16, 96), ratio, beta_m=-1)
    assert data["past_cost"] == sp.Rational(past.numerator, past.denominator)
    assert data["future_cost"] == sp.Rational(future.numerator, future.denominator)
    assert v == 16704 and past+future == F(2638797936917, 210000)
    assert gain-F(18, 5)-(past+future)/10**8-F(1313, 3500) > F(1, 3)


def test_geometric_source_dictionary_keeps_both_histories_and_all_signs():
    data = focusing.physical_enclosure(kappa=2, hbar=3, tau=5,
        cosmological_constant=7, other_eed_lower=-11)
    assert data["separate_source_penalty"] == 29
    assert data["sigma"] == 25*29
    assert data["exact_delta"] == sp.Rational(3, 100)/sp.pi**2
    assert data["rational_delta_upper"] == sp.Rational(1, 300)
    ratio = sp.Rational(1, 100)
    source_weight = sp.Rational(13, 35)*(1+ratio)
    assert source_weight == sp.Rational(1313, 3500)
    assert focusing.calibration()["source_sigma_at_most_one"]["source_cost"] == source_weight
    assert source_weight != sp.Rational(13, 35)


def test_quintic_weak_jet_bounds_and_positive_mollifier_history_margin():
    x = sp.Symbol("x", real=True)
    p = 10*x**3-15*x**4+6*x**5
    assert sp.expand(p-controls.switch(x)) == 0
    assert sp.factor(sp.diff(p, x)) == 30*x*x*(x-1)**2
    # Maxima of p5' and |p5''|, including critical points and endpoints.
    assert sp.diff(p, x).subs(x, sp.Rational(1, 2)) == sp.Rational(15, 8)
    critical = (sp.S.Zero, sp.S.One, (3-sp.sqrt(3))/6, (3+sp.sqrt(3))/6)
    assert all(sp.simplify(sp.diff(p, x, 2).subs(x, point)**2) <= sp.Rational(100, 3)
               for point in critical)
    assert sp.expand(10*x**3-p) == 15*x**4-6*x**5
    amplitude, width = F(31, 20), F(1, 1000)
    assert amplitude*(1-10*width**3) > F(3, 2)
    assert all(cap > value for cap, value in zip((2, 4, 16, 96),
        (amplitude, amplitude*F(15, 8), amplitude*6, amplitude*60), strict=True))


def test_complete_control_has_a_positive_late_scale_factor_and_infinite_affine_length():
    t, a, momentum = sp.symbols("t a P", positive=True)
    # After the compact Hubble transition the metric is exactly static.
    timelike_dlambda_dt = 1/sp.sqrt(1+momentum**2/a**2)
    null_dlambda_dt = a/momentum
    assert sp.limit(t*timelike_dlambda_dt, t, sp.oo) == sp.oo
    assert sp.limit(t*null_dlambda_dt, t, sp.oo) == sp.oo
    control = controls.calibration()
    assert control["future_timelike_and_null_complete_geometry"] is True
    assert control["actual_allowed_Maxwell_SEE_solution"] is False
    assert control["geometric_assumptions_alone_imply_the_conclusion"] is False


def test_constant_past_reference_requires_an_excluded_large_source_penalty():
    amplitude, delta = sp.Rational(31, 20), sp.Rational(1, 10**8)
    # R_UU*tau²=3A²; kappa*E_conf*tau²=-(31/60)delta A^4.
    required_penalty = 3*amplitude**2-sp.Rational(31, 60)*delta*amplitude**4
    data = controls.calibration()
    assert data["reference_vacuum_required_sigma_at_Lambda_zero"] == required_penalty > 1
    assert data["initial_pointwise_SEC_holds"] is False


def test_zero_history_or_large_extra_source_cannot_be_hidden_in_positive_test_flag():
    base = {
        "ratio": sp.Rational(1, 100), "contraction": sp.Rational(3, 2),
        "past_caps": (2, 4, 16, 96), "future_caps": (2, 4, 16, 96), "beta_m": 1,
        "delta": sp.Rational(1, 10**8),
    }
    assert focusing.theorem_constants(**base, sigma=1)["sufficient_focusing_test"] is True
    failed = focusing.theorem_constants(**base, sigma=10)
    assert failed["sufficient_focusing_test"] is False
    assert failed["actual_SEE_or_geometric_hypotheses_verified"] is False
