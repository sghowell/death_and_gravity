"""Independent actual-source, weighted-loss, two-frequency and index audit."""

from fractions import Fraction as Q

import pytest
import sympy as sp
from p8a_extended_qsei import focusing, reference, sampling


def exact(value):
    return sp.Rational(value.numerator, value.denominator)


def test_actual_quantum_eed_uses_full_see_and_subtracts_both_sources():
    x = sp.Symbol("x", real=True)
    a = sp.Function("a", positive=True)(x)
    c = sp.Function("c")(x)
    delta, time_scale, hbar, coupling = sp.symbols("delta T0 hbar kappa", positive=True)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    physical_hubble = h/(time_scale*a)
    proper = lambda value: sp.diff(value, x)/(time_scale*a)
    total = -3*(proper(physical_hubble)+physical_hubble**2)/coupling
    radiation = 3/(coupling*time_scale**2*a**4)
    conversion = hbar/(sp.pi**2*time_scale**4)
    external_rho = conversion*c/a**4
    external_p = conversion*(c-sp.diff(c, x)/h)/(3*a**4)
    external_eed = (external_rho+3*external_p)/2
    quantum = (total-radiation-external_eed).subs(
        coupling, 2880*sp.pi**2*delta*time_scale**2/hbar)/conversion
    code_free = reference.energy_dimensionless(a, h, u, delta)
    assert sp.simplify(quantum-code_free+(c-sp.diff(c, x)/(2*h))/a**4) == 0
    assert sp.simplify(quantum.subs({c: 0, sp.diff(c, x): 0})-code_free) == 0
    # During preparation a constant nonzero weighted source is still not
    # part of the quantum state. Its subtraction cannot be silently omitted.
    assert sp.simplify((quantum-code_free).subs({c: 1, sp.diff(c, x): 0})+a**-4) == 0


def test_signed_reference_loss_keeps_weight_and_correct_denominator():
    delta, sigma, distance, weight = Q(1, 10**14), Q(2*10**6), Q(9, 10**8), Q(9)
    b = Q(3, 10**12)
    # Integrating X, the damped h equation, then log(a) gives gains
    # sigma^-1,sigma^-2,sigma^-3. The a² logarithmic derivative is <=18.
    numerator_gain = 9/sigma+9/sigma**2+18*(b+Q(1, 4))/sigma**3
    loss = weight*distance*numerator_gain
    negative = loss/(960*delta*2**4)
    actual = reference.calibration()
    assert actual["numerator_weighted_pair"] == exact(numerator_gain)
    assert actual["numerator_loss_upper"] == exact(loss)
    assert actual["negative_EED_magnitude_upper"] == exact(negative)
    assert Q(1, 50) < negative < Q(1, 40)
    assert actual["actual_weighted_distance"] < exact(distance)
    assert actual["reference_EED_lower_dimensionless"] == -sp.Rational(1, 40)
    # Both tempting losses would be downward: forgetting exp(sigma L),
    # or dividing the negative numerator by the upper bound on a⁴.
    assert negative/weight < negative
    assert loss/(960*delta*3**4) < negative
    assert loss > delta/27  # The old positive-credit argument does not transfer.
    assert actual["positive_reference_credit_asserted"] is False


def test_positive_old_numerator_is_not_an_actual_reference_sign_claim():
    y, delta = sp.symbols("y delta", positive=True)
    f = 1-delta/y**4
    a2 = y**2*sp.sqrt(f)
    h2 = 1/(y**2*f**sp.Rational(3, 2))
    u = 2*delta/(y**6*f**sp.Rational(3, 2))
    numerator = sp.simplify(a2*(u+h2)-1)
    assert sp.simplify(numerator-3*delta/(y**4-delta)) == 0
    assert numerator.subs({y: 3, delta: sp.Rational(1, 10**14)}) > 0
    assert reference.lower_bound(time_scale=2, hbar=3) == -sp.Rational(3, 640)/sp.pi**2


def test_full_fraction_recalibration_uses_long_span_and_signed_credit():
    ell, history, m = Q(1, 10**6), Q(3), Q(1, 10**5)
    b, h0, h1, h2 = m*history, Q(1, 2), Q(1, 3), Q(1, 3)
    proper_h, proper_hdot = Q(1, 4), Q(1, 7)
    f0 = ell**2/4
    f1 = (ell+Q(3, 2)*proper_h*ell**2)/2
    f2 = 1+2*proper_h*ell+(Q(3, 4)*proper_h**2+Q(3, 2)*proper_hdot)*ell**2
    hf1 = proper_h*ell+(proper_hdot+proper_h**2/2)*ell**2
    q0, q1 = 2*b, 2*m+2*history*b*b
    a0, a1, a2 = 2*b*history, q0, q1
    c0 = (1+h0)*a0/2
    c1 = ((1+h0)*a1+h1*a0)/2
    c2 = ((1+h0)*a2+2*h1*a1+h2*a0)/2
    backward = (1+h0)/4
    p0, p1 = backward*q0, backward*q1+h1*q0/4
    remainder = m*(2*history*b)+b*(b+history*m+history**2*b*b)
    infrared = (1+h0)*b*history**2+2*b*history
    roots = {
        "forward": Q(2, 1000)*(c0*f2+2*c1*f1+c2*f0),
        "local": Q(2, 1000)*(p0*f1+p1*f0),
        "real_history_Parseval": 8*backward*m*f0,
        "history_remainder": 2*backward*history*remainder*f0,
        "infrared": 2*infrared*f0,
    }
    root = f2+Q(3, 2)*hf1+sum(roots.values())
    penalty = Q(2, 5)*ell**4
    d = sampling.calibration()
    assert d["sampler_norms"] == [exact(v) for v in (f0, f1, f2)]
    assert d["spectral_error_roots"] == {name: exact(value) for name, value in roots.items()}
    assert d["total_spectral_root"] == exact(root)
    assert d["signed_reference_penalty"] == exact(penalty) > 0
    assert d["derived_absolute_coefficient"] == exact(root**2+penalty)
    assert root < 1+Q(2, 10**6)
    assert root**2+penalty < 2
    mode = sampling.scattering_calibration()
    assert mode["history"] == 3
    assert mode["potential_cap"] == sp.Rational(3, 10**5)
    assert mode["derivative_cap"] == sp.Rational(1, 10**5)
    assert mode["hubble_caps"] == [sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 3)]


def test_sampler_clock_measure_and_reference_penalty_are_physical():
    s = sp.Symbol("s", real=True)
    a, psi = sp.Function("a", positive=True)(s), sp.Function("psi")(s)
    conformal = lambda value: a*sp.diff(value, s)
    f = psi/a**sp.Rational(3, 2)
    hubble = sp.diff(a, s)/a
    proper_hdot = sp.diff(hubble, s)
    assert sp.simplify(conformal(conformal(f))/sp.sqrt(a)-sp.diff(psi, s, 2)
                       +2*hubble*sp.diff(psi, s)
                       -(3*hubble**2/4-3*proper_hdot/2)*psi) == 0
    assert sp.simplify(conformal(sp.diff(a, s)*f)/sp.sqrt(a)-hubble*sp.diff(psi, s)
                       -(proper_hdot-hubble**2/2)*psi) == 0
    t0, ell, d2, hbar = sp.symbols("T0 ell D2 hbar", positive=True)
    negative_reference_magnitude = hbar/(40*sp.pi**2*t0**4)*(t0*ell**4*d2**2)
    normalized_penalty = hbar/(16*sp.pi**2)*(sp.Rational(2, 5)*ell**4)*(d2**2/t0**3)
    assert sp.simplify(negative_reference_magnitude-normalized_penalty) == 0


def test_two_frequency_history_estimate_uses_full_then_real_parseval():
    k, alpha = sp.symbols("k alpha", positive=True)
    forward_alpha = sp.integrate((alpha+k)**-4, (alpha, 0, sp.oo))
    assert sp.integrate(k*forward_alpha, (k, 1, sp.oo)) == sp.Rational(1, 3)
    real_half = sp.integrate(1/(1+4*k*k), (k, 0, sp.oo))
    complex_half = sp.integrate(1/(1+(2*k-1)**2), (k, 0, sp.oo))
    assert real_half == sp.pi/4
    assert complex_half == 3*sp.pi/8 > real_half
    # The first Fourier norm has a complex sampler-history product and
    # must use full Parseval; only the real u' history admits the half norm.
    amplitude, history, derivative, sampler_norm = sp.symbols("A T M F0", positive=True)
    two_parsevals = 2*sp.pi*amplitude**2*sampler_norm**2*(sp.pi*history*derivative**2/2)
    assert sp.expand(two_parsevals-sp.pi**2*amplitude**2*history*derivative**2*sampler_norm**2) == 0
    assert sp.integrate(k, (k, 0, 1)) == sp.Rational(1, 2)


def test_support_domain_keeps_original_source_switch_not_new_half():
    old, new = sp.Rational(1, 10**10), sp.Rational(1, 10**6)
    d = sampling.calibration()
    assert d["source_off_from"] == old/2
    assert d["source_free_conformal_width"] == new-old/2
    assert sampling.validate_support_interval(old, new/4) is True
    assert old < new/4 < new/2  # Valid samplers need not lie in the new half.
    for left, right in [(old/2, new/4), (old, new), (0, old/4)]:
        with pytest.raises(ValueError, match="source-free"):
            sampling.validate_support_interval(left, right)


def test_actual_index_curvature_and_q2_ratio_do_not_imply_focusing():
    x = sp.Symbol("x", real=True)
    a = sp.Function("a", positive=True)(x)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    ricci = 3*(sp.diff(h/a, x)/a+(h/a)**2)
    assert sp.simplify(ricci+3*(u+h*h)/a**2) == 0
    duration = sp.Symbol("D", positive=True)
    r = sp.Symbol("r", nonnegative=True)
    assert sp.integrate(duration-r, (r, 0, duration)) == duration**2/2
    lower = 3/duration-duration/8
    assert sp.diff(lower, duration) == -3/duration**2-sp.Rational(1, 8)
    span, delta = sp.Rational(3, 10**6), sp.Rational(1, 10**14)
    assert lower.subs(duration, span) > sp.Rational(3, 4)
    d = focusing.calibration()
    assert d["index_times_T0_lower"] == lower.subs(duration, span)
    assert d["normal_jacobian_lower"] == sp.Rational(8, 27)
    assert d["Q2_over_available_duration_squared_lower"] == 360*delta/span**2 == sp.Rational(2, 5)
    # Q2 is for this sufficient coefficient, not a lower bound on optimal QEIs.
    assert d["Q2_over_T0_squared"] == 180*delta*sampling.calibration()["rounded_absolute_coefficient"]
