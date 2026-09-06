"""Independent weighted full-map, unchanged-source and physical-barrier audits."""

from fractions import Fraction as Q

import sympy as sp


def test_causal_weighted_primitive_and_damped_metric_gain():
    t, sigma, damping = sp.symbols("t sigma damping", positive=True)
    s = sp.Symbol("s", nonnegative=True)
    primitive = sp.integrate(sp.exp(sigma*s), (s, 0, t))*sp.exp(-sigma*t)
    assert sp.simplify(primitive-(1-sp.exp(-sigma*t))/sigma) == 0
    damped = sp.integrate(sp.exp(-damping*(t-s))*sp.exp(sigma*s), (s, 0, t))
    assert sp.simplify(sp.exp(-sigma*t)*damped
                       -(1-sp.exp(-(sigma+damping)*t))/(sigma+damping)) == 0
    # The actual h-difference has positive time-dependent damping, so
    # comparison bounds it by the undamped primitive, not by exp(+t).
    h1, h2, u1, u2 = (sp.Function(name)(t) for name in ("h1", "h2", "u1", "u2"))
    derivative = sp.diff(h1-h2, t).subs({sp.diff(h1, t): -u1-h1**2,
                                       sp.diff(h2, t): -u2-h2**2})
    assert sp.expand(derivative+(h1+h2)*(h1-h2)+(u1-u2)) == 0


def test_logarithmic_metric_lipschitz_factors_keep_the_physical_caps():
    loga = sp.Symbol("loga", real=True)
    assert sp.diff(sp.exp(2*loga), loga).subs(loga, sp.log(3)) == 18
    assert -sp.diff(sp.exp(-2*loga), loga).subs(loga, sp.log(2)) == sp.Rational(1, 2)
    finite = -sp.Rational(19, 60)-(loga-sp.log(2))/2
    assert sp.diff(finite, loga) == -sp.Rational(1, 2)


def test_same_source_integration_by_parts_does_not_require_cutoff_jet_bounds():
    t = sp.Symbol("t", real=True)
    c, h, u = (sp.Function(name)(t) for name in ("c", "h", "u"))
    primitive_derivative = sp.diff(c/h, t)-c*(1+u/h**2)
    assert sp.simplify(primitive_derivative.subs(sp.diff(h, t), -u-h**2)
                       -sp.diff(c, t)/h) == 0
    # c is fixed on both metrics; the two nonzero variation contributions
    # are its endpoint inverse-h term and its history u/h² term.
    alpha, source, bound = sp.symbols("alpha C B", positive=True)
    endpoint = 8*source*9*alpha**2
    history = 8*source*alpha*(9*alpha+54*bound*alpha**2)
    assert sp.expand(endpoint+history-144*source*alpha**2-432*source*bound*alpha**3) == 0


def test_stiff_coefficient_is_resummed_but_rolling_mismatch_is_not_dropped():
    alpha, delta, b, ea = sp.symbols("alpha delta B Ea", positive=True)
    full = (9*alpha**2+18*b*alpha**4)/(60*delta)
    retained = (ea*alpha**2+18*b*alpha**4)/(60*delta)
    assert sp.factor(full-retained-(9-ea)*alpha**2/(60*delta)) == 0
    assert sp.diff(retained, ea) == alpha**2/(60*delta)
    # The source-free c/s² extraction does not change the fixed finite d.
    t, af = sp.symbols("t af", positive=True)
    a = sp.Function("a", positive=True)(t)
    d = -sp.Rational(19, 60)-sp.log(a/2)/2
    assert sp.simplify(sp.diff(d, t)+sp.diff(a, t)/(2*a)) == 0
    assert sp.simplify((d-d.subs(a, af))+sp.log(a/af)/2) == 0


def test_causal_prefix_bound_is_monotone_and_vanishes_at_zero():
    length, total = sp.symbols("ell T", positive=True)
    coefficient = length**2*(sp.Rational(5, 4)+sp.log(total/length)/2)
    derivative = sp.diff(coefficient, length)
    assert sp.simplify(derivative-length*(2+sp.log(total/length))) == 0
    assert sp.limit(coefficient, length, 0, dir="+") == 0
    # For ell<=L<=T this derivative is positive. Prefix input sup<=exp(sigma*t)
    # times the weighted norm; the output weight cancels that SAME factor.
    sigma, t = sp.symbols("sigma t", positive=True)
    assert sp.exp(-sigma*t)*sp.exp(sigma*t) == 1


def test_independent_fraction_full_map_and_physical_ball_margins():
    length, original, sigma = Q(1, 10**6), Q(1, 10**10), Q(2*10**6)
    alpha, delta, radius, exponential = 1/sigma, Q(1, 10**14), Q(2, 10**7), Q(9)
    source, b, v, qbar, pmax, qmax = Q(12, 10**9), Q(3, 10**12), Q(1, 10**10), Q(1, 10**8), Q(1, 10**5), Q(2, 10**8)
    mismatch, finite = delta+9*length, delta+length/4
    curvature = (b/2+Q(1, 120))*alpha**2+(b/30+Q(1, 60))*alpha**3
    source_pair = 144*source*alpha**2+432*source*b*alpha**3
    full_einstein = (9*alpha**2+18*b*alpha**4)/(60*delta)
    p_pair = full_einstein+curvature+source_pair
    q_pair = alpha*p_pair/4+pmax*alpha**4/2
    auxiliary = 9*q_pair+18*qmax*(alpha**2+alpha**3)
    remaining_einstein = (mismatch*alpha**2+18*b*alpha**4)/(60*delta)
    local = finite+(v+exponential*radius)*alpha**3/2+alpha/4+b*alpha**2/2
    m, total = Q(1, 10**5), Q(3)
    mode = m*total*length**2*(Q(5, 4)+10)+36*m*m*total**5*length
    rhs_pair = remaining_einstein+curvature+source_pair+auxiliary+local+mode
    inverse = Q(240, 769)
    center_p = 8*source*(3+length*(1+9*b))
    center = inverse*center_p*(1+9*length/4)
    contraction = inverse*rhs_pair
    fixed_point = center/(1-contraction)
    assert contraction < Q(3, 10**6)
    assert center+contraction*radius < radius
    assert 9*qbar+center_p+exponential*radius*p_pair < pmax
    assert qbar+length*pmax/4 < qmax
    assert exponential*fixed_point < Q(9, 10**7) < Q(1, 10**6)
    assert Q(2, 10**12)+exponential*radius*alpha < b
    assert v+exponential*radius < m
    assert length/original == 10000
    # The original cutoff is not length/2: the free interval starts at L0/2.
    assert original/2 < length/2


def test_weighted_ball_without_unweighting_is_not_a_physical_self_map():
    sigma, length, radius = sp.symbols("sigma L radius", positive=True)
    envelope = radius*sp.exp(sigma*length)
    assert sp.diff(envelope, length) == sigma*envelope
    # An explicit longer-window control: the same weighted radius can leave
    # the history derivative ball once the exponential loss is retained.
    assert Q(2, 10**7)*Q(8, 3)**20 > Q(1, 10**5)


def test_public_weighted_gate_keeps_original_cutoff_and_independent_margins():
    from p8a_extension import bounds
    from p8a_preparation import bounds as original

    d, old = bounds.calibration(), original.calibration()
    assert d["source_window"] == old["length"] == sp.Rational(1, 10**10)
    assert d["source_flat_end"] == old["length"]/4
    assert d["source_off_from"] == old["length"]/2
    assert d["length"] == sp.Rational(1, 10**6)
    assert d["source_off_from"] != d["length"]/2
    assert d["gate"]["fixed_point_weighted"] < sp.Rational(9, 10**8)
    assert d["gate"]["fixed_point_pointwise"] < original.RADIUS
    assert d["auxiliary_P_derived"] < d["auxiliary_P_cap"]
    assert d["auxiliary_q_derived"] < d["auxiliary_q_cap"]
    assert d["pairs"]["Einstein_full_P"] > 3
    assert d["pairs"]["Einstein_remainder"] > sp.Rational(3, 10**6)
    assert d["pairs"]["actual_mode_response"]["total"] > 0


def test_public_source_resummation_bridge_keeps_the_same_causal_operator():
    from p8a_continuation.resolvent import transfer
    from p8a_extension.map import full_forced_remainder

    s, delta, af = sp.symbols("s delta af", positive=True)
    beta = sp.Symbol("beta", real=True)
    base, rolling, anomaly, auxiliary, nonlinear, source = sp.symbols(
        "base rolling anomaly auxiliary nonlinear source", real=True)
    rhs = base+rolling/(60*delta)-anomaly+2*auxiliary-nonlinear+source
    actual = full_forced_remainder(base, rolling, anomaly, auxiliary, nonlinear, source, delta)
    block = (sp.log(s)+beta-af**2/(30*delta*s**2))/2
    assert sp.simplify(block*transfer(s, beta, af**2/(30*delta))*actual-rhs) == 0
