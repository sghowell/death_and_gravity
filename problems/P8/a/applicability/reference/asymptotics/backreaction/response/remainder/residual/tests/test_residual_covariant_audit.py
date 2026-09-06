"""Independent proper-time curvature, anomaly and finite-residual audits."""

import sympy as sp


def _geometry():
    eta = sp.Symbol("eta", positive=True)
    a = sp.Function("a", positive=True)(eta)
    h = sp.diff(a, eta)/a
    physical_h = sp.diff(a, eta)/a**2
    physical_h_dot = sp.diff(physical_h, eta)/a
    curvature = 6*(physical_h_dot+2*physical_h**2)
    kretschmann = 12*((physical_h_dot+physical_h**2)**2+physical_h**4)
    ricci_squared = 12*(physical_h_dot**2+3*physical_h_dot*physical_h**2+3*physical_h**4)
    box = lambda value: (sp.diff(value, eta, 2)+2*h*sp.diff(value, eta))/a**2
    anomaly = curvature**2/288+(kretschmann-ricci_squared)/720-box(curvature)/120
    return eta, a, h, -sp.diff(a, eta, 2)/a, box, anomaly


def _reconstructed(a, h, u, s, j, history, eta):
    # Strip only hbar/pi^2, not a^-4. No public reconstruction is used here.
    density = (h**4/960+(-h*sp.diff(s, eta)+(h**2-u)*s+j)/8-history/32)/a**4
    pressure = ((5*h**4+4*h**2*u)/2880
                +(sp.diff(s, eta, 2)-3*h*sp.diff(s, eta)+(3*h**2+u)*s+j)/24
                +(u**2-history)/96)/a**4
    effective = ((3*h**4+2*h**2*u)/960
                 +(sp.diff(s, eta, 2)-4*h*sp.diff(s, eta)+4*h**2*s+2*j)/16
                 +(u**2-2*history)/64)/a**4
    return density, pressure, effective


def test_actual_trace_from_independent_proper_time_curvature_invariants():
    eta, a, h, u, box, anomaly = _geometry()
    mode, kernel, j, history = (sp.Function(name)(eta) for name in ("Rmode", "K", "J", "L"))
    wick = (mode-kernel/2)/(4*a**2)
    s = mode-kernel/2+u/10
    density, pressure, effective = _reconstructed(a, h, u, s, j, history, eta)
    trace = -box(wick)/2-anomaly/4
    assert sp.simplify(density-3*pressure-trace) == 0
    assert sp.simplify(effective-(density+3*pressure)/2) == 0
    # Dropping the Box-R part of the anomaly changes a generic metric.
    wrong_s = mode-kernel/2
    wrong_density, wrong_pressure, _ = _reconstructed(a, h, u, wrong_s, j, history, eta)
    assert sp.simplify(wrong_density-3*wrong_pressure-trace) != 0


def test_two_distinct_retarded_histories_are_required_by_conservation():
    eta, a, h, u, _, _ = _geometry()
    s, j, history = (sp.Function(name)(eta) for name in ("S", "J", "L"))
    density, pressure, _ = _reconstructed(a, h, u, s, j, history, eta)
    conservation = sp.diff(density, eta)+3*h*(density+pressure)
    rules = {sp.diff(j, eta): sp.diff(u, eta)*s, sp.diff(history, eta): h*u**2}
    assert sp.simplify(conservation.subs(rules)) == 0
    assert sp.simplify(conservation.subs({**rules, sp.diff(j, eta): 0})) != 0
    assert sp.simplify(conservation.subs({**rules, sp.diff(history, eta): 0})) != 0


def test_common_past_fixes_radiation_constant_not_trace_alone():
    eta, a, h, u, _, _ = _geometry()
    s, j, history = (sp.Function(name)(eta) for name in ("S", "J", "L"))
    density, pressure, effective = _reconstructed(a, h, u, s, j, history, eta)
    scale = sp.Symbol("A", positive=True)
    past = {a: scale*eta, s: 0, j: 0, history: 0}
    assert sp.simplify(density.subs(past).doit()-1/(960*scale**4*eta**8)) == 0
    assert sp.simplify(pressure.subs(past).doit()-1/(576*scale**4*eta**8)) == 0
    assert sp.simplify(effective.subs(past).doit()-1/(320*scale**4*eta**8)) == 0
    constant = sp.Symbol("C", real=True)
    # Trace and conservation cannot reject this extra radiation by themselves.
    extra_density, extra_pressure = constant/a**4, constant/(3*a**4)
    assert extra_density-3*extra_pressure == 0
    assert sp.simplify(sp.diff(extra_density, eta)+3*h*(extra_density+extra_pressure)) == 0
    assert extra_density.subs(a, scale*eta) != 0


def test_exact_frozen_defect_and_literal_1922_plateau_enclosure():
    time, scale, d, epsilon = sp.symbols("t A d epsilon", positive=True)
    a = (4*scale**2*(time**2-4*epsilon*d))**sp.Rational(1, 4)
    hubble = sp.diff(a, time)/a
    kappa_radiation = 3*scale**2/a**4
    frozen_rho = sp.factor(3*hubble**2-kappa_radiation-3*epsilon*d/time**4)
    frozen_p = sp.factor(-2*sp.diff(hubble, time)-3*hubble**2-kappa_radiation/3
                         -5*epsilon*d/time**4)
    difference = (time**2-4*epsilon*d)**(-2)-time**(-4)
    assert sp.factor(frozen_rho-3*epsilon*d*difference) == 0
    assert sp.factor(frozen_p-5*epsilon*d*difference) == 0
    assert sp.factor((frozen_rho+3*frozen_p)/2-9*epsilon*d*difference) == 0
    z = sp.Symbol("z", positive=True)
    quotient = sp.factor(((1-z)**(-2)-1)/z)
    assert sp.factor(quotient-(2-z)/(1-z)**2) == 0
    assert sp.limit(quotient, z, 0) == 2
    assert sp.Rational(2048, 961)/2**12 == sp.Rational(1, 1922)
    # Both separately enlarged positive factors, not an endpoint-only sample.
    assert 2*(1-sp.Rational(1, 32))**(-2) == sp.Rational(2048, 961)


def test_coupling_and_reference_scales_use_the_same_proper_time():
    a_scale, eta_star, y, epsilon = sp.symbols("A eta_star y epsilon", positive=True)
    d, kappa, hbar = sp.symbols("d kappa hbar", positive=True)
    substitutions = {d: kappa*hbar/(46080*sp.pi**2)}
    amplitude = 16*epsilon*d/(a_scale**2*eta_star**4)
    assert sp.simplify(2880*amplitude*a_scale**2*eta_star**4
                       -epsilon*kappa*hbar/sp.pi**2).subs(substitutions).simplify() == 0
    time = a_scale*eta_star**2*y**2/2
    reference = hbar/(15360*sp.pi**2*time**4)
    assert sp.simplify(reference-hbar/(960*sp.pi**2*a_scale**4*eta_star**8*y**8)) == 0
    preparation_time = a_scale*eta_star**2/2
    assert sp.simplify(amplitude-4*epsilon*d/preparation_time**2) == 0
    t_star = sp.Symbol("t_star", positive=True)
    assert (4*d/t_star**2).subs(t_star, 2*10**7*sp.sqrt(d)) == sp.Rational(1, 10**14)


def test_public_reconstruction_matches_independent_invariants_and_past_freedom():
    from p8a_residual import reconstruction

    eta, a, h, u, box, anomaly = _geometry()
    mode, kernel, j, history = (sp.Function(name)(eta) for name in ("Rmode", "K", "J", "L"))
    hbar, gamma, constant = sp.symbols("hbar gamma C", real=True)
    s = mode-kernel/2+u/10
    actual = reconstruction.components(a, s, j, history, eta, hbar=hbar,
                                       gamma=gamma, radiation_constant=constant)
    wick = hbar*(mode-kernel/2)/(4*sp.pi**2*a**2)
    curvature = -6*u/a**2
    expected_trace = -box(wick)/2-hbar*anomaly/(4*sp.pi**2)-6*gamma*box(curvature)
    assert sp.simplify(actual["density"]-3*actual["pressure"]-expected_trace) == 0
    assert sp.simplify(reconstruction.trace(a, wick, eta, hbar=hbar, gamma=gamma)-expected_trace) == 0
    raw = _reconstructed(a, h, u, s, j, history, eta)
    unshifted = reconstruction.components(a, s, j, history, eta, hbar=hbar)
    for name, value in zip(("density", "pressure", "EED"), raw, strict=True):
        assert sp.simplify(unshifted[name]-hbar*value/sp.pi**2) == 0
    history_rules = {sp.diff(j, eta): sp.diff(u, eta)*s,
                     sp.diff(history, eta): h*u**2}
    conservation = sp.diff(actual["density"], eta)+3*h*(actual["density"]+actual["pressure"])
    assert sp.simplify(conservation.subs(history_rules)) == 0
    assert reconstruction.effective_s(mode, kernel, u) == s
    assert sp.simplify(reconstruction.fix_radiation_constant(a, constant/a**4, 0)-constant) == 0
