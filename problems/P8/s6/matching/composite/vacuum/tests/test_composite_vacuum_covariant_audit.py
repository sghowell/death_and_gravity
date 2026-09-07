"""Independent physical-metric, source, mass and elimination-boundary audit."""

import sympy as sp
from p8_composite_vacuum import cayley, matching, vacuum


def test_cayley_dictionary_on_nondiagonal_lorentz_metrics():
    eta = sp.diag(1, -1, -1, -1)
    frame = sp.Matrix([[1, sp.Rational(1, 3), 0, 0], [0, 1, sp.Rational(1, 4), 0],
                       [0, 0, 1, sp.Rational(1, 5)], [0, 0, 0, 1]])
    physical = frame.T*eta*frame
    delta = frame.inv()*sp.diag(sp.Rational(1, 5), sp.Rational(-1, 7),
                               sp.Rational(1, 9), sp.Rational(-1, 11))*frame
    unit = sp.eye(4)
    assert physical*delta != delta*physical  # No Euclidean-commuting shortcut.
    assert physical*delta == (physical*delta).T
    g = physical*(unit+delta)**2/4
    f = physical*(unit-delta)**2/4
    square_root = (unit-delta)*(unit+delta).inv()
    assert g == g.T and f == f.T
    assert g.inv()*f == square_root**2
    assert all(value > 0 for value in square_root.eigenvals())
    assert g+2*g*square_root+f == physical
    assert (unit-square_root)*(unit+square_root).inv() == delta
    assert sp.sqrt(-g.det()) == (unit+delta).det()/16
    assert sp.sqrt(-f.det()) == (unit-delta).det()/16


def test_exact_beta2_potential_includes_the_quartic_term_and_is_even():
    ds = sp.symbols("d0:4", real=True)
    t = sp.Symbol("t")
    numerator = sp.Poly(sp.prod((1+d)+t*(1-d) for d in ds), t).coeff_monomial(t**2)
    trace = sum(ds)
    expected = 6+sum(d*d for d in ds)-trace**2+6*sp.prod(ds)
    assert sp.expand(numerator-expected) == 0
    assert sp.expand(numerator-numerator.subs(dict.fromkeys(ds, 0))) != 0
    assert sp.expand(numerator-numerator.subs({d: -d for d in ds})) == 0
    assert sp.expand(numerator-(6+sum(d*d for d in ds)-trace**2)) == 6*sp.prod(ds)


def test_vacuum_energy_shifts_every_hr_beta_and_changes_the_physical_fp_mass():
    bare = tuple(map(sp.Rational, (0, 0, 1, 0, 0)))
    shifted = tuple(beta-sp.Rational(3, 8) for beta in bare)
    assert shifted == (sp.Rational(-3, 8), sp.Rational(-3, 8), sp.Rational(5, 8),
                       sp.Rational(-3, 8), sp.Rational(-3, 8))
    assert shifted[0]+3*shifted[1]+3*shifted[2]+shifted[3] == 0
    assert shifted[1]+3*shifted[2]+3*shifted[3]+shifted[4] == 0
    # The conventional g-clock mass is divided by four in physical G=4g.
    physical_mass_ratio = 2*(shifted[1]+2*shifted[2]+shifted[3])/4
    wrong_bare_ratio = 2*(bare[1]+2*bare[2]+bare[3])/4
    assert physical_mass_ratio == sp.Rational(1, 4)
    assert wrong_bare_ratio == 1 != physical_mass_ratio


def test_bianchi_i_kinetic_normalization_independently_recovers_mass():
    t = sp.Symbol("T", real=True)
    mass, planck = sp.symbols("m M", positive=True)
    d = sp.Function("d")(t)
    lapse = sp.Rational(1, 2)

    def eh_kinetic(sign):
        scales = [(1+sign*d)/2, (1-sign*d)/2, sp.Rational(1, 2)]
        hubbles = [sp.diff(a, t)/(lapse*a) for a in scales]
        cross = sum(hubbles[i]*hubbles[j] for i in range(3) for j in range(i+1, 3))
        return sp.factor(-planck**2*lapse*sp.prod(scales)*cross)

    kinetic = sp.factor(eh_kinetic(1)+eh_kinetic(-1))
    assert kinetic == planck**2*sp.diff(d, t)**2/2
    potential = -planck**2*mass**2*(2*d*d)/16
    lagrangian = kinetic+potential
    euler = sp.diff(sp.diff(lagrangian, sp.diff(d, t)), t)-sp.diff(lagrangian, d)
    assert sp.simplify(euler/planck**2-sp.diff(d, t, 2)-mass**2*d/4) == 0


def test_asymmetric_metric_matter_would_source_the_relative_field():
    ds = sp.symbols("d0:4", real=True)
    velocity = sp.Symbol("velocity", real=True)
    # A scalar coupled to g alone, with physical G=eta fixed.
    density_on_g = sp.prod(1+d for d in ds)*velocity**2/(8*(1+ds[0])**2)
    zero = dict.fromkeys(ds, 0)
    assert sp.diff(density_on_g, ds[0]).subs(zero) == -velocity**2/8
    # The stipulated physical-G scalar has no Delta dependence at all.
    density_on_physical = velocity**2/2
    assert all(sp.diff(density_on_physical, d) == 0 for d in ds)


def test_positive_algebraic_mass_is_not_invertibility_of_a_boundary_operator():
    t, mass = sp.symbols("T m", positive=True)
    mode = sp.sin(mass*t)
    assert sp.diff(mode, t, 2)+mass**2*mode == 0
    assert mode.subs(t, 0) == mode.subs(t, sp.pi/mass) == 0
    assert sp.diff(mode, t).subs(t, 0) == mass != 0
    # Dirichlet data on the resonant interval have a nonzero kernel;
    # this is not a zero-data retarded solution.


def test_even_stationarity_does_not_exclude_other_branches_or_loop_corrections():
    delta, parameter, x = sp.symbols("Delta parameter X", real=True)
    mass = sp.Symbol("m", positive=True)
    potential = delta**4/4-parameter*delta**2/2
    equation = sp.diff(potential, delta)
    assert equation.subs(delta, 0) == 0
    assert sp.factor(equation) == delta*(delta**2-parameter)
    assert sp.diff(equation, delta).subs({delta: 0, parameter: 0}) == 0
    determinant_term = sp.log(mass**2+x)/2
    assert sp.diff(determinant_term, x, 2).subs(x, 0) == -1/(2*mass**4)


def test_cd_physical_time_and_operator_gaps_cannot_be_removed_by_small_errors():
    t, tau = sp.symbols("T tau", real=True, positive=True)
    x = sp.Symbol("X", positive=True)
    scale = (1+(t/tau)**2)**2
    hubble = sp.factor(sp.diff(scale, t)/scale)
    assert hubble == 4*t/(t**2+tau**2)
    assert sp.diff(hubble, t).subs(t, 0) == 4/tau**2
    assert sp.diff(hubble, t).subs(t, tau/2) == sp.Rational(48, 25)/tau**2
    assert (hubble.subs(t, tau/2)-hubble.subs(t, -tau/2))*tau/2 == sp.Rational(8, 5)
    # Original dimensionless CD clock and physical-metric coefficient chart.
    denominator = (1+t*t)**3
    f2 = -sp.Rational(1, 2)+(1-x)/(2*denominator)
    a3 = 1/(denominator*x)
    assert sp.diff(f2, x).subs(t, 0) == -sp.Rational(1, 2)
    assert a3.subs({t: 0, x: 1}) == 1


def test_full_four_dimensional_fp_divergence_trace_and_five_polarizations():
    eta = sp.diag(1, -1, -1, -1)
    wave = sp.Matrix(sp.symbols("k0:4", real=True))
    raised = eta*wave
    values = sp.symbols("h00 h01 h02 h03 h11 h12 h13 h22 h23 h33", real=True)
    h = sp.Matrix([[values[0], values[1], values[2], values[3]],
                   [values[1], values[4], values[5], values[6]],
                   [values[2], values[5], values[7], values[8]],
                   [values[3], values[6], values[8], values[9]]])
    mass2 = sp.Symbol("m_FP_squared", positive=True)
    k2, trace = (wave.T*raised)[0], sp.trace(eta*h)
    divergence = h*raised
    double = (raised.T*h*raised)[0]
    kinetic = -k2*h-trace*(wave*wave.T)+wave*divergence.T+divergence*wave.T
    kinetic += eta*(-double+k2*trace)
    assert all(sp.expand(entry) == 0 for entry in kinetic*raised)
    equation = kinetic+mass2*(h-eta*trace)
    assert all(sp.expand(entry) == 0 for entry in equation*raised-mass2*(divergence-wave*trace))
    assert sp.expand(sp.trace(eta*equation)+2*(double-k2*trace)+3*mass2*trace) == 0
    # The four divergence constraints and the trace leave five modes in
    # a timelike rest frame, not six unconstrained spatial components.
    rest = dict(zip(wave, (1, 0, 0, 0), strict=True))
    constraints = [entry.subs(rest) for entry in divergence-wave*trace]+[trace]
    constraint_matrix, _ = sp.linear_eq_to_matrix(constraints, values)
    assert constraint_matrix.rank() == 5
    assert len(values)-constraint_matrix.rank() == 5


def test_physical_einstein_normalization_and_literal_tt_stress_source():
    t = sp.Symbol("T", real=True)
    scale, lapse = sp.Function("a")(t), sp.Symbol("N", positive=True)
    planck = sp.Symbol("M", positive=True)
    each_einstein = -3*planck**2*(scale/2)*(sp.diff(scale, t)/2)**2/(lapse/2)
    physical_einstein = -3*(planck**2/2)*scale*sp.diff(scale, t)**2/lapse
    assert sp.simplify(2*each_einstein-physical_einstein) == 0
    stress, gamma, source = sp.symbols("Pi gamma j", real=True)
    literal_probe = -stress*gamma/2
    assert (source*gamma/4).subs(source, -2*stress) == literal_probe


def test_primary_dictionary_matches_independent_mass_and_cd_gap_audits():
    c, v, m = cayley.derive(), vacuum.derive(), matching.derive()
    assert c["zero_relative_gravity_coefficient"] == -c["M2"]/4
    assert v["physical_Planck_squared"] == v["M2"]/2
    assert v["relative_Planck_squared"] == 2*v["M2"]
    assert v["physical_FP_mass_squared"] == v["m2"]/4
    assert v["physical_massive_source_response"] == 0
    assert m["target_F2X_at_centre"] == -sp.Rational(1, 2)
    assert m["target_A3_at_centre"] == 1
    assert m["zero_relative_F2X"] == m["zero_relative_A3"] == 0
    # Keeping the separate rolling free-M1 field costs another1/100 in
    # the normalized null equation, over the clock's nonnegative term.
    assert m["exact_null_residual_lower_magnitude"] == 8+sp.Rational(1, 100)
    assert m["endpoint_H_error_threshold"] == sp.Rational(8, 5)
