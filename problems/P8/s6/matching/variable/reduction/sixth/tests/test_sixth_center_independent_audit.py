"""Independent raw-action, canonical-clock and center cancellation audit.

No claimed tensor response on the full CD FLRW background is inferred from
the off-shell unit-volume restriction.  Time differentiation precedes all
center substitutions, including the canonical-clock pullback.
"""

from functools import cache

import sympy as sp
from p8_sixth_reduction import literal, tensor


@cache
def _implicit_root_jets():
    v, c = sp.symbols("v c", real=True)
    r1, r2 = sp.symbols("root_linear root_quadratic", real=True)
    d = 1+v
    w = c*(1-7*v)*d**4+12*v
    # Solve the original cleared root equation coefficient by coefficient,
    # not a fractional power or the primary quotient-jet recurrence.
    root_trial = 2+r1*v+r2*v*v
    equation = sp.Poly(sp.expand(root_trial**3*d**8*w-8*c*(1-v)), v)
    r1_value = sp.factor(sp.solve(equation.nth(1), r1)[0])
    r2_value = sp.factor(sp.solve(equation.nth(2).subs(r1, r1_value), r2)[0])
    root = 2+r1_value*v+r2_value*v*v
    abar = sp.series(root*c*d**10*(c*d**4-2)/(32*(1-v)), v, 0, 3).removeO().expand()
    y = 2/d**4
    hp = 4*(1-v)/d**2
    b1 = y**3*hp/(c*(c-y))
    b4 = 24*v/(c*c*d*d)-b1/y**3
    kbar = 2*(y**3/c-1)*hp-1/(100*d**12)
    return {"v": v, "c": c, "r1": r1_value, "r2": r2_value,
            "r_cubed": 8*c*(1-v)/(d**8*w), "b1": b1, "b4": b4,
            "A0": sp.factor(abar.coeff(v, 0)),
            "A2": sp.factor(2*abar.coeff(v, 1)),
            "A4": sp.factor(24*abar.coeff(v, 2)),
            "ell1": r1_value, "ell3": sp.factor(12*r2_value-3*r1_value**2),
            "k0": sp.factor(kbar.subs(v, 0)),
            "k2": sp.factor(2*sp.diff(kbar, v).subs(v, 0))}


@cache
def _raw_polynomial_action():
    """Direct even-background polynomial path through precisely needed jets.

    This uses ordinary polynomial derivatives and the exact diagonal ADM
    volume/rate expansion.  It does not call either primary jet derivative,
    primary B formula, or primary Euler implementation.
    """
    u = sp.Symbol("audit_u", real=True)
    a0, a2, a4, ell1, ell3, r0 = sp.symbols("a0 a2 a4 ell1 ell3 r0_squared", real=True)
    q1, q2, q3 = sp.symbols("audit_q1 audit_q2 audit_q3", real=True)

    def truncate(expression, degree=4):
        polynomial = sp.Poly(sp.expand(expression), u, q1, q2, q3)
        return sp.Add(*(coefficient*u**power[0]*q1**power[1]*q2**power[2]*q3**power[3]
                        for power, coefficient in polynomial.terms()
                        if power[0] <= degree and sum(power[1:]) <= 2))

    def product(*expressions, degree=4):
        result = sp.S.One
        for expression in expressions:
            result = truncate(result*expression, degree)
        return result

    abar = a0+a2*u*u/2+a4*u**4/24
    ell = ell1*u+ell3*u**3/6
    r_squared = r0*(1+ell1*u*u+(ell1**2/2+ell3/12)*u**4)
    schouten = (-5*q1*q1/sp.Integer(12)-2*sp.diff(ell, u)+ell*ell,
                -q2/2+q1*q1/12-ell*q1-ell*ell,
                q2/2+q1*q1/12+ell*q1-ell*ell,
                q1*q1/12-ell*ell)
    relative = tuple(product(abar, value) for value in schouten)
    derivative = tuple(truncate(sp.diff(value, u)+q2*sp.diff(value, q1)
                                +q3*sp.diff(value, q2), 3) for value in relative)
    volume1 = truncate((sum(relative[1:])-relative[0])/2, 3)
    volume2 = truncate(product(sum(relative[1:])-relative[0], sum(relative[1:])-relative[0], degree=3)/8
                       -(sum(product(value, value, degree=3) for value in relative[1:])
                         -product(relative[0], relative[0], degree=3))/4, 3)
    base_rates = (ell+q1/2, ell-q1/2, ell)
    pair0 = 3*ell*ell-q1*q1/4
    pair1 = sum(product(derivative[i+1], 3*ell-base_rates[i], degree=3) for i in range(3))/2
    pair2 = (sum(product(derivative[i+1], derivative[j+1], degree=3)
                 for i in range(3) for j in range(i+1, 3))/4
             -sum(product(relative[i+1], derivative[i+1], 3*ell-base_rates[i], degree=3)
                  for i in range(3))/2)
    einstein = -product(r_squared, pair2+product(volume1, pair1, degree=3)
                        +product(volume2, pair0, degree=3), degree=3)
    newton = (product(sum(schouten), sum(schouten), sum(schouten), degree=3)
              -6*product(sum(schouten), sum(product(value, value, degree=3) for value in schouten), degree=3)
              +5*sum(product(value, value, value, degree=3) for value in schouten))
    potential = product(r_squared, abar, abar, newton, degree=3)/24
    fourth = product(r_squared, abar,
                     sum(product(value, value, degree=3) for value in schouten)
                     -product(sum(schouten), sum(schouten), degree=3), degree=3)/4

    def center_euler(raw):
        poly = sp.Poly(sp.expand(raw), q1, q2, q3)
        a11, a12, a13, a22, a23, a33 = (poly.coeff_monomial(term) for term in
                                       (q1*q1, q1*q2, q1*q3, q2*q2, q2*q3, q3*q3))
        # Product-rule Euler coefficients, not frozen center Lagrangians.
        e2 = (-2*a11+sp.diff(a12, u)-3*sp.diff(a13, u, 2)
              +2*sp.diff(a22, u, 2)-sp.diff(a23, u, 3)).subs(u, 0)
        e4 = (2*a22-2*a13-sp.diff(a23, u)-6*sp.diff(a33, u, 2)).subs(u, 0)
        e6 = (-2*a33).subs(u, 0)
        return {2: sp.factor(e2), 4: sp.factor(e4), 6: sp.factor(e6)}

    raw = truncate(einstein+potential, 3)
    return {"u": u, "a0": a0, "a2": a2, "a4": a4, "ell1": ell1,
            "ell3": ell3, "r0_squared": r0, "q": (q1, q2, q3),
            "sixth_raw": raw, "E6": center_euler(raw),
            "EH_E": center_euler(einstein), "V3_E": center_euler(potential),
            "E4": center_euler(fourth)}


def _actual_point():
    data, profile = _raw_polynomial_action(), _implicit_root_jets()
    return {data["a0"]: profile["A0"], data["a2"]: profile["A2"], data["a4"]: profile["A4"],
            data["ell1"]: profile["ell1"], data["ell3"]: profile["ell3"], data["r0_squared"]: 4}


def test_cleared_root_and_all_required_actual_time_jets():
    data, primary = _implicit_root_jets(), tensor.center_jets()
    assert sp.factor(-data["b1"]/data["b4"]-data["r_cubed"]) == 0
    point = {data["c"]: tensor.C}
    for name, entry in zip(("A0", "A2", "A4"), primary["A_bar_u0_u2_u4"], strict=True):
        assert sp.factor(data[name].subs(point)-entry) == 0
    for name, entry in zip(("ell1", "ell3"), primary["log_r_u1_u3"], strict=True):
        assert sp.factor(data[name].subs(point)-entry) == 0
    assert sp.factor(data["A4"]-6*(22*data["c"]**3-39*data["c"]**2+32*data["c"]-16)/data["c"]) == 0


def test_center_profiles_match_both_primary_and_literal_routes():
    own, other = _implicit_root_jets(), literal.independent_center()
    for order, name in ((0, "A0"), (2, "A2"), (4, "A4")):
        coefficient = sp.diff(other["A_bar_series"], other["v"], order//2).subs(other["v"], 0)
        coefficient *= {0: 1, 2: 2, 4: 12}[order]
        assert sp.factor(coefficient-own[name].subs(own["c"], other["c"])) == 0


def test_ordinary_raw_action_center_Euler_formula():
    d = _raw_polynomial_action()
    a, a2, a4, ell1, ell3, r0 = (d[key] for key in ("a0", "a2", "a4", "ell1", "ell3", "r0_squared"))
    expected = {2: -r0*(a2*a2+22*a*a2*ell1+26*a*a*ell1**2+2*a*a*ell3+a*a4)/8,
                4: -r0*a*(8*a*ell1+7*a2)/8, 6: -r0*a*a/8}
    for n in (2, 4, 6):
        assert sp.factor(d["E6"][n]-expected[n]) == 0
        actual = sp.factor(d["E6"][n].subs(_actual_point()))
        assert sp.factor(actual.subs(_implicit_root_jets()["c"], tensor.C)
                         -tensor.center()["sixth_E_dimensionless"][n]) == 0


def test_Einstein_Hessian_cannot_be_replaced_by_the_potential_cubic():
    d = _raw_polynomial_action()
    assert sp.factor(d["V3_E"][2].subs(d["a0"], 0)) == 0
    assert sp.factor(d["EH_E"][2].subs(d["a0"], 0)
                     +d["r0_squared"]*d["a2"]**2/8) == 0
    assert d["E6"][2] != d["V3_E"][2]


def test_raw_action_time_parity_does_not_require_freezing_background_jets():
    d = _raw_polynomial_action()
    q1, _, q3 = d["q"]
    reflection = {d["u"]: -d["u"], q1: -q1, q3: -q3}
    assert sp.expand(d["sixth_raw"].subs(reflection, simultaneous=True)-d["sixth_raw"]) == 0
    assert all(tensor.center()["sixth_E_dimensionless"][n] == 0 for n in (0, 1, 3, 5))


def test_actual_canonical_clock_is_not_the_unit_time_label():
    d = _implicit_root_jets()
    c = d["c"]
    assert sp.factor(d["k0"]-(64/c-sp.Rational(801, 100))) == 0
    assert sp.factor(d["k2"]-(sp.Rational(1206, 25)-1920/c)) == 0
    assert d["k0"].subs(c, 3) != 1
    assert d["k2"].subs(c, 3) != 0


def test_nonlinear_canonical_clock_fourth_jet_round_trip():
    d = _implicit_root_jets()
    k0, k2 = d["k0"], d["k2"]
    # F is the same physical-component function, not a reindexed tensor.
    f2, f4 = d["A2"], d["A4"]
    f_varphi2 = f2/k0
    f_varphi4 = f4/k0**2-2*f2*k2/k0**3
    assert sp.factor(k0**2*f_varphi4+2*k2*f_varphi2-f4) == 0
    incorrectly_linear_clock = k0**2*f_varphi4
    assert sp.factor((incorrectly_linear_clock-f4).subs(d["c"], 3)) != 0
    assert sp.factor(k0*f_varphi2-f2) == 0


def test_dimensionful_clock_pullback_cancels_M_and_tau_correctly():
    d = _implicit_root_jets()
    m, tau = sp.symbols("M tau", positive=True)
    k0, k2 = d["k0"], d["k2"]
    phi_T = m*sp.sqrt(k0)/tau
    phi_TTT = m*k2/(2*tau**3*sp.sqrt(k0))
    A_phi2 = tau*tau*d["A2"]/(m*m*k0)
    A_phi4 = tau*tau/m**4*(d["A4"]/k0**2-2*d["A2"]*k2/k0**3)
    assert sp.factor(A_phi2*phi_T**2-d["A2"]) == 0
    assert sp.factor(A_phi4*phi_T**4+4*A_phi2*phi_T*phi_TTT-d["A4"]/tau**2) == 0


def test_all_density_and_Euler_terms_have_the_physical_dimensions():
    tau = sp.Symbol("audit_tau", positive=True)
    scale = {tensor.INV[n]: tau**(2-n)*tensor.INV[n] for n in range(10)}
    scale.update({tensor.LOG[n]: tau**(-1-n)*tensor.LOG[n] for n in range(10)})
    d = tensor.derive()
    for n in (1, 2, 3):
        actual = d[f"B{n}"].subs(scale, simultaneous=True)
        assert sp.factor(actual-tau**(2*n-2)*d[f"B{n}"]) == 0
    for order, coefficient in d["euler_coefficients"].items():
        assert sp.factor(coefficient.subs(scale, simultaneous=True)-tau**(order-2)*coefficient) == 0
    center = tensor.center()
    for order in (2, 4, 6):
        assert sp.factor(center["sixth_E_physical"][order]
                         -tensor.M**2*tensor.TAU**(order-2)*center["sixth_E_dimensionless"][order]) == 0
    assert center["sixth_E_physical"][2].has(tensor.M)
    assert not center["sixth_E_physical"][2].has(tensor.TAU)


def test_continuous_isolated_sixth_floor_is_not_a_sampled_bound():
    d, delta = _raw_polynomial_action(), sp.Symbol("delta", nonnegative=True)
    actual = sp.factor(d["E6"][2].subs(_actual_point()))
    expression = sp.expand(actual.subs(_implicit_root_jets()["c"], 2+delta))
    positive = sp.factor(-(expression+2-sp.Rational(7, 4)*delta)/delta**2)
    assert sp.factor(positive-(349*delta**2+572*delta+836)/128) == 0
    assert all(coefficient > 0 for coefficient in sp.Poly(positive, delta).all_coeffs())
    cap = sp.Rational(1, 100)
    floor = 2-sp.Rational(7, 4)*cap
    assert floor == sp.Rational(793, 400)
    assert sp.expand(-expression-floor-sp.Rational(7, 4)*(cap-delta)-delta**2*positive) == 0
    assert floor/5 == sp.Rational(793, 2000)


def test_fourth_order_is_derived_from_its_own_stationary_square():
    d, point = _raw_polynomial_action(), _actual_point()
    expected = d["r0_squared"]*(d["a2"]+2*d["a0"]*d["ell1"])/4
    assert sp.factor(d["E4"][2]-expected) == 0
    own = sp.factor(d["E4"][2].subs(point))
    assert sp.factor(own.subs(_implicit_root_jets()["c"], tensor.C)
                     -tensor.center()["fourth_E_dimensionless"][2]) == 0


def test_mandatory_S4_S6_cancellation_and_small_positive_combined_symbol():
    d, profile = _raw_polynomial_action(), _implicit_root_jets()
    e4, e6 = (sp.factor(d[key][2].subs(_actual_point())) for key in ("E4", "E6"))
    assert e4.subs(profile["c"], 2) == 2
    assert e6.subs(profile["c"], 2) == -2
    assert sp.factor((e4+e6).subs(profile["c"], 2)) == 0
    delta = sp.Symbol("delta", nonnegative=True)
    combined = sp.expand((e4+e6).subs(profile["c"], 2+delta))
    quotient = sp.cancel(combined/delta)
    assert quotient == (-349*delta**3-572*delta**2-692*delta+448)/128
    # On 0<delta<=.01 this combined contribution is O(delta), not O(1).
    assert quotient.subs(delta, sp.Rational(1, 100)) > sp.Rational(17, 5)
    assert sp.Poly(sp.Rational(7, 2)-quotient, delta).all_coeffs() == [
        sp.Rational(349, 128), sp.Rational(143, 32), sp.Rational(173, 32), 0]
    assert sp.factor(combined-tensor.center()["combined_E_dimensionless"][2].subs(tensor.C, 2+delta)) == 0


def test_freezing_center_coefficients_misses_the_low_symbol():
    data = tensor.center()
    frozen = -2*data["B_dimensionless"][1]
    assert frozen.subs(tensor.C, 2) == 0
    assert data["sixth_E_dimensionless"][2].subs(tensor.C, 2) == -2
    own = _raw_polynomial_action()
    assert sp.factor(own["E6"][2].subs(own["a0"], 0)) == -own["r0_squared"]*own["a2"]**2/8


def test_unit_volume_probe_is_not_full_FLRW_even_at_zero_Hubble():
    time, tau = sp.symbols("T tau", real=True, positive=True)
    a, r = sp.Function("a")(time), sp.Function("r")(time)
    scale_f, lapse_f = a*r, r
    ricci00 = -3*(sp.diff(scale_f, time, 2)/scale_f
                  -sp.diff(scale_f, time)*sp.diff(lapse_f, time)/(scale_f*lapse_f))
    scalar_times_r2 = -6*(sp.diff(scale_f, time, 2)/scale_f
                          +sp.diff(scale_f, time)**2/scale_f**2
                          -sp.diff(scale_f, time)*sp.diff(lapse_f, time)/(scale_f*lapse_f))
    full00 = ricci00-scalar_times_r2/6
    ell = sp.diff(r, time)/r
    probe00 = -2*sp.diff(ell, time)+ell**2
    hubble = sp.diff(a, time)/a
    assert sp.factor(full00-probe00+2*sp.diff(hubble, time)+hubble**2) == 0
    cd_scale = (1+(time/tau)**2)**2
    difference = (-2*sp.diff(hubble, time)-hubble**2).subs(a, cd_scale).doit()
    assert sp.factor(difference.subs(time, 0)+8/tau**2) == 0
    assert sp.factor(hubble.subs(a, cd_scale).doit().subs(time, 0)) == 0


def test_scope_keeps_probe_symbol_distinct_from_full_Xi_or_remainder():
    scope = tensor.calibration()["scope"]
    assert "unit-volume" in scope
    assert "not Xi" in scope
    assert "full CD propagation" in scope
    assert "all-order estimate" in scope
    assert "UV exclusion" in scope
    assert "2<c<=4" in tensor.center()["domain"]
