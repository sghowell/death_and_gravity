"""Independent coordinate and clock-jet audit of the stationary-f correction.

The audit does not assume that the formal stationary-f series approximates
the actual rolling solution.  Physical components and clock labels remain
distinct throughout.
"""

from functools import cache

import sympy as sp
from p8_variable_reduction import stationary


@cache
def _literal_conformal_schouten():
    u = sp.Symbol("u", real=True)
    a, r = sp.Function("a")(u), sp.Function("r")(u)
    metric = sp.diag(r*r, -r*r*a*a, -r*r*a*a, -r*r*a*a)
    inverse = metric.inv()

    def derivative(expr, index):
        return sp.diff(expr, u) if index == 0 else sp.S.Zero

    gamma = [[[sp.factor(sum(inverse[rho, lam]*(
        derivative(metric[lam, nu], mu)+derivative(metric[lam, mu], nu)
        -derivative(metric[mu, nu], lam)) for lam in range(4))/2)
        for nu in range(4)] for mu in range(4)] for rho in range(4)]
    ricci = sp.zeros(4)
    for mu in range(4):
        for nu in range(4):
            ricci[mu, nu] = sp.factor(sum(
                derivative(gamma[rho][mu][nu], rho)
                -derivative(gamma[rho][mu][rho], nu)
                +sum(gamma[rho][rho][lam]*gamma[lam][mu][nu]
                     -gamma[rho][nu][lam]*gamma[lam][mu][rho]
                     for lam in range(4)) for rho in range(4)))
    scalar = sp.factor(sp.trace(inverse*ricci))
    schouten = (ricci-metric*scalar/6).applyfunc(sp.factor)
    return u, a, r, metric, ricci, scalar, schouten


@cache
def _independent_profile():
    v, c = sp.symbols("v c", real=True)
    d = 1+v
    w = c*(1-7*v)*d**4+12*v
    r3 = 8*c*(1-v)/(d**8*w)
    # Differentiate the explicit logarithm, rather than the parent's r3.
    ell = -sp.Rational(2, 3)*(1/(1-v)+8/d+sp.diff(w, v)/w)
    kap = c*c*d*d*(c*d**4-2)/(16*w)
    hp, h2 = 4*(1-v)/d**2, 16*v/d**2
    schouten00 = -2*hp-h2-2*ell-4*v*sp.diff(ell, v)+v*ell**2
    correction00 = 4*kap*schouten00
    spatial = 4*kap*v*(4/d+ell)**2

    def even_jets(expr):
        return tuple(sp.factor(factor*sp.diff(expr, v, order).subs(v, 0))
                     for order, factor in ((0, 1), (1, 2), (2, 12)))

    # r²=4 exp[integral_0^v ell(z) dz], avoiding any symbolic root branch.
    r2_jets = (sp.Integer(4), sp.factor(8*ell.subs(v, 0)),
               sp.factor(48*(sp.diff(ell, v)+ell**2).subs(v, 0)))
    correction_jets = even_jets(correction00)
    remainder = (sp.factor(c*c-r2_jets[0]-correction_jets[0]),
                 sp.factor(-r2_jets[1]-correction_jets[1]),
                 sp.factor(-r2_jets[2]-correction_jets[2]))
    k0 = 64/c-sp.Rational(801, 100)
    k2 = sp.Rational(1206, 25)-1920/c
    return {"v": v, "c": c, "W": w, "r3": r3, "ell": ell,
            "kap": kap, "P00": schouten00, "f2": correction00,
            "spatial": spatial, "f2_jets": correction_jets,
            "r2_jets": r2_jets, "remainder": remainder,
            "spatial_jets": even_jets(spatial), "k0": k0, "k2": k2}


def test_literal_four_dimensional_curvature_sign_and_conformal_schouten():
    u, a, r, _, _, scalar, schouten = _literal_conformal_schouten()
    hubble, log_derivative = sp.diff(a, u)/a, sp.diff(r, u)/r
    expected00 = -2*(sp.diff(hubble, u)+sp.diff(log_derivative, u))-hubble**2+log_derivative**2
    assert sp.factor(schouten[0, 0]-expected00) == 0
    for index in (1, 2, 3):
        assert sp.factor(schouten[index, index]/a**2-(hubble+log_derivative)**2) == 0
    expected_scalar = -6*(sp.diff(hubble, u)+2*hubble**2)
    assert sp.factor(scalar.subs(r, 1).doit()-expected_scalar) == 0


def test_trace_subtraction_uses_full_not_half_schouten():
    _, _, _, metric, ricci, scalar, schouten = _literal_conformal_schouten()
    assert schouten == (ricci-metric*scalar/6).applyfunc(sp.factor)
    assert sp.factor(schouten[0, 0]/2-schouten[0, 0]) != 0


def test_profile_root_from_original_f_equation_not_rolling_ratio():
    data = _independent_profile()
    v, c = data["v"], data["c"]
    d = 1+v
    y = 2/d**4
    b1 = y**3*4*(1-v)/(d*d*c*(c-y))
    b4 = 24*v/(c*c*d*d)-b1/y**3
    assert sp.factor(-b1/b4-data["r3"]) == 0
    assert sp.factor(-1/(4*b4)-data["kap"]) == 0
    assert sp.factor(data["r3"]-y**3) != 0
    assert sp.factor(data["r3"].subs(v, 0)-8) == 0


def test_additive_component_and_relative_perturbation_are_different():
    m2, r, beta1, p00 = sp.symbols("M2 r beta1 P00", positive=True)
    relative = m2*r*p00/beta1
    additive = r*r*relative
    kappa = m2*m2*r**3/(4*beta1)
    assert sp.factor(additive-4*kappa*p00/m2) == 0
    data = _independent_profile()
    c = data["c"]
    assert sp.factor(data["f2_jets"][0]-4*(c-2)) == 0
    assert sp.factor(data["f2_jets"][0]/data["r2_jets"][0]-(c-2)) == 0
    relative_second = (data["f2_jets"][1]/4
                       -data["f2_jets"][0]*data["r2_jets"][1]/16)
    assert sp.factor(relative_second.subs(c, 2)-16) == 0
    assert sp.factor(data["f2_jets"][1].subs(c, 2)-64) == 0


def test_all_three_center_jets_from_independent_logarithmic_series():
    data, actual = _independent_profile(), stationary.center_jets()
    sub = {data["c"]: stationary.C}
    pairs = (("r2_jets", "r_squared_u0_u2_u4"),
             ("f2_jets", "f2_00_u0_u2_u4"),
             ("spatial_jets", "f2_space_over_a2_u0_u2_u4"),
             ("remainder", "metric00_remainder_u0_u2_u4"))
    for independent_key, primary_key in pairs:
        for left, right in zip(data[independent_key], actual[primary_key], strict=True):
            assert sp.factor(left.subs(sub)-right) == 0


def test_full_profile_api_against_different_polynomial_representation():
    data, primary = _independent_profile(), stationary.profile()
    sub = {data["c"]: primary["c"], data["v"]: primary["v"]}
    for left, right in (("r3", "r_cubed"), ("ell", "log_r_u_over_u"),
                        ("kap", "kappa_bar"), ("P00", "P00_bar"),
                        ("f2", "f2_00"), ("spatial", "f2_ii_over_a_squared")):
        assert sp.factor(data[left].subs(sub)-primary[right]) == 0


def test_original_clock_speed_jets_from_canonical_matter_formula():
    data = _independent_profile()
    v, c = data["v"], data["c"]
    d = 1+v
    kbar = 8*(8/(c*d**12)-1)*(1-v)/d**2-1/(100*d**12)
    assert sp.factor(kbar.subs(v, 0)-data["k0"]) == 0
    assert sp.factor(2*sp.diff(kbar, v).subs(v, 0)-data["k2"]) == 0


def test_nonlinear_clock_inverse_requires_the_k2_term():
    z, k0 = sp.symbols("varphi k0", positive=True)
    k2, f0, f2, f4 = sp.symbols("k2 f0 f2 f4", real=True)
    inverse_u = z/sp.sqrt(k0)-k2*z**3/(12*k0**sp.Rational(5, 2))
    field = f0+f2*inverse_u**2/2+f4*inverse_u**4/24
    fourth = 24*sp.expand(field).coeff(z, 4)
    assert sp.factor(fourth-f4/k0**2+2*f2*k2/k0**3) == 0
    assert sp.factor(fourth-f4/k0**2) != 0
    primary = stationary.center_jets()
    c = stationary.C
    r0, r2, r4 = primary["metric00_remainder_u0_u2_u4"]
    substitution = {k0: primary["kbar_center"], k2: primary["kbar_u2_center"],
                    f0: r0, f2: r2, f4: r4}
    assert sp.factor(fourth.subs(substitution)-primary["remainder00_varphi4"]) == 0
    assert sp.factor((fourth-f4/k0**2).subs(substitution).subs(c, 3)) != 0


def test_clock_label_derivative_is_not_a_tensor_coordinate_component():
    data = _independent_profile()
    c = data["c"]
    # f_TT(varphi) retains the physical-T basis; f_varphi,varphi gains tau²/kbar.
    tau = sp.Symbol("tau", positive=True)
    component = data["f2_jets"][0]
    transformed_tensor_component = tau*tau*component/data["k0"]
    assert sp.factor(transformed_tensor_component-component) != 0
    assert sp.factor((transformed_tensor_component-component).subs({tau: 1, c: 3})) != 0


def test_continuous_fourth_jet_lower_bound_has_positive_delta_polynomial():
    data = _independent_profile()
    c, delta = data["c"], sp.Symbol("delta", positive=True)
    cleared = sp.Poly(sp.cancel(c*c*(data["remainder"][2]-10752)).subs(c, 2+delta), delta)
    assert cleared.as_dict() == {(1,): 22272, (2,): 132288,
                                 (3,): 87552, (4,): 26832}
    assert all(value > 0 for value in cleared.coeffs())


def test_small_delta_first_correction_second_jet_exceeds_sixty():
    data = _independent_profile()
    c, delta = data["c"], sp.Symbol("delta", nonnegative=True)
    cleared = sp.cancel(c*(data["f2_jets"][1]-60)).subs(c, 2+delta).expand()
    assert cleared == 8-284*delta-184*delta**2-108*delta**3
    lower = cleared.subs(delta, sp.Rational(1, 100))
    assert lower > 0
    assert lower == stationary.center_jets()["f2_u2_above_60_for_delta_le_1over100_margin"]


def test_continuous_clock_bound_uses_three_separate_positive_factors():
    data = _independent_profile()
    c = data["c"]
    delta, eta = sp.symbols("delta eta", nonnegative=True)
    r2 = data["remainder"][1]
    normalized = sp.cancel(c*r2/(4*(c-2))).subs(c, 2+delta).expand()
    assert normalized == 27*delta**2+46*delta+64
    assert sp.factor(sp.Rational(2399, 100)-data["k0"]-32*(c-2)/c) == 0
    assert sp.factor(data["k0"]-sp.Rational(799, 100)-16*(4-c)/c) == 0
    assert sp.factor((-25*c*data["k2"]).subs(c, 4-eta)-43176-1206*eta) == 0
    # On2<c<=4: R2>0,k2<0,0<k0<2399/100, and R4>10752.
    lower = sp.Rational(10752)/(sp.Rational(2399, 100)**2)
    assert lower == sp.Rational(107520000, 5755201)
    assert stationary.center_jets()["clock4_uniform_lower_bound"] == lower


def test_nonzero_limits_and_small_lower_jet_controls():
    data = _independent_profile()
    c = data["c"]
    r0, r2, r4 = data["remainder"]
    assert sp.limit(data["f2_jets"][0], c, 2, dir="+") == 0
    assert sp.limit(data["f2_jets"][1], c, 2, dir="+") == 64
    assert sp.limit(data["f2_jets"][1]/data["k0"], c, 2, dir="+") == sp.Rational(6400, 2399)
    assert sp.limit(r0, c, 2, dir="+") == 0
    assert sp.limit(r2, c, 2, dir="+") == 0
    assert sp.limit(r4, c, 2, dir="+") == 10752
    assert sp.limit(r4/data["k0"]**2-2*r2*data["k2"]/data["k0"]**3,
                    c, 2, dir="+") == sp.Rational(107520000, 5755201)


def test_inner_center_remainder_polynomial_without_claiming_fixed_window_error():
    data = _independent_profile()
    c, delta, x = data["c"], sp.Symbol("delta", positive=True), sp.Symbol("x", real=True)
    r0, r2, r4 = data["remainder"]
    polynomial = r0+r2*delta*x*x/2+r4*delta**2*x**4/24
    leading = sp.limit(polynomial.subs(c, 2+delta)/delta**2, delta, 0, dir="+")
    assert leading == 1+64*x*x+448*x**4
    assert sp.factor(leading-1) != 0


def test_physical_units_of_fixed_component_clock_jets():
    m, tau = sp.symbols("M tau", positive=True)
    # phi=M varphi and T=tau u, without changing the tensor component basis.
    f2_phi2 = sp.Rational(6400, 2399)/m**2
    remainder_T4 = 10752/tau**4
    assert f2_phi2*m*m == sp.Rational(6400, 2399)
    assert remainder_T4*tau**4 == 10752


def test_primary_scope_does_not_promote_component_nonuniformity_to_EFT_exclusion():
    calibration = stationary.calibration()
    assert "No C2-small first correction or C4-small metric remainder" in calibration["remainder_scope"]
    assert "no EFT-wide exclusion" in calibration["remainder_scope"]
    assert "formal local series only" in calibration["stationary_branch_scope"]
    assert "physical g" in calibration["matter_scope"]
