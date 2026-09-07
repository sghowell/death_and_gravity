"""Independent physical variables, source measure and retarded-response audits."""

from fractions import Fraction as F
from itertools import product

import sympy as sp
from p8_composite_response import model
from p8_composite_response.intervals import Interval, Jet


def test_unscaled_two_metric_kinetics_gradients_and_stiffness_reproduce_flux_action():
    m2, sigma, eps, z, kap, ae, q = sp.symbols(
        "M_squared sigma epsilon z kappa A_e q", positive=True)
    e = eps*z
    k1 = m2*ae**3*e**2*(e+kap)/(8*(1+e)**3)
    k2 = m2*ae**3*(e+kap)/(8*kap*(1+e)**3)
    kr = k1*k2/(k1+k2)
    cg, cf = (1+e)**2/(e+kap)**2, kap**2*(1+e)**2/(e+kap)**2
    mass = (1-e)*(1-kap)*(1+kap*e**2)/((1+e)*(e+kap)**2)
    scale = m2*eps**2*sigma**2/8
    derived = {
        "A_g": k1*sigma**2/scale,
        "A_f": k2*eps**2*sigma**2/scale,
        "V_g": k1*sigma**2*cg*q/(ae**2*scale),
        "V_f": k2*eps**2*sigma**2*cf*q/(ae**2*scale),
        "W": kr*mass*sigma**2/scale,
    }
    expected = {
        "A_g": ae**3*z**2*(e+kap)/(1+e)**3,
        "A_f": ae**3*(e+kap)/(kap*(1+e)**3),
        "V_g": ae*z**2*q/((1+e)*(e+kap)),
        "V_f": ae*kap*q/((1+e)*(e+kap)),
        "W": ae**3*z**2*(1-e)*(1-kap)/((1+e)**4*(e+kap)),
    }
    d = model.derive()
    mapping = {eps: d["epsilon"], z: d["z"], kap: d["kappa"], ae: d["A_e"], q: d["q"]}
    for name, value in derived.items():
        assert sp.cancel(value-expected[name]) == 0
        assert sp.cancel(expected[name].subs(mapping, simultaneous=True)-d[name]) == 0


def test_composite_projection_and_physical_source_measure_include_both_rescalings():
    eps, z, sigma, m2, p, g, f = sp.symbols("epsilon z sigma M_squared p g f", real=True)
    e = eps*z
    physical = (e*sigma*g+eps*sigma*f)/(1+e)
    source = m2*eps*sigma*p/4
    action_scale = m2*eps**2*sigma**2/8
    assert sp.cancel(source*physical/action_scale-2*p*(z*g+f)/(1+e)) == 0
    d = model.derive()
    assert sp.cancel(d["output_g"]-d["z"]/(1+d["epsilon"]*d["z"])) == 0
    assert sp.cancel(d["output_f"]-1/(1+d["epsilon"]*d["z"])) == 0
    # Proper TT stress sign in one polarization: delta g_ij=-Ae² Gamma_ij.
    ae, mass, probe = sp.symbols("A_e m Pi", positive=True)
    normalized_u_source = -ae**3*probe/(2*mass)
    assert sp.cancel(-2*mass*normalized_u_source/ae**3-probe) == 0
    # The displayed model divides the literal u-action by m. Its source
    # Jbar therefore differs from the literal J_u by exactly this m.
    jbar = normalized_u_source/mass
    assert sp.cancel(-2*mass**2*jbar/ae**3-probe) == 0


def test_external_TT_probe_is_covariantly_conserved_without_a_time_profile_equation():
    t, x = sp.symbols("t x", real=True)
    scale = sp.Function("a", positive=True)(t)
    pulse = sp.Function("Pi")(t)
    momentum = sp.Symbol("k", real=True)
    metric = sp.diag(1, -scale**2, -scale**2, -scale**2)
    inverse = metric.inv()

    def derivative(expr, index):
        return sp.diff(expr, t) if index == 0 else sp.diff(expr, x) if index == 1 else sp.S.Zero

    connection = [[[sum(inverse[i, n]*(derivative(metric[n, j], k)
                        +derivative(metric[n, k], j)-derivative(metric[j, k], n))/2
                        for n in range(4))
                    for k in range(4)] for j in range(4)] for i in range(4)]
    stress = sp.zeros(4)
    stress[2, 2] = scale**2*pulse*sp.exp(sp.I*momentum*x)
    stress[3, 3] = -stress[2, 2]
    assert sp.trace(inverse*stress) == 0
    for b in range(4):
        divergence = sum(inverse[a, c]*(derivative(stress[a, b], c)
            -sum(connection[d][c][a]*stress[d, b]+connection[d][c][b]*stress[a, d]
                 for d in range(4))) for a, c in product(range(4), repeat=2))
        assert sp.simplify(divergence) == 0


def test_literal_locked_projection_has_no_residual_source_or_relative_potential():
    d = model.derive()
    eps, g, f, gp, fp, p = sp.symbols("eps g f gp fp p", real=True)
    lagrangian = d["A_g"]*gp**2+d["A_f"]*fp**2-d["V_g"]*g*g-d["V_f"]*f*f
    lagrangian -= d["W"]*(g-d["epsilon"]*f)**2
    lagrangian += 2*p*(d["output_g"]*g+d["output_f"]*f)
    ell, velocity = sp.symbols("ell velocity", real=True)
    locked = lagrangian.subs({g: d["epsilon"]*ell, f: ell,
        gp: d["epsilon"]*velocity, fp: velocity}, simultaneous=True)
    assert sp.cancel(locked-(d["C_lock"]*velocity**2-d["V_lock"]*ell**2+2*p*ell)) == 0
    assert sp.cancel(d["C_lock"]-d["A_f"]-d["epsilon"]**2*d["A_g"]) == 0
    assert sp.cancel(d["V_lock"]-d["V_f"]-d["epsilon"]**2*d["V_g"]) == 0
    assert eps not in locked.free_symbols


def test_flux_matrix_is_the_full_forced_Euler_Lagrange_system():
    d = model.derive()
    g, pg, f, pf, ell, pl, p = sp.symbols("g Pg f Pf L PL pulse", real=True)
    state = sp.Matrix([g, pg, f, pf, ell, pl])
    rhs = sp.Matrix(d["matrix"])*state+p*sp.Matrix(d["forcing"])
    expected = sp.Matrix([pg/d["A_g"],
        -d["V_g"]*g-d["W"]*(g-d["epsilon"]*f)+p*d["output_g"],
        pf/d["A_f"], -d["V_f"]*f+d["epsilon"]*d["W"]*(g-d["epsilon"]*f)+p*d["output_f"],
        pl/d["C_lock"], -d["V_lock"]*ell+p])
    assert all(sp.cancel(value) == 0 for value in rhs-expected)


def test_fixed_old_density_requires_parameter_dependent_initial_positive_root():
    eps = sp.Symbol("epsilon", nonnegative=True)
    root = sp.sqrt(1+3*(1+eps)**3)
    assert sp.cancel(3*(root*root-1)/(1+eps)**3-9) == 0
    assert sp.diff(root, eps).subs(eps, 0) == sp.Rational(9, 4)
    assert sp.cancel(3*(sp.Integer(2)**2-1)/(1+eps)**3-9) != 0


def test_full_root_variable_flow_preserves_original_density_and_positive_null_source():
    d = model.derive()
    e, eta, j, lam = (d[key] for key in ("e", "eta", "j", "lambda"))
    density_flow = -3*j*(eta*e+2/(1+e)**2)+2*eta*lam
    explicit_eta_flow = 6*d["R"]*d["Rprime"]/(1+e)**3
    explicit_eta_flow += 9*(d["R"]**2-1)*e*lam/(1+e)**4
    assert sp.cancel(explicit_eta_flow-density_flow) == 0
    assert sp.cancel(d["null_over_e"]-(eta*e+2/(1+e)**2)) == 0


def test_limiting_potential_green_comparison_has_a_uniform_nonzero_momentum_margin():
    assert F(1, 4)-F(140, 19)*F(1, 8)**2 == F(41, 304) > F(1, 8)
    u, s = sp.symbols("u s", real=True)
    kernel = sp.sqrt(8)*sp.sinh((u-s)/sp.sqrt(8))
    assert sp.simplify(sp.diff(kernel, u, 2)-kernel/8) == 0
    assert kernel.subs(u, s) == 0
    assert sp.diff(kernel, u).subs(u, s) == 1


def test_positive_pulse_endpoint_weight_and_finite_error_gate():
    s = sp.Symbol("s", real=True)
    pulse = s*s*(1-s)**2
    for endpoint in (0, 1):
        assert pulse.subs(s, endpoint) == sp.diff(pulse, s).subs(s, endpoint) == 0
    assert sp.integrate((1-s)*pulse, (s, 0, 1)) == sp.Rational(1, 60)
    # R0=L0+H0, H0>=I, 0<=L0<=I. Errors in R,L each <=I/10.
    # (2/3)R-L >= (2/3)H0-(1/3)L0-(5/3)err >= I/3-(5/3)err.
    # No false upper bound on the potentially large heavy response is used.
    i = F(1, 60)
    err = F(1, 600)
    assert i/F(3)-F(5, 3)*err == i/6 > 0


def test_small_individual_metric_amplitudes_do_not_mean_unsuppressed_absolute_output():
    eps, sigma, z, g, f = sp.symbols("epsilon sigma z g f", positive=True)
    metric_g, metric_f = sigma*g, eps*sigma*f
    composite = (eps*z*metric_g+metric_f)/(1+eps*z)
    assert sp.limit(composite, eps, 0) == 0
    assert sp.cancel(sp.limit(composite/(eps*sigma), eps, 0)-z*g-f) == 0


def test_outward_radical_grid_enclosures_protect_both_endpoints():
    for value in (F(1, 3), F(2), F(1, 10**12), F(12345, 6789)):
        root = Interval(value).sqrt()
        assert root.lo**2 <= value <= root.hi**2
    for a, b in ((F(-2, 3), F(1, 7)), (F(-3), F(-1, 2)), (F(1, 11), F(7, 3))):
        interval = Interval(a, b)
        assert interval.lo <= a <= b <= interval.hi
        squared = interval*interval
        for point in (a, (a+b)/2, b):
            assert squared.contains(point*point)


def test_full_formula_AD_enclosures_contain_independent_symbolic_partials():
    d = model.derive()
    eps, u, z, root, q = (F(1, 65536), F(1, 2), F(3), F(1), F(1, 128))
    ad = model.expressions(Jet.variable(eps, 0), Jet.constant(u), Jet.variable(z, 1),
                          Jet.variable(root, 2), Jet.constant(q))
    point = {d[key]: sp.Rational(value.numerator, value.denominator)
             for key, value in zip(("epsilon", "u", "z", "R", "q"),
                                   (eps, u, z, root, q), strict=True)}
    for key in ("X_bar", "eta", "lambda", "kappa", "C_lock", "output_g"):
        value = sp.simplify(d[key].subs(point, simultaneous=True))
        assert value.is_Rational
        assert ad[key].value.contains(F(int(value.p), int(value.q)))
        for index, name in enumerate(("epsilon", "z", "R")):
            derivative = sp.simplify(sp.diff(d[key], d[name]).subs(point, simultaneous=True))
            assert derivative.is_Rational
            assert ad[key].gradient[index].contains(F(int(derivative.p), int(derivative.q)))
