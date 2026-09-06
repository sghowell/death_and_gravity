"""Independent physical-clock/source and differential-polynomial controls.

The algebra below does not call the main map's identity generator. Its
polynomial data are nonphysical controls, never sampled proof of existence.
"""

import sympy as sp


def test_physical_clock_and_source_weight_restore_trace_coefficient():
    length, normalization, delta, coupling, hbar = sp.symbols(
        "eta_star A delta kappa hbar", positive=True)
    a, h, u, cp = sp.symbols("a h u cprime", real=True, nonzero=True)
    a_phys, h_phys = normalization*length*a, h/length
    c_phys_prime = cp/length**5
    theta_ext = hbar*c_phys_prime/(sp.pi**2*h_phys*a_phys**4)
    # q_phys=q/(A^2 eta_star^4); physical auxiliary equation has
    # source 8pi^2 a_phys^2 Theta_ext/hbar. Restore dimensionless q''.
    source_rhs = 8*sp.pi**2*a_phys**2*theta_ext/hbar
    assert sp.simplify(source_rhs*normalization**2*length**6-8*cp/(h*a**2)) == 0
    einstein = 48*sp.pi**2*(u/length**2)/(coupling*hbar)
    assert sp.simplify((einstein*normalization**2*length**6).subs(
        coupling*hbar, 2880*sp.pi**2*delta*normalization**2*length**4)-u/(60*delta)) == 0


def test_constraint_directly_from_conserved_density_and_trace():
    x = sp.Symbol("x", real=True)
    a = 2+x+x**3/100
    rho = 1+x**2/7
    coupling = sp.Rational(13, 11)
    h = sp.diff(a, x)/a
    pressure = -rho-sp.diff(rho, x)/(3*h)
    scalar = 6*sp.diff(a, x, 2)/a**3
    density_defect = 3*sp.diff(a, x)**2/a**4-coupling*rho
    trace_defect = scalar-coupling*(rho-3*pressure)
    assert sp.cancel(sp.diff(a**4*density_defect, x)-a**4*h*trace_defect) == 0


def test_density_weight_integration_by_parts_on_nonradiation_metric():
    x = sp.Symbol("x", real=True)
    a, c = 2+x+x**4/50, 1-3*x**2+2*x**3
    h = sp.diff(a, x)/a
    u = -sp.diff(a, x, 2)/a
    assert sp.cancel(sp.diff(c/h, x)-c*(1+u/h**2)-sp.diff(c, x)/h) == 0
    assert sp.cancel(sp.diff(c/h, x)-sp.diff(c, x)/h) != 0


def test_radiation_past_is_off_shell_until_source_is_included():
    x, delta = sp.symbols("x delta", positive=True)
    a = 1+x
    h = 1/a
    quantum_density_bracket = h**4/960
    c = -quantum_density_bracket
    forced_constraint = 3*sp.diff(a, x)**2-3-2880*delta*(quantum_density_bracket+c)
    assert sp.simplify(forced_constraint) == 0
    assert sp.simplify(forced_constraint+2880*delta*c) != 0
    # A free radiation metric has u=q=0. Its conserved preparation
    # trace exactly cancels the anomaly in the auxiliary equation.
    assert sp.simplify(-h**4/(30*a**2)+8*sp.diff(c, x)/(h*a**2)) == 0


def test_flat_start_is_a_real_endpoint_hypothesis():
    x = sp.Symbol("x", positive=True)
    # Direct D_gammaE[t] calculation: -1/2 int_0^t log(t-r) dr in
    # the 4pi²-normalized convention. Its derivative diverges at zero.
    direct = -x*(sp.log(x)-1)/2
    assert sp.diff(direct, x) == -sp.log(x)/2
    assert sp.limit(sp.diff(direct, x), x, 0, dir="+") == sp.oo
