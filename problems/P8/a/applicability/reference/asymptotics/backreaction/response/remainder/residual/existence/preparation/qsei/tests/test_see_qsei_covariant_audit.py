"""Independent two-frequency, physical-source and short-index controls."""

import sympy as sp


def test_linear_potential_born_split_with_both_scattering_directions():
    x, k = sp.symbols("x k", positive=True)
    h = sp.Symbol("h", real=True)
    phase = sp.exp(-2*sp.I*k*x)
    aa = x**2/2
    bb = phase*(-x/(2*sp.I*k)+1/(4*k**2))-1/(4*k**2)
    assert sp.simplify(sp.diff(bb, x)-phase*x) == 0
    first_mode = -sp.exp(-sp.I*k*x)*(sp.exp(2*sp.I*k*x)*bb-aa)/(2*sp.I*k)
    direct = sp.diff(first_mode, x)-h*first_mode
    split = (sp.exp(-sp.I*k*x)*(-sp.Rational(1, 2)-h/(2*sp.I*k))*aa
             +sp.exp(sp.I*k*x)*(-sp.Rational(1, 2)+h/(2*sp.I*k))*bb)
    assert sp.simplify(direct-split) == 0
    history = (1-phase)/(2*sp.I*k)
    assert sp.simplify(bb-(-phase*x+history)/(2*sp.I*k)) == 0
    assert sp.simplify(split-sp.exp(-sp.I*k*x)*(-sp.Rational(1, 2)-h/(2*sp.I*k))*aa) != 0


def test_second_parseval_uses_real_history_not_complex_mode_product():
    k = sp.Symbol("k", nonnegative=True)
    # Real history exp(-s), s>=0, has transform 1/(1+2ik).
    real_integral = sp.integrate(1/(1+4*k**2), (k, 0, sp.oo))
    assert real_integral == sp.pi/4  # (pi/2) times its L2 norm squared.
    # The same half-line equality fails for complex exp((-1+i)s).
    complex_integral = sp.integrate(1/(1+(2*k-1)**2), (k, 0, sp.oo))
    assert complex_integral == 3*sp.pi/8
    assert complex_integral > real_integral


def test_exact_solution_quantum_credit_subtracts_ordinary_radiation():
    x = sp.Symbol("x", real=True)
    coupling, hbar, scale, length, delta = sp.symbols("kappa hbar A eta_star delta", positive=True)
    a = sp.Function("a", positive=True)(x)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    hubble = h/(scale*length**2*a)
    proper = lambda value: sp.diff(value, x)/(scale*length**2*a)
    total_eed = -3*(proper(hubble)+hubble**2)/coupling
    radiation_eed = 3/(coupling*scale**2*length**4*a**4)
    actual_quantum = total_eed-radiation_eed
    expected = hbar*(a**2*(u+h**2)-1)/(960*sp.pi**2*delta*scale**4*length**8*a**4)
    restored = actual_quantum.subs(coupling, 2880*sp.pi**2*delta*scale**2*length**4/hbar)
    assert sp.simplify(restored-expected) == 0
    assert sp.simplify((total_eed-actual_quantum)-radiation_eed) == 0
    assert radiation_eed != 0


def test_exact_reference_margin_does_not_borrow_old_stress_sign():
    delta, length, distance = sp.Rational(1, 10**14), sp.Rational(1, 10**10), sp.Rational(72, 10**9)
    b = sp.Rational(2, 10**12)
    loss = (9*length+sp.Rational(9, 2)*length**2+(3*b+sp.Rational(3, 4))*length**3)*distance
    assert loss < delta/54
    assert (delta/54)/(960*delta*81) == sp.Rational(1, 4199040)


def test_short_window_cannot_meet_normal_index_threshold():
    span = sp.Rational(3, 2*10**10)
    lower = (3-span**2/8)/span
    assert lower > sp.Rational(3, 4)
    assert sp.Rational(360, 10**14)/span**2 == 160000000
    # Normal comoving Jacobi volume cannot vanish while 2<=a<=3.
    assert sp.Rational(2, 3)**3 == sp.Rational(8, 27)
