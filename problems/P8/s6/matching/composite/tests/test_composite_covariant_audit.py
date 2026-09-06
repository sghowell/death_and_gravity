"""Independent positive-root metric, physical-clock and CD-shape checks.

These controls do not call the production identity generators. They are
exact background tests, never a surrogate perturbation or UV certificate.
"""

import sympy as sp


def test_covariant_positive_root_gives_both_affine_physical_components():
    a, b, ng, nf, alpha, beta = sp.symbols("a b Ng Nf alpha beta", positive=True)
    g = sp.diag(ng**2, -a**2, -a**2, -a**2)
    f = sp.diag(nf**2, -b**2, -b**2, -b**2)
    root = sp.diag(nf/ng, b/a, b/a, b/a)
    assert sp.simplify(root**2-g.inv()*f) == sp.zeros(4)
    effective = sp.simplify(g*(alpha*sp.eye(4)+beta*root)**2)
    assert effective == sp.diag((alpha*ng+beta*nf)**2,
                                -(alpha*a+beta*b)**2,
                                -(alpha*a+beta*b)**2,
                                -(alpha*a+beta*b)**2)


def test_actual_free_clock_bounce_and_no_cd_shape_at_any_time_scale():
    y = sp.Symbol("y", positive=True)
    x = sp.sqrt(y*(2+5*y)/3)
    z = -sp.sqrt((5+2*y)/(3*y**2))
    denominator = (y+5)*x-(5*y+1)*z
    speed = 6*y*(1+y)*x*z/denominator
    scale = ((1+y)**2/(4*y))**sp.Rational(1, 6)
    h = sp.factor(sp.diff(scale, y)*speed/scale)
    h1 = sp.diff(h, y)*speed
    h2 = sp.diff(h1, y)*speed
    h3 = sp.diff(h2, y)*speed
    first, third = sp.simplify(h1.subs(y, 1)), sp.simplify(h3.subs(y, 1))
    assert sp.simplify(h.subs(y, 1)) == 0
    assert first == sp.Rational(7, 36)
    assert third == -sp.Rational(73, 216)
    assert third+sp.Rational(3, 2)*first**2 == -sp.Rational(9, 32)
    # Every CD time scale has the vanishing dimensionless shape invariant.
    time, tau = sp.symbols("T tau", real=True, positive=True)
    cd = (1+time**2/tau**2)**2
    a2 = sp.diff(cd, time, 2).subs(time, 0)
    a4 = sp.diff(cd, time, 4).subs(time, 0)
    assert sp.simplify(a4-sp.Rational(3, 2)*a2**2) == 0


def test_double_pressure_root_is_not_zero_ratio_velocity():
    y = sp.Symbol("y", positive=True)
    pressure = 2*y/(1+y)**2
    at_bounce_pressure = 2*y-sp.Rational(1, 2)*(1+y)**2
    assert sp.factor(at_bounce_pressure) == -(y-1)**2/2
    assert sp.diff(pressure, y).subs(y, 1) == 0
    assert sp.diff(pressure, y, 2).subs(y, 1) == -sp.Rational(1, 4)
    assert -sp.diff(pressure, y, 2).subs(y, 1)*sp.Rational(7, 3)/(6*pressure.subs(y, 1)) == sp.Rational(7, 36)


def test_local_reconstruction_clock_without_dividing_effective_hubble():
    y, x, z, h = sp.symbols("y X Y h", real=True)
    denominator = x-y*z
    ng, nf = ((1+y)*h-y*z)/denominator, (x-(1+y)*h)/denominator
    speed = y*(1+y)*(x*z-h*(x+z))/denominator
    assert sp.cancel(ng+nf-1) == 0
    assert sp.cancel(h-speed/(1+y)-ng*x) == 0
    assert sp.cancel(h+speed/(y*(1+y))-nf*z) == 0
    at_bounce = {y: 1, x: sp.sqrt(sp.Rational(7, 3)), z: -sp.sqrt(sp.Rational(7, 3)), h: 0}
    assert ng.subs(at_bounce) == nf.subs(at_bounce) == sp.Rational(1, 2)
