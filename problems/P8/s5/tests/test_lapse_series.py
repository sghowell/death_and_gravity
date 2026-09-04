import sympy as sp
from p8_s5 import lapse_series as ls
from p8_s5 import nonlinear_d as m


def test_generic_constraint_and_stationary_expansion():
    assert all(value == 0 for value in ls.generic_checks().values())


def test_exact_invariant_degrees():
    result = ls.derive()
    for group, first in (("lapse", 1), ("hamiltonian", 0)):
        for degree, expr in enumerate(result[group], start=first):
            for key, coefficient in ls.coefficients(expr).items():
                p, r, w = map(int, key.split(","))
                assert coefficient != 0
                assert p+r+2*w == degree


def test_bounce_quadratic_lapse_and_tensor_invariant():
    result = ls.derive()
    assert sp.expand(result["lapse"][0].subs(m.u, 0)-m.rho/8) == 0
    expected = -m.sigma**2/3+2*m.shear2+m.rho**2/32
    assert sp.expand(result["hamiltonian"][2].subs(m.u, 0)-expected) == 0


def test_missing_second_order_lapse_changes_quartic_answer():
    # h=-2 n^2+eps n+eps^2 n. Exact stationary value=(eps+eps^2)^2/8.
    A = (0, 0, sp.Integer(-4), 0, 0)
    L = (0, sp.Integer(1), 0, 0)
    Q = (0, sp.Integer(1), 0)
    result = ls.stationary_series(A, L, Q)
    assert result["hamiltonian"][4] == sp.Rational(1, 8)
    assert result["hamiltonian"][4] != 0  # naive n=n1 gives zero here


def test_quartic_does_not_drop_shear_squared_channel():
    quartic = ls.coefficients(ls.derive()["hamiltonian"][4])
    assert "0,0,2" in quartic
    assert quartic["0,0,2"].subs(m.u, 0) != 0


def test_bounce_curvature_series_against_literal_nonlinear_action():
    # Independent literal bounce slice: F(0,X)=-19+20X-5X^2.
    # This is an invariant-direction regression, not a scalar Fourier amplitude.
    n, eps = sp.symbols("n eps")
    N = 1+n
    omega = (N**-4-1)/2
    F = -19+20/N**2-5/N**4
    h = -N*sp.exp(3*omega)*F-N*sp.exp(omega)*eps/2
    truncated = sp.series(h, n, 0, 5).removeO()
    n_star = eps/8+97*eps**2/128+9459*eps**3/1024
    force = sp.Poly(sp.diff(truncated, n).subs(n, n_star), eps)
    assert all(force.nth(k) == 0 for k in range(4))
    reduced = sp.Poly(truncated.subs(n, n_star), eps)
    assert reduced.nth(2) == sp.Rational(1, 32)
    assert reduced.nth(3) == sp.Rational(29, 256)
    assert reduced.nth(4) == sp.Rational(25553, 24576)
