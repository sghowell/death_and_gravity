"""Independent finite matrices and full-rank sign/rank controls."""
import pytest
import sympy as sp
from p8_affine_kinetic import scalar as old
from p8_affine_traces import rank_two as r


@pytest.mark.parametrize("name", tuple(r.checks()))
def test_exact_identities(name):
    value = r.checks()[name]
    assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


@pytest.mark.parametrize("parameters", ((1, 1, 0), (2, 3, 1), (5, 1, 2), (2, 2, -1)))
def test_direct_temporal_solve_after_exact_matrix_choice(parameters):
    rr, ss, tt = parameters
    z = r.exact_matrix(rr, ss, tt)
    assert r.classify(*parameters)["chart"] == "positive_transverse_test_passed"
    bad = sp.Matrix([1, 1])
    threshold = 2*(bad.T*z.inv()*bad)[0]/9
    q = threshold+1
    velocity = sp.Matrix(sp.symbols("a b", real=True))
    time = sp.Matrix(sp.symbols("c d", real=True))
    lag = ((time.T*r.D.inv()*time)[0]+q*((velocity-time).T*z*(velocity-time))[0])/2
    solution = sp.solve([sp.diff(lag, item) for item in time], tuple(time))
    kinetic = sp.hessian(lag.subs(solution), tuple(velocity))/2
    assert kinetic.det() < 0
    assert kinetic == (r.D+z.inv()/q).inv()/2


@pytest.mark.parametrize("parameters", ((1, -1, 0), (-1, -1, 0), (1, 1, 2)))
def test_invertible_nonpositive_forms_fail_transverse(parameters):
    assert r.classify(*parameters)["chart"] == "negative_transverse_direction"


@pytest.mark.parametrize("parameters", ((1, 4, 2), (-1, -4, -2), (0, 3, 0), (0, -3, 0)))
def test_rank_one_forms_are_not_inverted_as_rank_two(parameters):
    data = r.classify(*parameters)
    assert data["rank"] == 1
    assert data["coupling"]*data["coefficient"]*data["coefficient"].T == r.exact_matrix(*parameters)


def test_zero_form_is_auxiliary():
    assert r.classify(0, 0, 0) == {"rank": 0, "chart": "unchanged_auxiliary"}


@pytest.mark.parametrize("time", (sp.Rational(1, 2), -sp.Rational(1, 2)))
def test_full_actual_metric_matter_and_temporal_vector_elimination(time):
    bg = old.background()
    mapping = {symbol: bg[name].subs(old.u, time) for symbol, name in (
        (old.theta, "theta"), (old.J, "J"), (old.lam, "lam"),
        (old.ell, "ell"), (old.w, "w"))}
    mapping[old.q] = 2
    h, hubble = bg["h"].subs(old.u, time), bg["H"].subs(old.u, time)
    d, e = 3*hubble*sp.Matrix([7, 3])/(8*h), sp.Matrix([1, -1])/h
    sigma = sp.Matrix(sp.symbols("sigma1 sigma2", real=True))
    dot = sp.Matrix(sp.symbols("dot1 dot2", real=True))
    temporal = sp.Matrix(sp.symbols("temporal1 temporal2", real=True))
    z = sp.Matrix([[2, 1], [1, 3]])
    extra = ((temporal+d*old.n).T*r.D.inv()*(temporal+d*old.n))[0]/2
    extra -= ((sigma+e*old.n).T*r.D.inv()*(sigma+e*old.n))[0]
    extra += ((dot-temporal).T*z*(dot-temporal))[0]
    lag = old.action()["base"].subs(mapping)+extra
    auxiliary = (old.n, old.shift, *temporal)
    solved = sp.solve([sp.diff(lag, item) for item in auxiliary], auxiliary)
    assert len(solved) == 4
    reduced = lag.subs(solved)
    kinetic = (sp.hessian(reduced, tuple(dot))/2).applyfunc(sp.factor)
    expected = (r.D+z.inv()/2).inv()/2
    assert kinetic == expected
    assert kinetic.det() < 0


@pytest.mark.parametrize("bad", (True, 1.0, "1", sp.Float(1), sp.I, sp.oo, sp.nan, sp.Symbol("x")))
def test_exact_rank_input_guards(bad):
    with pytest.raises((TypeError, ValueError)):
        r.classify(bad, 1, 0)
