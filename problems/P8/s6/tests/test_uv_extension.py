import pytest
import sympy as sp
from p8_s5 import nonlinear_d
from p8_uv import extension as e


@pytest.mark.parametrize("point,expected", [(-1, 0), (0, 0), (sp.Rational(1, 4), 0),
                                          (sp.Rational(1, 2), sp.Rational(1, 2)),
                                          (sp.Rational(3, 4), 1), (sp.Rational(9, 10), 1),
                                          (1, 1), (sp.Rational(11, 10), 1), (2, 1)])
def test_switch_exact_samples_and_endpoints(point, expected):
    assert sp.simplify(e.clock_switch(point)-expected) == 0


def test_all_extension_residuals():
    assert len(e.checks()) >= 25
    assert set(e.checks().values()) == {0}


def test_jet_preservation_through_quartic():
    # The open-branch identity proves arbitrary order. This is a finite regression.
    old = nonlinear_d.functions()
    for name, key in (("F", "F"), ("A3", "a3"), ("A4", "a4"), ("A5", "a5")):
        difference = e.branch(1)[name]-old[key].subs(nonlinear_d.X, e.X)
        for n in range(5):
            assert sp.cancel(sp.diff(difference, e.X, n).subs(e.X, 1)) == 0


def test_both_vacuum_signs_have_same_quadratic_health():
    fv = e.branch(0)["F"]
    for sign in (-1, 1):
        at = {e.X: 0, e.u: e.u0, e.coupling: sign}
        assert sp.diff(fv, e.X).subs(at) == e.Z/2
        assert sp.diff(fv, e.u, 2).subs(at) == -e.Z*e.mass**2
        assert sp.diff(fv, e.X, 2).subs(at) == 2*sign*e.Z**2


def test_dependents_must_not_be_linearly_spliced():
    old_a3 = nonlinear_d.functions()["a3"]
    wrong_a4 = e.switch*(-old_a3+e.X**2*old_a3**2/4)
    defect = sp.factor(e.family()["A4"]-wrong_a4)
    assert sp.cancel(defect-e.X**2*old_a3**2*(e.switch**2-e.switch)/4) == 0
    assert defect.subs({e.switch: sp.Rational(1, 2), e.u: 0, e.X: sp.Rational(1, 2)}) != 0


def test_invalid_branch_and_derivative_count_rejected():
    with pytest.raises(ValueError, match="constant switch"):
        e.branch(sp.Rational(1, 2))
    with pytest.raises(ValueError, match="nonnegative integer"):
        e.flat_polynomials(-1)
