"""Independent dense mass-update and complete-source checks."""
import pytest
import sympy as sp
from p8_affine import connection as old
from p8_affine_retuned import geometry as g


def test_all_exact_geometry_and_complement_identities():
    for value in g.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


@pytest.mark.parametrize("p", (sp.Rational(49, 100), sp.Rational(1, 2), sp.Rational(51, 100)))
def test_independent_full_dense_inverse_and_determinant(p):
    data = g.update()
    original, changed = [data[key].subs(g.P, p) for key in ("M_old", "M")]
    inverse = changed.inv(method="DM")
    assert changed*inverse == sp.eye(60)
    old_inverse = original.inv(method="DM")
    candidate = old_inverse-old_inverse*data["N"].T*data["middle"].subs(g.P, p)*data["N"]*old_inverse
    assert candidate == inverse
    assert inverse.norm(sp.oo) < 6900
    assert data["N"]*inverse*data["N"].T == data["D"].subs(g.P, p)
    assert sp.factor(changed.det(method="domain-ge")/original.det(method="domain-ge")
                     -data["determinant_ratio"].subs(g.P, p)) == 0
    assert changed.rank() == 60


def test_literal_source_centering_includes_both_linear_term_and_counterterm():
    values = {old.P: sp.Rational(49, 100), old.S: sp.Rational(5, 4),
              old.PP: sp.Rational(2, 7), old.PX: -sp.Rational(3, 11),
              old.CP: sp.Rational(5, 13), old.CX: sp.Rational(7, 17), old.F3: -sp.Rational(1, 19)}
    values.update({item: sp.Rational(2*i-7, 23) for i, item in enumerate(old.H_SYMBOLS)})
    centered = g.source_centering()
    original = old.quadratic()
    k = sp.Matrix(old.KAPPA)
    n = g.update()["N_full"]
    tstar = centered["Tstar"].subs(values)
    delta = n*k-tstar
    addition = g.MU*(delta.T*g.ETA*delta)[0]/2
    added_matrix = sp.hessian(addition, tuple(k))
    added_source = sp.Matrix([sp.diff(addition, item) for item in k]).subs(dict.fromkeys(k, 0))
    assert added_matrix+original["hessian"].subs(values) == centered["M_full"].subs(values)
    assert added_source+original["source"].subs(values) == centered["source_new"].subs(values)
    assert addition.subs(dict(zip(k, centered["stationary"].subs(values), strict=True))) == 0
    assert addition.subs(dict.fromkeys(k, 0)) == centered["counterterm"].subs(values)
    assert centered["counterterm"].subs(values) != 0


def test_mass_flip_is_not_an_unchanged_member_of_the_closed_curl_family():
    data = g.update()
    assert data["D_old"].subs(g.P, sp.Rational(1, 2)) == -sp.Rational(9, 2)*g.ETA
    assert data["D"].subs(g.P, sp.Rational(1, 2)) == g.ETA
    assert g.MU == sp.Rational(11, 9)
    assert data["M"] != data["M_old"]
