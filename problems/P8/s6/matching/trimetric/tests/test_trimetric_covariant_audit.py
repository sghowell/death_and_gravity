"""Independent full-matrix, matter-frame and vacuum-extension regressions.

Noncommuting fixtures test the ambient vierbein variation algebra; they
are not asserted to solve the Lorentz constraints or the field equations.
"""

import pytest
import sympy as sp
from p8_trimetric import matter, model


def zero(matrix):
    assert all(sp.cancel(value) == 0 for value in matrix)


def pair(left, right):
    return sp.trace(left.T*right)


def frames():
    e = sp.eye(4)
    e[0, 1], e[2, 3] = sp.Rational(1, 3), sp.Rational(1, 5)
    v = sp.diag(2, 3, 1, 4)
    v[1, 0], v[3, 2] = sp.Rational(1, 7), sp.Rational(1, 2)
    u = sp.eye(4)
    u[0, 2], u[1, 3], u[2, 1] = sp.Rational(1, 4), sp.Rational(1, 6), sp.Rational(1, 8)
    return e, v, u


@pytest.mark.parametrize("index,key", [(0, "E_e"), (1, "E_v"), (2, "E_u")])
def test_literal_noncommuting_directional_variations(index, key):
    values = list(frames())
    direction = sp.Matrix([[1, 2, 0, -1], [0, -1, 3, 1], [2, 0, 1, 0], [1, -2, 0, 2]])
    t = sp.Symbol("variation", real=True)
    changed = list(values)
    changed[index] += t*direction
    b, pg, pf, bg, bf = -9, 2, 3, -1, 5
    e, v, u = changed
    # Reconstruct the density directly from its determinant and trace,
    # without calling the primary potential function.
    literal = -2*u.det()*(b+sp.trace(u.inv()*(pg*e+pf*v)))-2*bg*e.det()-2*bf*v.det()
    expected = sp.diff(literal, t).subs(t, 0)
    actual = model.euler_maps(*values, b=b, pg=pg, pf=pf, bg=bg, bf=bf, epsilon=0)
    assert sp.cancel(expected-pair(actual[key], direction)) == 0


def test_einstein_density_coframe_dictionary_in_both_sectors():
    e, v, u = frames()
    direction = sp.Matrix([[0, 1, 2, 0], [3, 0, 0, 1], [0, 1, -1, 0], [2, 0, 0, 1]])
    einstein = sp.Matrix([[2, 1, 0, 0], [1, 3, 1, 0], [0, 1, -1, 2], [0, 0, 2, 4]])
    actual = model.euler_maps(e, v, u, b=0, pg=0, pf=0, g=7, f=11,
                              einstein_g=einstein, einstein_f=einstein, epsilon=0)
    for frame, coefficient, key in ((e, 7, "E_e"), (v, 11, "E_v")):
        metric_variation = direction.T*model.ETA*frame+frame.T*model.ETA*direction
        expected = coefficient*frame.det()*pair(einstein, metric_variation)/2
        assert sp.cancel(expected-pair(actual[key], direction)) == 0


def test_arbitrary_source_trace_reversal_is_full_matrix_identity():
    e, v, u = frames()
    ju = sp.Matrix([[1, 2, 0, 3], [0, 1, -2, 0], [3, 0, 1, 2], [1, 0, 0, -1]])
    epsilon = sp.Rational(2, 7)
    actual = model.euler_maps(e, v, u, b=-5, pg=2, pf=3, matter_gradient=ju, epsilon=epsilon)
    tmat = u*ju.T*u/(2*u.det())
    calt = sp.trace(u.inv()*tmat)*u/3-tmat
    c = 2*e+3*v-sp.Rational(5, 3)*u-epsilon*calt
    reconstructed = 2*u.det()*u.inv().T*(c.T-sp.trace(u.inv()*c)*u.T)*u.inv().T
    zero(actual["E_u"]-reconstructed)


def test_nonzero_link_gradient_cannot_vanish_for_any_auxiliary_source():
    e, v, u = frames()
    ju = sp.ones(4)*17
    actual = model.euler_maps(e, v, u, b=1001, pg=2, pf=-3, matter_gradient=ju, epsilon=19)
    for link, key in ((2, "E_e"), (-3, "E_v")):
        zero(u.T*actual[key]+2*link*u.det()*sp.eye(4))
        assert sp.cancel(actual[key].det()-16*link**4*u.det()**3) == 0
        assert actual[key].det() != 0


def test_one_disconnected_link_does_not_remove_the_other_obstruction():
    e, v, u = frames()
    actual = model.euler_maps(e, v, u, pg=0, pf=2, epsilon=0)
    zero(actual["E_e"])
    assert actual["E_v"].det() != 0


@pytest.mark.parametrize("potential", [0, 7])
def test_separately_named_cosmological_extension_has_exact_flat_vacuum(potential):
    ident = sp.eye(4)
    q, epsilon = sp.Integer(3), sp.Rational(2, 5)
    ju = -potential*ident
    actual = model.euler_maps(ident, ident, ident, pg=q, pf=q,
                              b=-6*q-epsilon*potential/2, bg=-q, bf=-q,
                              matter_gradient=ju, epsilon=epsilon)
    for key in ("E_e", "E_v", "E_u"):
        zero(actual[key])


def test_canonical_matter_gradient_in_nonorthogonal_frame():
    _, _, u = frames()
    k = sp.Matrix([3, 1, -1, 0])
    direction = sp.Matrix([[1, 1, 0, 2], [0, 1, 2, 0], [1, 0, -1, 1], [0, 2, 0, 1]])
    t = sp.Symbol("variation", real=True)
    changed = u+t*direction
    h = changed.T*model.ETA*changed
    literal = changed.det()*((k.T*h.inv()*k)[0]/2-5)
    expected = sp.diff(literal, t).subs(t, 0)
    actual = matter.canonical_source(u, k, 5)
    assert sp.cancel(expected-pair(actual["J_u"], direction)) == 0


def test_actual_physical_metric_has_a_nonzero_first_order_correction():
    e, v = sp.diag(1, 2, 3, 4), sp.eye(4)
    k = sp.Matrix([1, 1, 0, 0])
    actual = matter.first_order(e, v, k, b=-6, pg=1, pf=1)
    expected = -(3*k*k.T-actual["Y"]*actual["h0"]/2)/6
    zero(actual["h1"]-expected)
    assert actual["Y"] == sp.Rational(5, 9)
    assert actual["h1"][0, 1] == -sp.Rational(1, 2)
    anti = e.T*model.ETA*actual["calT"]-actual["calT"].T*model.ETA*e
    assert anti[0, 1] != 0


def test_constant_potential_shift_is_exact_before_elimination():
    e, v, u = frames()
    potential, epsilon = sp.Integer(7), sp.Rational(2, 3)
    source = -potential*u.det()*u.inv().T
    first = model.euler_maps(e, v, u, b=-8, pg=2, pf=3,
                             matter_gradient=source, epsilon=epsilon)
    second = model.euler_maps(e, v, u, b=-8+epsilon*potential/2,
                              pg=2, pf=3, epsilon=0)
    for key in ("E_e", "E_v", "E_u"):
        zero(first[key]-second[key])
