"""All-component alignment, nonlinear remainder and field-map controls."""
import sympy as sp
from p8_affine_aligned import alignment as a
from p8_affine_aligned import dynamics as d
from p8_affine_retuned import degeneracy as previous


def test_exact_complete_source_alignment():
    for value in a.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_source_starts_at_second_order_but_is_not_a_nonlinear_spectator():
    data = a.rolling()
    assert data["zero_order"] == data["first_order"] == 0
    assert data["second_order"] != 0
    assert a.normal_source()["nonlinear_source_is_not_identically_zero"] is True
    assert a.quadratic_mass()["third_order_source_coupling"] != 0


def test_normal_source_keeps_the_primary_lapse_null():
    s = a.affine.S
    shifted = previous.chart()["Tstar_normal"]-s*a.coefficient().subs(a.X, -s**2)
    assert not shifted.has(previous.H0)
    assert sp.diff(shifted, previous.KHAT) == 4*a.affine.P**2-1
    assert shifted.has(a.affine.F3)


def test_temporal_field_map_has_unit_jacobian_and_preserves_joint_constraints():
    old = d.old
    shift = d.previous.DD
    before = d.scalar()["literal_before_temporal_shift"]
    after = d.scalar()["after_shift"]
    mapping = {old.temporal: old.temporal-shift*old.n}
    assert sp.factor(before.subs(mapping)-after) == 0
    jacobian = sp.Matrix([old.n, old.temporal-shift*old.n]).jacobian((old.n, old.temporal))
    assert jacobian.det() == 1
    prior = sp.Matrix([sp.diff(before, old.n), sp.diff(before, old.temporal)]).subs(mapping)
    new = sp.Matrix([sp.diff(after, old.n), sp.diff(after, old.temporal)])
    assert (new-jacobian.T*prior).applyfunc(sp.factor) == sp.zeros(2, 1)


def test_wrong_sign_or_unshifted_curl_does_not_cancel_the_source():
    old = d.old
    shift = d.previous.DD
    wrong = (old.action()["base"]+(old.temporal+shift*old.n)**2/2-old.q*old.sigma**2/2
             +d.ZETA*old.q*(old.sigmad-old.temporal+shift*old.n)**2/2)
    changed = wrong.subs(old.temporal, old.temporal-shift*old.n)
    assert sp.diff(changed-d.scalar()["after_shift"], old.sigmad, old.n) != 0
