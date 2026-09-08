"""Primary null and spatial-chart controls, not a secondary Dirac theorem."""
import pytest
import sympy as sp
from p8_affine import connection as old
from p8_affine_retuned import bounds
from p8_affine_retuned import degeneracy as d
from p8_affine_retuned import geometry as g


def test_all_exact_primary_null_and_chart_identities():
    for value in d.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


@pytest.mark.parametrize("p", (sp.sqrt(sp.Rational(9, 40)), sp.Rational(1, 2), sp.sqrt(sp.Rational(11, 40))))
def test_mass_and_trace_legendre_bounds_include_closed_endpoints(p):
    data = g.update()
    gt, gs = [sp.simplify(data[key].subs(old.P, p)) for key in ("gamma_t", "gamma_s")]
    assert sp.Rational(18, 19) < gt < sp.Rational(9, 8)
    assert sp.Rational(18, 19) < gs < sp.Rational(36, 35)
    assert 0 <= sp.simplify(d.chart()["trace_mass_relative_kinetic_addition"].subs(old.P, p)) < sp.Rational(19, 1080)
    bounds.require_domain(p, sp.Rational(1, 2000), 1, 1)


def test_source_centered_square_preserves_the_original_null_for_arbitrary_temporal_vector():
    temporal, coefficient = sp.symbols("T_normal positive_mass", real=True)
    data = d.trace_kinetic()
    action = data["action"]+coefficient*(temporal-data["Tstar_kinetic"])**2/2
    hessian = sp.hessian(action, (d.H0, d.T))
    assert (hessian*data["null"]).applyfunc(sp.factor) == sp.zeros(2, 1)
    assert d.chart()["Tstar_normal"].has(old.F3)
    assert not d.chart()["Tstar_normal"].has(d.H0)
    assert d.chart()["nonlinear_secondary_constraint_rank_claim"] is False


def test_missing_curvature_time_boundary_breaks_the_common_primary_null():
    data = d.trace_kinetic()
    wrong = data["action"]-2*d.H0*d.T/d.h
    wrong_hessian = sp.hessian(wrong, (d.H0, d.T))
    # Omitting this boundary flips the old null, not its rank by itself.
    assert sp.factor(wrong_hessian.det()) == 0
    assert (wrong_hessian*data["null"]).applyfunc(sp.factor) != sp.zeros(2, 1)
    wrong_updated = wrong+data["Tstar_kinetic"]**2/2
    assert sp.factor(sp.hessian(wrong_updated, (d.H0, d.T)).det()) != 0


def test_continuous_bounds_are_all_exact_and_true():
    assert all(value is True for value in bounds.proof_checks().values())


@pytest.mark.parametrize("p", (sp.sqrt(sp.Rational(9, 40)), sp.Rational(1, 2), sp.sqrt(sp.Rational(11, 40))))
def test_all_ten_unreduced_velocities_have_the_stated_rank(p):
    data = d.velocity_rank()
    matrix = data["Hessian"].subs({old.P: p, data["zeta"]: sp.Rational(1, 2000)})
    assert matrix.rank() == 10
    assert sp.simplify(matrix[0, 0]) < 0
    assert all(sp.simplify(matrix[1:i, 1:i].det()) > 0 for i in range(2, 11))
    assert data["physical_inertia_not_inferred"] is True
