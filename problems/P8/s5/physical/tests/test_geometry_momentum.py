import pytest
import sympy as sp
from p8_physical import geometry, momentum
from p8_physical import jets as j
from p8_physical.vertices import tensor_basis

WAVES = ((1, 0, 1), (0, 1, 1), (1, 1, -1), (-2, -2, -1))


def test_coordinate_curvature_against_independent_conformal_identity():
    c = j.Context(WAVES)
    v = sum(c.leg(i) for i in range(4))
    geo = geometry.derive(c, v, j.zeros(c))
    lap = sum(v.derivative(i).derivative(i) for i in range(3))
    grad2 = sum(v.derivative(i)**2 for i in range(3))
    target = -4*lap*(1+2*v).power(-2)+6*grad2*(1+2*v).power(-3)
    assert (geo["curvature"]-target).is_zero()
    product = j.mmul(geo["metric"], geo["inverse"])
    assert all((product[i][k]-int(i == k)).is_zero() for i in range(3) for k in range(3))


def test_generic_york_projector_not_a_one_dimensional_inverse():
    out = momentum.projector_checks()
    assert out["inverse"] == sp.zeros(3)
    assert out["determinant"] != 0


def test_all_constraints_through_third_order_with_scalar_and_tensor_legs():
    c = j.Context(WAVES)
    v, p = c.leg(0), c.leg(1)
    E, F = sp.Matrix(tensor_basis(WAVES[2])[0]), sp.Matrix(tensor_basis(WAVES[3])[1])
    tensor = [[E[i, k]*c.leg(2) for k in range(3)] for i in range(3)]
    tensor_p = [[F[i, k]*c.leg(3) for k in range(3)] for i in range(3)]
    geo = geometry.derive(c, v, tensor)
    out = momentum.derive(c, geo, v, tensor, p, tensor_p, sp.Rational(6, 5))
    assert len(out["checks"]) == 9
    assert all(out["checks"].values())
    assert any(not item.is_zero() for item in out["vector_orders"][1])
    assert any(not item.is_zero() for item in out["vector_orders"][2])


def test_zero_momentum_source_fails_instead_of_being_discarded():
    c = j.Context([(1, 0, 0), (-1, 0, 0)])
    source = [c.leg(0)*c.leg(1), c.jet(), c.jet()]
    with pytest.raises(momentum.ZeroMomentumConstraint):
        momentum.solve_vector(c, source)


def test_scalar_inputs_generate_transverse_shift_at_second_order():
    c = j.Context([(1, 2, 0), (0, 1, 3), (-1, -3, -3)])
    v, p = c.leg(0), c.leg(1)
    geo = geometry.derive(c, v, j.zeros(c))
    out = momentum.derive(c, geo, v, j.zeros(c), p, j.zeros(c), sp.Rational(6, 5))
    vector = out["vector_orders"][1]
    curl = [vector[(i+2) % 3].derivative((i+1) % 3)
            - vector[(i+1) % 3].derivative((i+2) % 3) for i in range(3)]
    assert any(not item.is_zero() for item in curl)


def test_rational_tensor_basis_transverse_tracefree_orthogonal():
    k = sp.Matrix((1, 2, 3))
    plus, cross = map(sp.Matrix, tensor_basis(tuple(k)))
    assert plus*k == cross*k == sp.zeros(3, 1)
    assert sp.trace(plus) == sp.trace(cross) == sp.trace(plus*cross) == 0


def test_flat_metric_pullback_has_zero_curvature_through_fourth_order():
    c = j.Context(WAVES)
    displacement = [c.leg(0)+2*c.leg(3), c.leg(1)-c.leg(3), c.leg(2)+c.leg(3)]
    jac = [[displacement[k].derivative(i) for k in range(3)] for i in range(3)]
    perturbation = [[jac[i][k]+jac[k][i]+sum(jac[i][a]*jac[k][a] for a in range(3))
                     for k in range(3)] for i in range(3)]
    geo = geometry.derive(c, c.jet(), perturbation)
    assert geo["curvature"].is_zero()
