"""Joint nonlinear constraints, local operator inverse and independent solves."""
import pytest
import sympy as sp
from p8_affine_nonlinear import adm, constraints, trace
from p8_affine_retuned import geometry


def test_exact_joint_legendre_and_constraint_rank_checks():
    for value in trace.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0
    assert all(value is True for value in constraints.proof_checks().values())
    assert constraints.checks()["four_auxiliary_constraint_block_inverse"] == sp.zeros(4)


@pytest.mark.parametrize("delta", (-sp.Rational(1, 10), 0, sp.Rational(1, 10)))
@pytest.mark.parametrize("volume", (1, 3))
def test_independent_simultaneous_trace_and_temporal_solve(delta, volume):
    d = trace.legendre()
    values = {d["a"]: -sp.Rational(volume, 3), d["volume"]: volume, d["gamma"]: 1,
              d["delta"]: delta, d["b"]: 2, d["f"]: 3, d["c"]: -1}
    L = d["L"].subs(values)
    K, T, p, j = [d[key] for key in ("K", "T", "p", "j")]
    answer = sp.solve((sp.diff(L, K)-p, sp.diff(L, T)+j), (K, T))
    assert len(answer) == 2
    assert sp.factor(answer[K]-d["K_joint"].subs(values)) == 0
    assert sp.factor(answer[T]-d["T_joint"].subs(values)) == 0
    H = (p*K-L-T*j).subs(answer, simultaneous=True)
    assert sp.factor(H-d["H_after_temporal"].subs(values)) == 0


@pytest.mark.parametrize("point", (-sp.Rational(1, 2), -sp.Rational(1, 4), 0, sp.Rational(1, 4), sp.Rational(1, 2)))
def test_complete_background_auxiliary_Jacobian_is_nonsingular(point):
    D = trace.clock_auxiliary()["jacobian"].subs(adm.u, point)
    assert D[0, 1] == D[1, 0] == 0
    assert D[0, 0] < -sp.Rational(1, 20)
    assert D[1, 1] == -1
    assert D.det() > sp.Rational(1, 20)


@pytest.mark.parametrize("p", (sp.sqrt(sp.Rational(9, 40)), sp.Rational(1, 2), sp.sqrt(sp.Rational(11, 40))))
def test_effective_temporal_mass_keeps_the_trace_Legendre_map_regular(p):
    gamma = geometry.update()["gamma_t"].subs(geometry.P, p)
    # The independent delta envelope is conservative on the whole tube.
    effective = gamma-3*sp.Rational(1, 100)/(8*p**2)
    assert effective > sp.Rational(9, 10)
    d = trace.legendre()
    isolated = sp.diff(d["H_before_temporal"], d["T"], 2)+d["volume"]/d["gamma"]
    assert sp.factor(isolated) != 0


def test_dense_constraint_inverse_does_not_require_secondary_bracket_inverse():
    D = sp.Matrix([[2, 1], [1, 3]])
    E = sp.Matrix([[0, 7], [-7, 0]])
    zero = sp.zeros(2)
    C = sp.BlockMatrix([[zero, -D], [D.T, E]]).as_explicit()
    inverse = sp.BlockMatrix([[D.T.inv()*E*D.inv(), D.T.inv()], [-D.inv(), zero]]).as_explicit()
    assert C*inverse == sp.eye(4)
    assert C.inv() == inverse
    assert constraints.block_inverse()["does_not_require_an_inverse_of_secondary_secondary_operator"] is True


def test_block_identity_control_is_materialized_before_subtraction():
    # Do not rely on repeated block_collapse of block_diag(I,I)-Identity(4).
    block = sp.BlockMatrix([[sp.Identity(2), sp.ZeroMatrix(2, 2)],
                           [sp.ZeroMatrix(2, 2), sp.Identity(2)]])
    assert block.as_explicit()-sp.eye(4) == sp.zeros(4)
    assert constraints.block_inverse()["residual"] == sp.zeros(4)


@pytest.mark.parametrize("D", (sp.diag(0, -1), sp.diag(-1, 0), sp.zeros(2)))
def test_either_missing_auxiliary_pivot_destroys_the_four_constraint_inverse(D):
    E = sp.Matrix([[0, 7], [-7, 0]])
    C = sp.BlockMatrix([[sp.zeros(2), -D], [D.T, E]]).as_explicit()
    assert C.det() == 0
    assert C.rank() < 4
