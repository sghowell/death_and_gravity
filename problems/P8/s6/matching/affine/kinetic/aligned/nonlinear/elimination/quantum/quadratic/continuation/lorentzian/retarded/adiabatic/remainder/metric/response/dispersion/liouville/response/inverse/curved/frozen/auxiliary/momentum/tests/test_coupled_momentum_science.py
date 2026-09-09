"""Independent scientific checks without the recursively expensive parent replay."""
import json
from fractions import Fraction

import pytest
import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_momentum import (
    canonical,
    constraints,
    fixtures,
    projector,
    quadratic,
    verify,
)
from p8_physical import jets as j


def test_every_symbolic_time_and_scale_residual():
    values=verify.residuals()
    assert len(values)==147
    for value in values.values():
        assert all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0


def test_all_scientific_audit_gates():
    assert len(verify.gates())==17
    assert all(value is True for value in verify.gates().values())


@pytest.mark.parametrize("n",(2,3,4))
def test_complete_mixed_constraint_recursion(n):
    data=fixtures.build(n)
    assert len(data["result"]["vector_orders"])==n-1
    assert len(data["result"]["checks"])==3*(n-1)
    assert all(data["result"]["checks"].values())


def test_exact_three_component_York_inverse():
    data=projector.data()
    assert (data["operator"]*data["inverse"]).applyfunc(sp.cancel)==sp.eye(3)


def test_full_tensor_lift_bound_is_not_scalar_only():
    data=projector.data()
    k=sp.Matrix((1,2,3))
    mapping=dict(zip(data["wavevector"],k))
    R=data["inverse"].subs(mapping)
    L=data["real_York_lift_without_common_i"].subs(mapping)
    norm=(R.T*L.T*L*R).applyfunc(sp.factor)
    assert norm==data["full_lift_inverse_norm_squared"].subs(mapping)
    assert set(norm.eigenvals())=={sp.Rational(3,28),sp.Rational(1,7)}


def test_nonvanishing_rational_TT_frames_and_dual_pairing():
    for n in (2,3,4):
        data=fixtures.build(n)
        for k,E in zip(data["context"].momenta,data["polarizations"]):
            assert E*sp.Matrix(k)==sp.zeros(3,1)
            assert sp.trace(E)==0
            assert sp.trace(E*E)==2
            assert sp.trace(E*(E/2))==1


@pytest.mark.parametrize("name",(
    "omitted_vector_divergence_constraint_defect",
    "omitted_background_matter_linear_constraint_solution_defect"))
def test_actual_omissions_have_nonzero_defects(name):
    assert any(value!=0 for value in fixtures.negative_controls()[name])


def test_both_canonical_matter_boundary_shifts():
    d=canonical.data()
    assert d["old_curvature_momentum"]==canonical.p-3*canonical.P0*canonical.chi
    assert d["old_matter_momentum"]==canonical.pm-3*canonical.P0*canonical.v
    assert canonical.checks()["matter_boundary_one_form_identity"]==0
    assert canonical.checks()["cross_momentum_bracket_zero"]==0


def test_background_comoving_matter_momentum_not_time_frozen():
    bg=model.coefficients()["background"]
    assert sp.factor(bg["a"]**3*bg["ell"]-sp.Rational(1,10))==0
    assert sp.factor(sp.diff(bg["a"]**3*bg["ell"],model.u))==0
    assert sp.diff(bg["ell"],model.u)!=0


def test_actual_linear_lapse_force_keeps_matter_and_curvature():
    data=quadratic.lapse_coefficients()
    bg=model.coefficients()["background"]
    expected={str(model.dp):3*bg["theta"],str(model.dc):-bg["w"],
              str(model.j):0,str(model.curvature):-bg["lam"]/2}
    assert all(sp.factor(data[key]-value)==0 for key,value in expected.items())


def test_compact_scalar_Hamiltonian_exactly_replays_prior_block():
    data=quadratic.compact_scalar()
    assert data["prior_regular_scalar_replay"]==0
    assert data["J_e"].subs(model.u,0)==sp.Rational(1199,800)+sp.Rational(4,10**6)


def test_all_105_independent_quadratic_pairs():
    pairs=[quadratic.pair(left,right) for n,left in enumerate(quadratic.CHANNELS)
           for right in quadratic.CHANNELS[n:]]
    assert len(pairs)==105
    assert all(row["actual"]==row["expected"] and row["residual"]==0 for row in pairs)


def test_full_fourteen_phase_channel_adjoint():
    data=quadratic.phase_matrix()
    assert data.shape==(14,14)
    assert data==data.T.subs(constraints.wave_scale,-constraints.wave_scale)


def test_vector_mass_and_longitudinal_gauss_energy_retained():
    k=constraints.wave_scale
    assert quadratic.pair("Px","Px")["actual"]==10**6+k*k
    assert quadratic.pair("Py","Py")["actual"]==10**6
    assert quadratic.pair("Wx","Wx")["actual"]==1
    assert sp.factor(quadratic.pair("Wy","Wy")["actual"]-1-k*k/10**6)==0


def test_both_tensor_polarizations_are_canonical():
    k=constraints.wave_scale
    for field,momentum in (("tensor","tensor_p"),("tensor2","tensor2_p")):
        assert quadratic.pair(field,field)["actual"]==k*k/2
        assert quadratic.pair(momentum,momentum)["actual"]==2


def test_actual_homogeneous_constraint_source_rejected():
    ctx=constraints.context(((1,0,0),(-1,0,0)))
    source=[j.Jet(ctx,{ctx.full:ctx.one_coefficient}),ctx.jet(),ctx.jet()]
    with pytest.raises(constraints.ZeroMomentumConstraint):
        constraints.old.solve_vector(ctx,source)


def test_native_validation_before_cache():
    fixtures.build(2)
    quadratic.pair("curvature","curvature")
    for bad in (True,False,2.0,sp.Integer(2),"2",None):
        with pytest.raises((ValueError,TypeError)):
            fixtures.build(bad)
    for bad in (True,sp.true,1.0,"1"):
        with pytest.raises((ValueError,TypeError)):
            constraints.context(((bad,0,0),(-1,0,0)))


def test_exact_Fraction_and_rational_momentum_components():
    data=constraints.context(((Fraction(1,2),0,0),(sp.Rational(-1,2),0,0)))
    assert data.momenta==((sp.Rational(1,2),0,0),(sp.Rational(-1,2),0,0))


def test_full_fail_closed_control_inventory():
    assert verify.controls()["rejected_inputs"]==65


def test_exact_serialization_round_trip_without_float():
    data=verify.serialize({"canonical":canonical.data(),"norms":projector.data()})
    assert json.loads(json.dumps(data))==data
    assert json.loads(json.dumps(list(quadratic.CHANNELS)))==list(quadratic.CHANNELS)
    for value in (1.0,sp.Float(1),sp.oo,sp.nan):
        with pytest.raises((ValueError,TypeError)):
            verify.serialize(value)


def test_fourth_order_tracefree_correction_drops_at_background_only():
    ctx=constraints.context(((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)))
    W=[ctx.leg(0,1),ctx.leg(0,2),ctx.leg(0,3)]
    lift=constraints.old.york(W)
    assert j.trace(lift).is_zero()
    # The same tensor is not zero and can couple to a perturbed metric.
    assert any(not entry.is_zero() for row in lift for entry in row)
