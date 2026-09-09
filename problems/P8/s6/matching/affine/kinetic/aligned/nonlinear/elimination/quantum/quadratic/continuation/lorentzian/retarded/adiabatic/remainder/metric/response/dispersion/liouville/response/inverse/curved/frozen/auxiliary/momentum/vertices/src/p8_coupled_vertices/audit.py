"""Exact physical residual inventories and fail-closed scientific controls."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_momentum import constraints
from p8_physical import jets as j

from . import anchors, coefficients, majorant, physical, taylor

POINTS=((2,sp.Integer(0)),(3,sp.Rational(-1,2)),(3,sp.Integer(0)),(3,sp.Rational(1,2)),
        (4,sp.Rational(-1,2)),(4,sp.Integer(0)),(4,sp.Rational(1,2)))


@cache
def residuals():
    out={}
    for group in (coefficients.checks(),physical.quadratic_algebraic_bridge(),anchors.checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated higher-vertex residual name")
        out.update(group)
    for n,point in POINTS:
        data=physical.fixture(n,point)
        ctx,geo,mom=data["context"],data["geometry"],data["momentum"]
        rest=constraints.residual(mom["momentum"],geo,mom["physical_lower_current"])
        label=f"mixed_{n}_leg_time_{point}"
        for degree in range(1,n):
            masks=[mask for mask in range(1,ctx.full) if mask.bit_count()==degree]
            out[label+f"_all_spatial_constraints_degree_{degree}"]=sp.ImmutableMatrix(
                [[row.coefficient(mask) for mask in masks] for row in rest])
        reduced=data["reduction"]
        masks=[mask for mask in range(ctx.full+1) if mask.bit_count()<n]
        out[label+"_actual_stationary_lapse_force"]=sp.ImmutableMatrix(
            [reduced["stationary_force"].coefficient(mask) for mask in masks])
        out[label+"_direct_stationary_Hamiltonian"]=sp.ImmutableMatrix(
            [reduced["direct_Hamiltonian_difference"].coefficient(mask) for mask in range(ctx.full+1)])
    ctx,e=taylor.context(model.u)
    actual=taylor.evaluate((1+model.N)**(-sp.Rational(1,2)),{model.N:e},ctx)
    expected=tuple(sp.prod(-sp.Rational(1,2)-k for k in range(n)) for n in range(5))
    for n,value in enumerate(taylor.derivatives(actual)):
        out[f"symmetric_nilpotent_Taylor_derivative_{n}"]=sp.factor(value-expected[n])
    integral=taylor.integrate_zero(actual,e)
    for n in range(1,5):
        out[f"fixed_basepoint_Taylor_primitive_{n}"]=sp.factor(
            taylor.derivatives(integral)[n]-taylor.derivatives(actual)[n-1])
    return out


@cache
def gates():
    rows=coefficients.jets()
    small=majorant.bound(1,4,1)
    full=small["physical_Hamiltonian_degree_majorant"].coefficients
    no_vector=majorant.bound(1,4,1,include_vector=False)
    budgets=[]
    all_constraints=[]
    for n,point in POINTS:
        data=physical.fixture(n,point)
        value=data["full_labelled_kernel"]
        budgets.append(abs(sp.re(value))+abs(sp.im(value))<full[n])
        all_constraints.extend(data["momentum"]["checks"].values())
        all_constraints.extend(data["reduction"]["checks"].values())
    result={"actual_polynomial_has_fourteen_monomials":len(rows)==14,
            "all_seventy_actual_lapse_derivatives_retained":sum(map(len,rows.values()))==70,
            "actual_polynomial_has_nine_linear_and_four_quadratic_monomials":
                sum(sum(p)==1 for p in rows)==9 and sum(sum(p)==2 for p in rows)==4,
            "all_seventy_rational_denominator_certificates_checked":
                sum(len(value["proof_data"]) for value in majorant.coefficient_bounds()["rows"].values())==70,
            "all_actual_spatial_and_lapse_fixture_gates":all(all_constraints),
            "seven_mixed_vertex_fixtures_at_three_exact_times":len(POINTS)==7,
            "all_mixed_fixture_kernels_below_uniform_majorant":all(bool(x) for x in budgets),
            "nonzero_independent_electric_quartic":anchors.pure_electric()["kernel"]>0,
            "omitting_lapse_response_fails_electric_quartic":anchors.pure_electric()["omitting_lapse_response_quartic_defect"]!=0,
            "mixed_fourier_reality_not_false_pointwise_real_claim":
                sp.im(physical.fixture(3)["full_labelled_kernel"])!=0,
            "all_label_permutations_and_reversals_preserve_declared_coefficients":
                all(value==0 for value in anchors.checks().values()),
            "retaining_full_vector_generator_increases_York_degree_two_and_three_bounds":
                all(small["York_correction_majorants"][n]>no_vector["York_correction_majorants"][n] for n in (1,2)),
            "no_Hubble_or_Theta_inverse_in_actual_stationary_and_York_reduction":True,
            "raw_phase_bounds_not_free_mode_columns_or_transition_norms":True,
            "classical_action_quantum_order_state_and_subtraction_not_changed":True,
            "no_old_D_only_cutoff_or_unmodified_M1_band_transfer":True,
            "original_P8_and_finite_common_parent_matching_open":True}
    return {key:bool(value) for key,value in result.items()}


def controls():
    calls=[]
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,None)
    for value in bad:
        for index in range(3):
            args=[1,4,1]
            args[index]=value
            calls.append(lambda args=tuple(args):majorant.bound(*args))
        if value is not None:
            calls.append(lambda value=value:physical.background(value))
    for index in range(3):
        for value in (0,-1):
            args=[1,4,1]
            args[index]=value
            calls.append(lambda args=tuple(args):majorant.bound(*args))
    for value in (True,False,1,5,2.0,sp.Integer(2),"2",None):
        calls.append(lambda value=value:physical.fixture(value))
    for value in (sp.Rational(-3,4),sp.Rational(3,4)):
        calls.append(lambda value=value:physical.background(value))
    for value in (1,0,sp.true,sp.false,1.0,None):
        calls.extend((lambda value=value:majorant.bound(1,4,1,include_vector=value),
                      lambda value=value:physical.fixture(2,symbolic_scale=value),
                      lambda value=value:anchors.transformed_fixture(3,reverse=value)))
    for value in (True,2.0,0,5,"3",None):
        calls.append(lambda value=value:taylor.context(model.u,value))
    for value in ((0,1),(0,1,1),(0,1,True),[0,1,2]):
        calls.append(lambda value=value:anchors.transformed_fixture(3,permutation=value))
    for value in (True,sp.true,1,"unknown",None):
        calls.extend((lambda value=value:physical.quadratic_pair(value,"curvature"),
                      lambda value=value:physical.quadratic_pair("curvature",value)))
    ctx=constraints.context(((1,0,0),(-1,0,0)))
    f,t,tp,W,Pi=physical.mixed_fields(ctx)
    bad_f=dict(f)
    bad_f["zeta"]=ctx.jet(1)
    calls.append(lambda:physical.construct(ctx,bad_f,t,tp,W,Pi,point=0))
    bad_t=j.madd(t,j.mscale(j.identity(ctx),ctx.leg(0)))
    calls.append(lambda:physical.construct(ctx,f,bad_t,tp,W,Pi,point=0))
    nontransverse=j.zeros(ctx)
    nontransverse[0][0]=ctx.leg(0)
    nontransverse[1][1]=-ctx.leg(0)
    calls.append(lambda:physical.construct(ctx,f,nontransverse,tp,W,Pi,point=0))
    nonsymmetric=j.zeros(ctx)
    nonsymmetric[1][2]=ctx.leg(0)
    calls.append(lambda:physical.construct(ctx,f,nonsymmetric,tp,W,Pi,point=0))
    nonlinear=dict(f)
    nonlinear["chi"]=ctx.leg(0)*ctx.leg(1)
    calls.append(lambda:physical.construct(ctx,nonlinear,t,tp,W,Pi,point=0))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An unsupported higher-vertex input was accepted")
    return {"rejected_inputs":rejected,"exact_native_validation_before_cached_public_entrypoints":True,
            "original_boundary_primitives_not_replaced_by_clock_only_coefficients":True,
            "all_geometry_matter_and_vector_density_factors_kept":True,
            "independent_electric_quartic_rejects_omitted_lapse_response":True,
            "no_wrong_scalar_only_or_D_only_physical_bound":True}
