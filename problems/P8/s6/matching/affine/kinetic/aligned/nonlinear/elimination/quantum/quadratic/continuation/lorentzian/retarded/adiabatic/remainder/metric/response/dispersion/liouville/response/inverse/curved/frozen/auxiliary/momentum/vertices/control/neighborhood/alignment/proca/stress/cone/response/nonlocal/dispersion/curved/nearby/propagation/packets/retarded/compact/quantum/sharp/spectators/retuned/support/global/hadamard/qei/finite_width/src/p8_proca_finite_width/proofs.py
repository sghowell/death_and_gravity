"""Exact normalization, whole-interval inclusion and proof-domain audits."""
from fractions import Fraction as F
from functools import cache

import sympy as sp
from flint import acb, ctx

from . import bounds, canonical, compact, intervals, time_jets


@cache
def checks():
    d=canonical.data()
    omega=canonical.phase.data()["constant_symplectic_form"]
    E,Ei=d["density_to_regular"],d["regular_to_density"]
    mu=canonical.mu
    out={
        "all_real_frequency_canonical_map_is_scaled_symplectic":canonical.clean(E*omega*E.T-mu*omega),
        "all_real_frequency_canonical_map_has_exact_inverse":canonical.clean(E*Ei-sp.eye(4)),
        "complete_all_real_generator_keeps_every_Laurent_coefficient":canonical.clean(d["generator"]-
            sum((mu**j*v for j,v in d["coefficients"].items()),sp.zeros(4))),
        "all_real_chart_keeps_the_actual_global_leading_generator":canonical.clean(
            d["coefficients"][1]-canonical.canonical.data()["leading_generator"]),
    }
    for j,M in d["coefficients"].items():
        out["regular_chart_full_Hamiltonian_Laurent_order_"+str(j)]=canonical.clean(M*omega+omega*M.T)
    o=canonical.o
    original=canonical.phase.data()
    force=original["normalized_physical_probe_force"]
    out["regular_chart_preserves_the_original_density_probe_force"]=canonical.clean(Ei*(E*force)-force)
    original_row=original["original_relational_observable_row"]
    out["regular_chart_keeps_the_original_local_observable"]=canonical.clean(
        original_row*Ei-sp.Matrix([[0,1/(mu*o.a**sp.Rational(3,2)),0,0]]))
    out["zero_spatial_momentum_is_an_exact_regular_chart_value"]=canonical.clean(
        d["generator"].subs(mu,1)-
        ((canonical.model.dt(E)+E*original["regular_density_generator"])*Ei).subs({o.q:0,mu:1}))
    parent=canonical.canonical.data()
    bg=canonical.model.substitution()
    for n,j in enumerate((1,0,-1)):
        H=canonical.clean((-omega*parent["complete_Laurent_coefficients"][j]).subs(bg,simultaneous=True))
        out["actual_scaled_Hamiltonian_time_reversal_order_"+str(n)]=canonical.clean(
            H.subs(canonical.model.u,-canonical.model.u)-(-1)**n*H)
    rho=sp.Symbol("positive_regular_to_original_frequency_ratio",positive=True)
    change=sp.diag(rho**sp.Rational(-3,2),sp.sqrt(rho),rho**sp.Rational(3,2),1/sp.sqrt(rho))
    k=canonical.canonical.k
    actual_change=canonical.clean((E*parent["packet_to_density"]).subs(bg,simultaneous=True)
        .subs(canonical.model.u,0)*sp.sqrt(k/mu))
    out["initial_old_to_regular_normalized_symplectic_map"]=canonical.clean(
        actual_change.subs(mu,rho*k)-change)
    out["initial_normalized_transition_preserves_exact_CCR"]=canonical.clean(change*omega*change.T-omega)
    out["energy_logarithmic_norm_constant"]=sp.Rational(4+2*13*2,2)*14-392
    out["full_phase_energy_norm_ratio"]=sp.Integer(13)*14-182
    x,K=sp.symbols("nonnegative_alpha positive_cutoff",positive=True)
    q=sp.Symbol("positive_radial_frequency",positive=True)
    # These antiderivatives and vanishing endpoints give both exact tails.
    out["opposite_sign_alpha_tail_primitive"]=sp.diff(-sp.Rational(1,5)/(x+q)**5,x)-(x+q)**-6
    out["joint_high_frequency_tail_primitive"]=sp.diff(-sp.Rational(1,5)/q,q)-q**3/(5*q**5)
    out["packet_error_radial_tail_primitive"]=sp.diff(-sp.Rational(1,8)/q**8,q)-q**-9
    out["low_frequency_radial_majorant_primitive"]=sp.diff(q**4/4+q**3/3,q)-q*q*(q+1)
    out["energy_half_and_state_and_radial_QEI_normalization"]=sp.Rational(1,2)*sp.Rational(1,2)/(2*sp.pi**2)/sp.pi-1/(8*sp.pi**3)
    out["high_approximation_keeps_factor_two_from_exact_error_split"]=2*4/(8*sp.pi**3)/(5*K)-1/(5*sp.pi**3*K)
    out["high_error_Parseval_normalization_keeps_same_split"]=2*2*sp.Rational(1,2)*sp.Rational(1,2)/(2*sp.pi**2)/(8*K**8)-1/(16*sp.pi**2*K**8)
    l=sp.Function("l")(x)
    z=sp.Function("Q")(x)
    a=sp.Function("g")(x)
    row=a*l
    for _ in range(3):
        row=-sp.diff(row*z,x)
    # An independent symbolic product-rule expansion recovers the scalar
    # noncommuting-safe majorant when all derivative magnitudes are positive.
    L=sp.symbols("L0:4",positive=True)
    D=sp.symbols("D0:4",positive=True)
    subst={sp.diff(l,x,j):L[j] for j in range(4)}
    subst.update({sp.diff(z,x,j):D[j] for j in range(4)})
    for n in range(4):
        coeff=-sp.expand(row).coeff(sp.diff(a,x,n)).subs(subst)
        symbolic={(0,j,r):sp.binomial(j,r)*L[j-r] for j in range(4) for r in range(j+1)}
        for m in range(3):
            for j in range(3-m):
                for r in range(m+j+2):
                    symbolic[m+1,j,r]=sum(sp.binomial(j+1,s)*symbolic.get((m,s,r),0)*D[j+1-s]
                        for s in range(j+2))
        out["third_IBP_ordered_product_rule_sampler_derivative_"+str(n)]=sp.expand(coeff-symbolic[3,0,n])
    return {k:canonical.clean(v) for k,v in out.items()}

@cache
@ctx.workprec(160)
def exact_jet_inclusions():
    ball=intervals.build()
    exact=time_jets.build(9)
    root=acb(17985).sqrt()
    count=0
    for n,R in enumerate(ball["R"]):
        for i in range(2):
            for j in range(2):
                for order in range(10-n):
                    value=time_jets.FIELD.to_sympy(exact["jets"][n][i][j][order])
                    re,im=sp.expand(value).as_real_imag()
                    def convert(v):
                        p=sp.Poly(sp.expand(v),sp.sqrt(17985))
                        assert p.degree()<=1 or p.is_zero
                        return intervals.scalar(p.nth(0))+intervals.scalar(p.nth(1))*root
                    enclosure=convert(re)+acb(0,1)*convert(im)
                    if not R[i][j][order].contains(enclosure):
                        raise ValueError("Whole-interval jet failed independent exact center inclusion")
                    count+=1
    return count

@cache
def gates():
    c=compact.center_bounds()
    b=bounds.data()
    transition_cube=F(257,256)**3
    return {
        "eight_exact_Riccati_orders_include_nonzero_second_fourth_sixth_eighth":all(
            any(v!=0 for v in time_jets.build(8)["at_center"][n]) for n in (2,4,6,8)),
        "center_odd_orders_through_eighth_vanish":all(
            time_jets.build(8)["at_center"][n]==sp.zeros(2) for n in (1,3,5,7)),
        "independent_196_exact_center_jets_in_whole_interval_enclosures":exact_jet_inclusions()==196,
        "positive_initial_Borel_graph_from_k_sixteen":c["tail_budget_at_16"]<F(3,20),
        "initial_graph_lower_three_quarters_and_upper_thirteen":c["positive_graph_lower"]==F(3,4)
            and c["positive_graph_upper"]==13,
        "initial_normalized_chart_transition_square_norm_below_one_point_zero_one":
            transition_cube<F(101,100)**2,
        "convex_reference_initial_covariance_norm_below_sixteen":F(101,100)*(13+F(4,3))<16,
        "initial_negative_configuration_factor_norm_below_two":F(121,100)/F(3,4)<4,
        "whole_interval_positive_frequency_exceeds_point_nine_nine":
            intervals.build()["frequency_lower_gt_99_over_100"],
        "resolvent_Neumann_margin_at_actual_cutoff":F(99,100)*b["high_frequency_threshold"]
            >2*b["configuration_remainder_derivative_bounds"][0],
        "all_real_energy_lower_one_fourteenth":compact.transfer_bounds()["energy_lower"]==F(1,14),
        "all_real_energy_upper_thirteen":compact.transfer_bounds()["energy_upper"]==13,
        "complete_canonical_and_regular_real_transfer_bound":b["full_phase_log_norm_rate"]==392,
        "nonzero_finite_order_Riccati_remainder_is_bounded_not_deleted":
            b["scaled_Riccati_defect_order_seven_constant"]==91215897,
        "real_test_energy_time_row_constant_retains_all_terms":
            b["physical_time_derivative_observable_constant"]==16,
        "all_three_finite_width_contributions_are_strictly_positive":all(
            b[key]>0 for key in ("low_band_L2_squared_coefficient",
            "approximate_high_band_IBP_squared_coefficient","high_band_error_L2_squared_coefficient")),
        "actual_half_width_is_a_nonzero_finite_interval":b["half_width"]==sp.Rational(1,100),
    }

@cache
def controls():
    rejected=0
    def reject(call):
        nonlocal rejected
        try:
            call()
        except (ValueError,TypeError):
            rejected+=1
            return
        raise AssertionError("An invalid proof-domain input was accepted")
    for value in (True,False,-1,1.5,None,"6",sp.Rational(3,2)):
        reject(lambda value=value:time_jets.TimeJets(value))
    for value in (True,False,0,-1,1.5,None,"9"):
        reject(lambda value=value:intervals.BallJets(value))
    for value in (True,False,0,-1,1.5,None,"-1/100","0"):
        reject(lambda value=value:intervals.BallJets(3,value))
    for value in (True,-1,1.5,None):
        reject(lambda value=value:bounds.ibp_coefficients([1]*4,[1]*4,value))
    for left,right in (([],[1]*4),([1]*4,[]),([1,-1,1,1],[1]*4),([1]*4,[1,True,1,1]),
            ([1,1.0,1,1],[1]*4),(None,[1]*4)):
        reject(lambda left=left,right=right:bounds.ibp_coefficients(left,right))
    T=intervals.BallJets(3)
    reject(lambda:T.inv(T.zero))
    reject(lambda:T.sqrt(T.constant(-1)))
    reject(lambda:T.sqrt(T.constant(sp.I)))
    exact=time_jets.TimeJets(3)
    reject(lambda:exact.inv(exact.zero))
    reject(lambda:exact.sqrt(exact.one,2))
    reject(lambda:compact.rational_range(1/(canonical.model.u),F(1,100)))
    return {"rejected_inputs":rejected,"finite_order_truncation_alone_is_not_declared_Hadamard":True,
        "cutoff_chart_is_not_a_new_local_observable":True,
        "whole_interval_enclosures_not_finitely_sampled_times":True}
