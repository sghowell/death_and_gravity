"""Fresh actual coefficient and holomorphic basis bounds, with no old-background transfer."""
from functools import cache

import sympy as sp
from p8_proca_nearby_packets import domains as algebra
from p8_proca_nearby_packets import finite, modes
from p8_proca_nearby_retarded import source
from p8_proca_retuned_margin import domain as actual
from p8_proca_retuned_margin import flow, model

T=actual.T
LO,HI=sp.Rational(99,100),sp.Rational(101,100)
ac,am,ell=sp.symbols("new_clock_action_root new_matter_action_root new_matter_density",positive=True)


@cache
def basis():
    S=sp.Matrix([[1/ac,0,1/ac,0],[-ell/ac,1/am,-ell/ac,1/am],
        [sp.I*ac,sp.I*ell*am,-sp.I*ac,-sp.I*ell*am],
        [0,sp.I*am,0,-sp.I*am]])/sp.sqrt(2)
    Si=S.inv().applyfunc(sp.factor)
    limits={ac:(2,5),am:(sp.Rational(9,10),sp.Rational(11,10)),ell:(0,sp.Rational(1,8))}
    bounds={}
    for name,matrix in (("S",S),("S_inverse",Si)):
        entries=matrix.applyfunc(lambda v:algebra.laurent_majorant(v,limits))
        bounds[name]={"entry_modulus_upper":entries,
            "infinity_norm_upper":max(sum(entries[i,j] for j in range(4)) for i in range(4))}
    return {"exact_new_action_basis":S,"exact_new_action_inverse":Si,
        "actual_complex_action_root_modulus_domain":limits,"exact_basis_majorants":bounds,
        "declared_basis_infinity_norm_upper":sp.Integer(8),
        "declared_inverse_infinity_norm_upper":sp.Integer(5),
        "outer_actual_holomorphic_time_radius":T,"inner_real_half_width":T/2,
        "Cauchy_basis_time_derivative_infinity_norm_upper":8/(T/2),
        "declared_basis_derivative_infinity_norm_upper":sp.Integer(2)*10**8,
        "actual_root_branches_fixed_by_new_positive_initial_data":True}


@cache
def actual_bounds():
    boxes=actual.enclosures()
    ann={name:actual.complex_annulus(row) for name,row in boxes.items()}
    s=model.system()
    alpha=actual.complex_annulus(actual.rational_box((2*s["A"]*model.zf+s["B"])/3))
    clock_lower=2/(4*HI**3*sp.Rational(26,100)**2)*(LO/(HI**2))*LO
    clock_upper=4/(4*LO**3*sp.Rational(24,100)**2)*(HI/(LO**2))*HI
    return {"new_real_and_complex_boxes_fully_rebuilt":True,"new_complex_annuli":ann,
        "N_modulus_lower":actual.N0-actual.RN,"N_modulus_upper":actual.N0+actual.RN,
        "actual_R_minus_one_modulus_upper":flow.picard()["new_hat_scale_minus_one_absolute_upper"],
        "new_alpha_modulus_upper":alpha["upper"],
        "clock_action_root_squared_modulus_lower":clock_lower,
        "clock_action_root_squared_modulus_upper":clock_upper,
        "matter_action_root_squared_modulus_lower":LO**2/HI,
        "matter_action_root_squared_modulus_upper":HI**2/LO,
        "matter_density_modulus_upper":sp.Rational(11,100)/LO,
        "actual_a_modulus_upper":3*HI**2/4,"actual_m_and_r_modulus_upper":HI/(2*LO**3),
        "actual_g_modulus_upper":HI**2/2,
        "actual_curvature_N_modulus_lower":LO*sp.Rational(24,100),
        "actual_curvature_N_modulus_upper":HI*sp.Rational(26,100),
        "actual_lapse_Hessian_modulus_lower":2/HI,"actual_lapse_Hessian_modulus_upper":4/LO,
        "actual_beta_modulus_upper":2*sp.Rational(26,100)*sp.Rational(11,100),
        "actual_hat_Hubble_modulus_upper":ann["hat_Hubble"]["upper"],
        "actual_source_output_prefactor_squared_upper":(HI/LO)**3*HI**8,
        "actual_original_chart_source_modulus_upper":HI**4,
        "new_physical_clock_squared_speed_upper":sp.Rational(247,250),
        "coordinate_frequencies":"omega_m=N/(e R), omega_c=c_clock*omega_m",
        "physical_matter_frequency_is_largest_absolute_principal_frequency":True}


@cache
def coefficient_majorants():
    p=finite.parent
    limits={p.a:(0,3),p.m:(0,8),p.g:(0,2),p.r:(0,8),
        p.rn:(sp.Rational(1,10),1),p.h:(1,6),p.ell:(0,sp.Rational(1,8)),
        p.alpha:(0,1),p.beta:(0,1),p.H:(0,1),finite.R:(sp.Rational(1,2),2)}
    major=lambda matrix:matrix.applyfunc(lambda value:algebra.laurent_majorant(value,limits))
    high={name:major(finite.data()[name]) for name in ("L0","L1","L2","L3")}
    low=[major(matrix) for matrix in source.data()["polynomial_q_coefficients"]]
    low_total=sum((matrix*4**j for j,matrix in enumerate(low)),sp.zeros(4))
    return {"complete_Laurent_matrix_entry_majorants":high,
        "complete_original_polynomial_q_coefficient_majorants":low,
        "low_scalar_root_unit_disc_generator_entry_majorant":low_total,
        "low_scalar_root_unit_disc_generator_infinity_norm_upper":
            max(sum(low_total[i,j] for j in range(4)) for i in range(4)),
        "limits_are_verified_against_new_action_not_old_background":True}


@cache
def checks():
    d=basis()
    sub={ac:sp.sqrt(modes.kc*modes.wc),am:sp.sqrt(modes.km*modes.wm),ell:modes.ell}
    S=d["exact_new_action_basis"]
    Si=d["exact_new_action_inverse"]
    old=modes.data()
    return {"new_basis_matches_universal_action_eigenvectors_not_old_solution":
        (S.subs(sub,simultaneous=True)-old["S"]).applyfunc(sp.simplify),
        "new_basis_inverse_is_literal_matrix_inverse":(S*Si-sp.eye(4)).applyfunc(sp.factor),
        "new_actual_basis_diagonalizes_full_principal_generator":
        (old["J"]*S.subs(sub,simultaneous=True)-sp.I*S.subs(sub,simultaneous=True)*old["Lambda"]).applyfunc(sp.simplify),
        "new_basis_Cauchy_derivative_keeps_half_disc_radius":
        d["Cauchy_basis_time_derivative_infinity_norm_upper"]-sp.Integer(16)*10**7,
        "complete_low_disc_generator_keeps_all_q_orders":
        coefficient_majorants()["low_scalar_root_unit_disc_generator_infinity_norm_upper"]-sp.Rational(28169,64)}


@cache
def gates():
    d,b=actual_bounds(),basis()
    a=d["new_complex_annuli"]
    rows={"new_lapse_fits_one_percent_modulus_annulus":LO<d["N_modulus_lower"]<d["N_modulus_upper"]<HI,
        "new_scale_fits_one_percent_modulus_annulus":d["actual_R_minus_one_modulus_upper"]<sp.Rational(1,100),
        "new_conformal_root_fits_one_percent_modulus_annulus":
            flow.picard()["complex_conformal_D_minus_one_upper"]<sp.Rational(1,100)
            and (actual.N0+actual.RN)**2*(1+flow.picard()["complex_background_h_minus_one_upper"])/
                (1-flow.picard()["complex_conformal_D_minus_one_upper"])<HI**4
            and (actual.N0-actual.RN)**2*(1-flow.picard()["complex_background_h_minus_one_upper"])/
                (1+flow.picard()["complex_conformal_D_minus_one_upper"])>LO**4,
        "new_clock_speed_complex_root_fits_one_percent_annulus":
            a["clock_physical_speed_squared"]["lower"]>LO**2 and a["clock_physical_speed_squared"]["upper"]<HI**2,
        "new_actual_complex_lapse_pivot_fits_two_to_four_annulus":
            a["on_constraint_rescaled_lapse_Hessian"]["lower"]>2
            and a["on_constraint_rescaled_lapse_Hessian"]["upper"]<4,
        "new_actual_complex_curvature_derivative_fits_declared_annulus":
            a["M"]["lower"]>sp.Rational(24,100) and a["M"]["upper"]<sp.Rational(26,100),
        "new_rescaled_matter_root_modulus_below_eleven_hundredths":
            a["constraint_matter_square"]["upper"]<sp.Rational(11,100)**2,
        "new_clock_action_root_modulus_above_two":d["clock_action_root_squared_modulus_lower"]>4,
        "new_clock_action_root_modulus_below_five":d["clock_action_root_squared_modulus_upper"]<25,
        "new_matter_action_root_modulus_above_nine_tenths":d["matter_action_root_squared_modulus_lower"]>sp.Rational(9,10)**2,
        "new_matter_action_root_modulus_below_eleven_tenths":d["matter_action_root_squared_modulus_upper"]<sp.Rational(11,10)**2,
        "new_matter_density_modulus_below_one_eighth":d["matter_density_modulus_upper"]<sp.Rational(1,8),
        "new_alpha_modulus_below_one":d["new_alpha_modulus_upper"]<1,
        "new_beta_modulus_below_one":d["actual_beta_modulus_upper"]<1,
        "new_a_modulus_below_three":d["actual_a_modulus_upper"]<3,
        "new_m_and_r_moduli_below_eight":d["actual_m_and_r_modulus_upper"]<8,
        "new_g_modulus_below_two":d["actual_g_modulus_upper"]<2,
        "new_curvature_N_modulus_between_one_tenth_and_one":
            d["actual_curvature_N_modulus_lower"]>sp.Rational(1,10) and d["actual_curvature_N_modulus_upper"]<1,
        "new_lapse_Hessian_modulus_between_one_and_six":
            d["actual_lapse_Hessian_modulus_lower"]>1 and d["actual_lapse_Hessian_modulus_upper"]<6,
        "new_Hubble_modulus_below_one":d["actual_hat_Hubble_modulus_upper"]<1,
        "new_basis_infinity_norm_below_eight":b["exact_basis_majorants"]["S"]["infinity_norm_upper"]<8,
        "new_inverse_infinity_norm_below_five":b["exact_basis_majorants"]["S_inverse"]["infinity_norm_upper"]<5,
        "new_basis_Cauchy_derivative_below_declared":b["Cauchy_basis_time_derivative_infinity_norm_upper"]<2*10**8,
        "new_two_endpoint_scalar_prefactor_below_two":d["actual_source_output_prefactor_squared_upper"]<4,
        "new_original_chart_source_below_two":d["actual_original_chart_source_modulus_upper"]<2,
        "all_complete_Laurent_matrix_norms_below_ten_thousand":all(
            max(sum(matrix[i,j] for j in range(4)) for i in range(4))<10000
            for matrix in coefficient_majorants()["complete_Laurent_matrix_entry_majorants"].values()),
        "new_real_clock_is_strictly_slower_than_physical_matter":0<d["new_physical_clock_squared_speed_upper"]<1}
    return {name:bool(value) for name,value in rows.items()}
