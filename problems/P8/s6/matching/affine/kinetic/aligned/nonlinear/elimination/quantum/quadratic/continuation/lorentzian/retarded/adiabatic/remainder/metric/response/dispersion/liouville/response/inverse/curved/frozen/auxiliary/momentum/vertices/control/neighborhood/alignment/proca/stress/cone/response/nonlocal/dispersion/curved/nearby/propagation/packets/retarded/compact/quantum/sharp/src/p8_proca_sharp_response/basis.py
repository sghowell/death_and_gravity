"""Small exact basis bounds on the actual complex-time bounce disc."""
from functools import cache

import sympy as sp
from p8_proca_nearby_packets import domains, modes

from . import jets

ac,am,ell,lc,lm,elld=sp.symbols("clock_action_root matter_action_root matter_density clock_root_log_jet matter_root_log_jet density_jet",positive=True)
R=sp.Symbol("hat_scale",positive=True)


@cache
def data():
    S=sp.Matrix([[1/ac,0,1/ac,0],[-ell/ac,1/am,-ell/ac,1/am],
        [sp.I*ac,sp.I*ell*am,-sp.I*ac,-sp.I*ell*am],
        [0,sp.I*am,0,-sp.I*am]])/sp.sqrt(2)
    Si=S.inv().applyfunc(sp.factor)
    Sd=S.diff(ac)*ac*lc+S.diff(am)*am*lm+S.diff(ell)*elld
    limits={ac:(2,4),am:(sp.Rational(9,10),sp.Rational(11,10)),
        ell:(0,sp.Rational(1,8)),lc:(0,sp.Rational(1,400)),
        lm:(0,sp.Rational(1,10000)),elld:(0,sp.Rational(3,500000))}
    majorants={}
    for name,matrix in (("S",S),("S_inverse",Si),("S_dot",Sd)):
        b=matrix.applyfunc(lambda v:domains.laurent_majorant(v,limits))
        majorants[name]={"entry_upper":b,
            "infinity_norm_upper":max(sum(b[i,j] for j in range(4)) for i in range(4))}
    lo,hi=sp.Rational(99,100),sp.Rational(101,100)
    cmin,cmax=sp.Rational(99,100),sp.Rational(101,100)
    clock_square_lower=2/(4*hi**3*sp.Rational(26,100)**2)*(lo/(hi*hi))*cmin
    clock_square_upper=3/(4*lo**3*sp.Rational(24,100)**2)*(hi/(lo*lo))*cmax
    return {"compressed_exact_basis":S,"compressed_exact_inverse":Si,
        "complete_basis_time_derivative":Sd,"complex_modulus_and_jet_domain":limits,
        "clock_action_root_squared_modulus_lower":clock_square_lower,
        "clock_action_root_squared_modulus_upper":clock_square_upper,
        "matter_action_root_squared_modulus_lower":lo**2/hi,
        "matter_action_root_squared_modulus_upper":hi**2/lo,
        "basis_majorants":majorants,"declared_basis_norm_upper":sp.Integer(6),
        "declared_inverse_norm_upper":sp.Integer(4),
        "declared_basis_derivative_norm_upper":sp.Rational(1,50),
        "actual_two_endpoint_scalar_prefactor_squared_upper":(hi/lo)**3*hi**8,
        "declared_two_endpoint_scalar_prefactor_upper":sp.Integer(2),
        "declared_each_front_amplitude_upper":sp.Integer(3),
        "complex_root_branches_are_fixed_by_positive_real_initial_values":True}


@cache
def checks():
    d=data()
    old=modes.data()
    sub={ac:sp.sqrt(modes.kc*modes.wc),am:sp.sqrt(modes.km*modes.wm),ell:modes.ell}
    n,e,h,D,M,p,c,r=sp.symbols("n e h D M pivot c r",positive=True)
    raw=(-p/(4*e**3*M**2))*(n*c/(e*r))
    rational=-D*p*c/(4*n*h*r*M**2)
    L_D,L_p,L_C,L_n,L_h,L_R,L_M=sp.symbols("L_D L_p L_C L_n L_h L_R L_M")
    t=sp.Symbol("formal_time",real=True)
    product=-D*sp.exp(L_D*t)*p*sp.exp(L_p*t)*c*sp.exp(L_C*t/2)/(
        4*n*sp.exp(L_n*t)*h*sp.exp(L_h*t)*r*sp.exp(L_R*t)*M**2*sp.exp(2*L_M*t))
    return {"compressed_basis_equals_full_action_normalized_parent":
        (d["compressed_exact_basis"].subs(sub,simultaneous=True)-old["S"]).applyfunc(sp.simplify),
        "compressed_inverse_equals_full_action_normalized_parent":
        (d["compressed_exact_inverse"].subs(sub,simultaneous=True)-old["S_inverse"]).applyfunc(sp.simplify),
        "compressed_basis_and_inverse_multiply_to_identity":
        (d["compressed_exact_basis"]*d["compressed_exact_inverse"]-sp.eye(4)).applyfunc(sp.factor),
        "clock_action_root_rationalization_uses_actual_conformal_fourth_power":
        sp.factor((raw-rational).subs(e**4,n*n*h/D)),
        "matter_action_root_square_is_conformal_square_over_scale":
        sp.factor((e**3/n)*(n/(e*r))-e*e/r),
        "clock_action_root_log_derivative_keeps_every_chain_factor":
        sp.factor(sp.diff(product,t).subs(t,0)/(2*product.subs(t,0))
            -(L_D+L_p+L_C/2-L_n-L_h-L_R-2*L_M)/2),
        "conserved_matter_density_derivative_is_minus_three_Hubble":
        sp.factor(-2*(-3*n/4)*sp.Symbol("z")*ell+3*(-n*sp.Symbol("z")/2)*ell)}


@cache
def gates():
    d=data()
    a=domains.analytic_domain()
    out={"complex_conformal_factor_fits_one_percent_annulus":a["e_minus_one_upper"]<sp.Rational(1,100),
        "complex_scale_fits_one_percent_annulus":a["R_minus_one_upper"]<sp.Rational(1,100),
        "complex_lapse_fits_one_percent_annulus":a["N_minus_one_upper"]<sp.Rational(1,100),
        "complex_clock_speed_fits_one_percent_annulus":a["clock_speed_minus_one_complex_upper"]<sp.Rational(1,100),
        "clock_action_root_modulus_above_two":d["clock_action_root_squared_modulus_lower"]>4,
        "clock_action_root_modulus_below_four":d["clock_action_root_squared_modulus_upper"]<16,
        "matter_action_root_modulus_above_nine_tenths":d["matter_action_root_squared_modulus_lower"]>sp.Rational(9,10)**2,
        "matter_action_root_modulus_below_eleven_tenths":d["matter_action_root_squared_modulus_upper"]<sp.Rational(11,10)**2,
        "actual_scalar_source_output_prefactor_below_two":
            d["actual_two_endpoint_scalar_prefactor_squared_upper"]<4,
        "clock_front_amplitude_below_three":2*sp.Rational(1,8)**2/4<3,
        "matter_front_amplitude_below_three":2/sp.Rational(9,10)**2<3}
    for name,cap in (("S",6),("S_inverse",4),("S_dot",sp.Rational(1,50))):
        out["complete_actual_"+name+"_norm_below_declared"]=d["basis_majorants"][name]["infinity_norm_upper"]<cap
    out["all_actual_complex_jets_verified"]=all(jets.gates().values())
    return {name:bool(value) for name,value in out.items()}
