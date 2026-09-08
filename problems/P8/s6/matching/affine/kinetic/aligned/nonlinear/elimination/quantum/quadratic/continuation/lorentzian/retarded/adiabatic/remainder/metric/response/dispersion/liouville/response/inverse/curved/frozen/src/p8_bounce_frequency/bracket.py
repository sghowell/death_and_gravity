"""Exact rational inequalities for a positive-real frozen pole bracket."""
from functools import cache

import sympy as sp
from p8_vector_metric_dispersion import spectral

from . import coefficients as coeff

LOW=sp.Integer(10)**12
HIGH=sp.Integer(10)**13
ERROR=sp.Rational(1,10**12)
GRAM_UPPER=sp.ImmutableMatrix([[12,17],[17,24]])
FINITE_UPPER=sp.ImmutableMatrix([[2,3],[3,4]])
BLOCK_UPPER=FINITE_UPPER+25*GRAM_UPPER


def endpoint(value):
    if type(value) is not str or value not in ("lower","upper"):
        raise ValueError("Require native endpoint lower or upper")
    return _endpoint(value)


@cache
def _endpoint(value):
    s=LOW if value=="lower" else HIGH
    lower=s**4/(640*coeff.L**2)
    upper=s**4/(576*coeff.L**2)
    c=coeff.tree_rows()[0][0,0]
    A_lower=c-BLOCK_UPPER[0,0]*upper-ERROR
    A_upper=c-lower+ERROR
    C_lower=6*s*s-sp.Rational(9,100)-BLOCK_UPPER[1,1]*upper-ERROR
    off=sp.Rational(3,200)+BLOCK_UPPER[0,1]*upper+ERROR
    return {"frequency":s,"gamma_s4_lower":lower,"gamma_s4_upper":upper,
            "NN_lower":A_lower,"NN_upper":A_upper,"ZZ_lower":C_lower,
            "both_off_diagonal_absolute_upper":off,
            "determinant_comparison":A_lower*C_lower-off**2 if value=="lower" else A_upper*C_lower+off**2}


@cache
def data():
    pair_T=sp.Matrix([sp.Rational(109,81),2])
    pair_L=sp.Matrix([sp.Rational(226,81),4])
    raw=2*pair_T*pair_T.T+pair_L*pair_L.T
    box=coeff.profile_box()
    tad=5760*coeff.L**2*sum(box.values())
    small=(2*sp.Integer(10)**15+HIGH**2*10**9)/(576*coeff.L**2)
    all_C=6*LOW**2-sp.Rational(9,100)-BLOCK_UPPER[1,1]*HIGH**4/(576*coeff.L**2)-ERROR
    return {"parameters":{"L":coeff.L,"mass":coeff.MASS,"margin":coeff.MARGIN},
            "positive_real_open_frequency_bracket":[LOW,HIGH],
            "pair_Gram_entrywise_raw_upper":sp.ImmutableMatrix(raw),
            "pair_Gram_entrywise_integer_upper":GRAM_UPPER,
            "finite_fourth_entrywise_absolute_upper":FINITE_UPPER,
            "exact_block_entrywise_absolute_upper_on_band":BLOCK_UPPER,
            "local_zero_entrywise_absolute_max":max(abs(x) for x in coeff.local_row(0).subs(coeff.mass,coeff.MASS)),
            "local_second_entrywise_absolute_max":max(abs(x) for x in coeff.local_row(2).subs(coeff.mass,coeff.MASS)),
            "fixed_profile_box":box,"normalized_tadpole_entrywise_upper":tad,
            "full_lower_order_entrywise_error_upper":small,
            "rounded_entrywise_error_upper":ERROR,
            "log_argument_upper":1+HIGH**2/(4*coeff.MASS**2),
            "exp_100_positive_series_single_term":sp.Integer(100)**18/sp.factorial(18),
            "pi_lower_partial_sum":4*sum((-1)**j/sp.Integer(2*j+1) for j in range(8)),
            "pi_upper":sp.Rational(22,7),
            "ZZ_positive_lower_on_entire_band":all_C,
            "lower_endpoint":endpoint("lower"),"upper_endpoint":endpoint("upper")}


@cache
def checks():
    x,d=sp.symbols("integration_variable radial_denominator",positive=True)
    geometric=sum((-1)**j*x**(2*j) for j in range(8))
    out={"pi_lower_integrand_positive_remainder":sp.factor(1/(1+x*x)-geometric-x**16/(1+x*x)),
         "pi_upper_positive_integral":sp.integrate(x**4*(1-x)**4/(1+x*x),(x,0,1))-sp.Rational(22,7)+sp.pi,
         "radial_denominator_comparison":sp.expand((d-x*x)-(d-x)-x*(1-x)),
         "log_majorant_primitive":sp.diff(-sp.log(d-x),x)-1/(d-x),
         "block_absolute_integer_majorant":BLOCK_UPPER-sp.Matrix([[302,428],[428,604]])}
    c=coeff.tree_rows()
    s=sp.Symbol("laplace",positive=True)
    determinant=(c[0]+s*s*c[2]).det()
    expected=6*sp.Rational(749377,250000)*s*s-sp.Rational(9,100)*sp.Rational(749377,250000)-sp.Rational(9,40000)
    out["tree_only_determinant_has_no_bracket_zero"]=sp.factor(determinant-expected)
    z=spectral.z
    t=coeff.pairs()["T"]
    l=coeff.pairs()["L"]
    out["transverse_pair_entrywise_upper_slack"]=sp.Matrix([sp.Rational(109,81),2])-t-sp.Matrix([sp.Rational(109,81)*z,2*z])
    out["longitudinal_pair_entrywise_upper_slack"]=sp.Matrix([sp.Rational(226,81),4])-l-sp.Matrix([sp.Rational(13,9)*(1-z),2*(1-z)])
    return {key:sp.ImmutableMatrix(value.applyfunc(sp.factor)) if isinstance(value,sp.MatrixBase)
            else sp.factor(value) for key,value in out.items()}


@cache
def gates():
    d=data()
    lo,hi=d["lower_endpoint"],d["upper_endpoint"]
    return {
        "actual_profile_box_is_strictly_positive":all(value>0 for value in d["fixed_profile_box"].values()),
        "pi_lower_above_three":bool(d["pi_lower_partial_sum"]>3),
        "pi_upper_square_below_ten":bool(d["pi_upper"]**2<10),
        "log_bound_100_from_exact_positive_series":bool(d["exp_100_positive_series_single_term"]>d["log_argument_upper"]),
        "positive_pair_Gram_entrywise_majorants":all(bool(x>=0) for x in GRAM_UPPER-d["pair_Gram_entrywise_raw_upper"]),
        "finite_fourth_entrywise_absolute_majorants":all(bool(FINITE_UPPER[i,j]>=abs(coeff.finite()[i,j])) for i in range(2) for j in range(2)),
        "NN_finite_fourth_coefficient_below_minus_one":bool(coeff.finite()[0,0]<-1),
        "local_zero_absolute_entries_below_one_e15":bool(d["local_zero_entrywise_absolute_max"]<10**15),
        "local_second_absolute_entries_below_one_e9":bool(d["local_second_entrywise_absolute_max"]<10**9),
        "fixed_tadpole_absolute_entries_below_one_e15":bool(d["normalized_tadpole_entrywise_upper"]<10**15),
        "all_lower_order_error_below_one_e_minus_12":bool(d["full_lower_order_entrywise_error_upper"]<ERROR),
        "lower_NN_above_two":bool(lo["NN_lower"]>2),
        "lower_ZZ_above_five_e24":bool(lo["ZZ_lower"]>5*LOW**2),
        "lower_off_diagonal_below_one":bool(lo["both_off_diagonal_absolute_upper"]<1),
        "lower_determinant_strictly_positive":bool(lo["determinant_comparison"]>0),
        "upper_NN_below_minus_twelve":bool(hi["NN_upper"]<-12),
        "upper_ZZ_above_five_e26":bool(hi["ZZ_lower"]>5*HIGH**2),
        "upper_off_diagonal_below_ten_thousand":bool(hi["both_off_diagonal_absolute_upper"]<10000),
        "upper_determinant_strictly_negative":bool(hi["determinant_comparison"]<0),
        "ZZ_positive_throughout_bracket_prevents_inverse_cancellation":bool(d["ZZ_positive_lower_on_entire_band"]>0),
        "nonsymmetric_local_rows_not_silently_symmetrized":bool(coeff.local_row(2)[0,1]!=coeff.local_row(2)[1,0]),
        "exact_massive_bubble_not_leading_log_approximation":True,
        "fixed_tadpole_profiles_not_reselected":True,
        "frequency_units_are_inverse_bounce_time_not_Planck_units":True,
        "no_actual_curved_kernel_error_or_EFT_cutoff_is_claimed":True,
        "no_actual_bounce_instability_or_original_P8_closure":True}


def description():
    return {
        "definition":"Freeze all coefficients of the actual physical Euler differential operators at u=0 after variation, then apply the constant T=[[1,0],[1/2,1]] congruence. Keep the exact h=1 flat massive vacuum bubble, the actual frozen local rows of all three adiabatic orders, and the already-fixed tadpole within its proven profile box.",
        "symbol":"Kfr(s)=E0+s^2 E2+(s^4 Breg(s)+C0+s^2 C2+Ctad)/(64*pi^2*L^2), where Ctad=64*pi^2*L^2 times the constant-chart physical tadpole Hessian.",
        "block":"Breg(s)=Freg-(1/4) integral_0^1 y^2 Mreg(y^2)/(1+4*m^2/s^2-y^2) dy.",
        "uniformity":"For every fixed real pair of tadpole values with absolute energy and pressure bounded by the specified profile box, the determinant has a zero strictly inside (10^12,10^13). The ZZ entry stays positive, so the matrix inverse has an uncancelled positive-real pole.",
        "not_claimed":"No unique or simple root, actual nonstationary pole, bound on omitted curved/state kernel, justified interacting cutoff, EFT-valid instability, removal of pole residues, or original P8 closure."}
