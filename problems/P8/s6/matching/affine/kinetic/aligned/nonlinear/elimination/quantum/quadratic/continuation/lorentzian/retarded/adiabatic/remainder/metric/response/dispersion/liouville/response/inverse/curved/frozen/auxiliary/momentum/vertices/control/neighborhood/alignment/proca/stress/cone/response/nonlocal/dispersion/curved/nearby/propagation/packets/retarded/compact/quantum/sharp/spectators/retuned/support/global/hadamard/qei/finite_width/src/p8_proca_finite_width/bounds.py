"""Explicit finite-width reference QEI majorant, with separate errors."""
from fractions import Fraction as F
from functools import cache
from math import comb

import sympy as sp
from flint import ctx
from p8_proca_global_hadamard import canonical
from p8_proca_global_support import model

from .compact import ceil, transfer_bounds
from .compact import matrix_norm as rational_matrix_norm
from .intervals import build, matrix_norm


def ibp_coefficients(L,D,order=3):
    if isinstance(order,bool) or not isinstance(order,int) or order<0:
        raise ValueError("Require a nonnegative integer integration-by-parts order")
    if not isinstance(L,(list,tuple)) or not isinstance(D,(list,tuple)):
        raise TypeError("Require finite derivative-majorant sequences")
    if any(isinstance(v,bool) or not isinstance(v,int) or v<0 for v in list(L)+list(D)):
        raise ValueError("Require nonnegative exact integer derivative majorants")
    if len(L)<=order or len(D)<=order:
        raise ValueError("Insufficient nonnegative derivative majorants")
    c={(0,j,r):comb(j,r)*L[j-r] for j in range(order+1) for r in range(j+1)}
    for n in range(order):
        for j in range(order-n):
            for r in range(n+j+2):
                c[n+1,j,r]=sum(comb(j+1,t)*c.get((n,t,r),0)*D[j+1-t] for t in range(j+2))
    return [c[order,0,r] for r in range(order+1)]

@cache
@ctx.workprec(160)
def data():
    K=32
    h=F(1,100)
    j=build()
    T=j["T"]
    bc=j["B_coefficient_derivative_upper_integers"]
    B=[ceil(sum(F(row[d],K**n) for n,row in bc.items())) for d in range(4)]
    omega=j["omega_derivative_upper_integers"]
    assert j["frequency_lower_gt_99_over_100"]
    assert F(99,100)*K>2*B[0]
    # s=alpha+k; Q=(A-i alpha I)^-1. Bounds Q^(j)<=D[j]/s.
    D=[3]
    C=[ceil(F(omega[d])+F(B[d],K)) for d in range(4)]
    for n in range(1,4):
        D.append(D[0]*sum(comb(n,r)*C[r]*D[n-r] for r in range(1,n+1)))
    defect=ceil(sum(F(v,K**(n-7)) for n,v in j["defect_coefficient_upper_integers"].items()))
    # Norms of the physical observable rows in the moving V frame.
    u=model.u
    d=1+u*u
    ell=model.background()["ell"]
    f=T.matrix([[T.rational(-ell/d**3),T.rational(1/d**3)],[T.zero,T.zero]])
    space=T.matrix([[T.rational(-ell/d**5),T.rational(1/d**5)],[T.zero,T.zero]])
    freq=T.matrix([[j["omega"][0],T.zero],[T.zero,j["omega"][1]]])
    time_terms={0:T.mscale(T.mmul(f,freq),-sp.I),1:T.mderivative(f)}
    for n,b in j["B_coefficients"].items():
        time_terms[n+1]=T.madd(time_terms.get(n+1,T.mzero()),T.mmul(f,b))
    lt=[ceil(sum(F(matrix_norm(v,r),K**n) for n,v in time_terms.items())) for r in range(4)]
    ls=[matrix_norm(space,r) for r in range(4)]
    tt,ts=ibp_coefficients(lt,D),ibp_coefficients(ls,D)
    transfer=transfer_bounds(h)
    E=transfer["energy_upper"]/transfer["energy_lower"]
    old=canonical.data()
    actual={n:v.subs(model.substitution(),simultaneous=True) for n,v in old["complete_Laurent_coefficients"].items()}
    b_can=ceil(sum(max(rational_matrix_norm(actual[n],h),rational_matrix_norm(actual[n].T,h))
        for n in (0,-1)))
    assert b_can<=transfer["remainder_norm_upper"]
    jnorm=ceil(rational_matrix_norm(actual[1],h))
    hnorm=ceil(F(3,2)*rational_matrix_norm(sp.Matrix([[model.background()["H"]]]),h))
    Ct=jnorm+transfer["remainder_norm_upper"]+hnorm
    assert Ct<=17
    assert matrix_norm(j["V"])<=2 and matrix_norm(j["Vi"])<=2
    # Initial negative-graph configuration factor Vi(0) B(0,k)^(-1/2) has norm<2.
    assert F(11,10)**2/F(3,4)<4
    # All constants below multiply hbar/kappa. Exponentials are exact symbolic values.
    H=sp.Rational(h.numerator,h.denominator)
    ER=sp.Rational(E.numerator,E.denominator)
    rate=transfer["norm_log_rate"]
    # Conservative sqrt(ER)<14; full transfer within [-h,h] uses duration<=2h.
    assert E<14**2
    Ce=14*sp.exp(2*rate*H)*(sp.Rational(1200002,K**2)+4*H*defect*sp.exp(B[0]*H))
    low=4*(Ct**2+1)*ER*sp.exp(2*rate*H)/sp.pi**2*(sp.Rational(K**4,4)+sp.Rational(K**3,3))
    # N=3 => integral_{K}^inf integral_0^inf k^3/(alpha+k)^6 = 1/(5K).
    # Include the factor 2 from the actual = approximate + error split here too.
    app=sp.exp(2*B[0]*H)/(5*sp.pi**3*K)
    # Error receives the extra factor 2 from |actual|^2<=2|approx|^2+2|error|^2.
    error=(Ct**2+1)*Ce**2/(16*sp.pi**2*K**8)
    return {"half_width":H,"high_frequency_threshold":K,"configuration_remainder_derivative_bounds":B,
        "frequency_derivative_bounds":omega,"inverse_derivative_bounds":D,
        "time_observable_row_derivative_bounds":lt,"spatial_observable_row_derivative_bounds":ls,
        "time_IBP_Sobolev_L1_coefficients":tt,"spatial_IBP_Sobolev_L1_coefficients":ts,
        "scaled_Riccati_defect_order_seven_constant":defect,
        "full_phase_energy_ratio":ER,"full_phase_log_norm_rate":rate,
        "physical_time_derivative_observable_constant":Ct,
        "normalized_packet_mode_error_constant":Ce,
        "low_band_L2_squared_coefficient":low,
        "approximate_high_band_IBP_squared_coefficient":app,
        "high_band_error_L2_squared_coefficient":error,
        "final_bound_description":"B_ref <= (hbar/kappa)*[(C_low+C_error)*||g||_2^2 + C_app*((sum T_j||g^(j)||_1)^2+(sum S_j||g^(j)||_1)^2)]"}

if __name__=="__main__":
    for k,v in data().items():
        print(k,v,flush=True)
