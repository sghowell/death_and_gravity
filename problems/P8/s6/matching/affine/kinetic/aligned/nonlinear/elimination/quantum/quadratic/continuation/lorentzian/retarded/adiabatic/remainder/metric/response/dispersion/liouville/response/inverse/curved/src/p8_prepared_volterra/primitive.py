"""Exact time-primitive identities underlying the no-loss Volterra class."""
from functools import cache

import sympy as sp

t,s=sp.symbols("output_time source_time",real=True)
lag=sp.Symbol("positive_lag",positive=True)


def order(value):
    if type(value) is not int or value not in (0,1,2,3,4):
        raise ValueError("Require native derivative order zero through four")
    return value


def local_kernel(value):
    value=order(value)
    return _local_kernel(value)


@cache
def _local_kernel(value):
    c=sp.Function("local_coefficient")(s)
    start=1 if value==4 else 0
    return sp.expand(sum((-1)**j*sp.binomial(value,j)*sp.diff(c,s,j)
        *(t-s)**(3-value+j)/sp.factorial(3-value+j) for j in range(start,value+1)))


def primitive_shape(value):
    value=order(value)
    if value==4:
        return 1/(24*lag)
    n=value+1
    # d^4 of lag^(4-n) log(lag) has only lag^-n; fix its exact scalar.
    trial=lag**(4-n)*sp.log(lag)
    coefficient=sp.simplify(sp.diff(trial,lag,4)*lag**n)
    return trial/coefficient


@cache
def checks():
    c=sp.Function("local_coefficient")(s)
    out={}
    for j in range(4):
        out["local_order_"+str(j)+"_independent_IBP"]=sp.simplify(
            local_kernel(j)-(-1)**j*sp.diff((t-s)**3*c/6,s,j))
    # Fourth derivative: distributional endpoint c(t)*f(t), plus this ordinary kernel.
    out["fourth_local_regular_IBP_kernel"]=sp.simplify(
        local_kernel(4)-sp.diff((t-s)**3*c/6,s,4))
    for j in range(5):
        out["singular_power_"+str(j+1)+"_four_primitive"]=sp.simplify(
            sp.diff(primitive_shape(j),lag,4)-lag**(-j-1))
    out["fourth_local_instantaneous_coefficient"]=sp.diff(lag**3/6,lag,3)-1
    out["fourth_primitive_subleading_is_log_not_inverse_lag"]=sp.simplify(
        primitive_shape(3)+sp.log(lag)/6)
    f=sp.Function("coefficient")
    # Derivatives of a multiplier around a convolution retain every cross term.
    b=sp.Function("convolution_kernel")(t-s)
    direct=sp.diff(f(t)*b,t,4)
    terms=sum(sp.binomial(4,j)*sp.diff(f(t),t,j)*sp.diff(b,t,4-j) for j in range(5))
    out["all_fourth_multiplier_commutators_retained"]=sp.simplify(direct-terms)
    T=sp.Symbol("short_interval",positive=True)
    out["log_kernel_small_interval_integral"]=sp.simplify(
        sp.integrate(1-sp.log(lag),(lag,0,T))-T*(2-sp.log(T)))
    return out


def description():
    return {"local_fourth_derivative":"I4(c f'''')=c(t)f(t)+integral local_kernel(4)(t,s)f(s)ds for prepared f.",
            "subleading_kernel_class":"For each diagonal derivative order n, norm((partial_t+partial_s)^n V(t,s))<=C_n*(1+abs(log(t-s))) on the compact causal triangle.",
            "uniform_C0_short_interval_bound":"C_0*T*(2-log(T)) for 0<T<=1, with C_0 finite but not numerically evaluated.",
            "inverse_route":"Rpre=Rinf(t) multiplication plus a uniformly L1 Volterra kernel; Rpre V remains a Volterra kernel with an integrable scalar majorant.",
            "weighted_norm_route":"Choose lambda>0 so integral_0^1 exp(-lambda*r)*w(r)dr<1, for the derived integrable majorant w. Dominated convergence gives such a lambda; the full finite-interval Neumann inverse follows."}
