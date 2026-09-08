"""Physical bilocal current vertices and their dimensional leading symbol."""
from functools import cache

import sympy as sp
from p8_vector_metric_dispersion import principal, spectral
from p8_vector_metric_local import canonical, jets, variation

h=spectral.h
k, a, p = sp.symbols("comoving_k frozen_scale physical_momentum",positive=True)
D=sp.Symbol("spatial_dimension")
At,Bt,As,Bs,bt,bs,Ut,Us=sp.symbols("A_t B_t A_s B_s pump_t pump_s U_t U_s",real=True)


def clean(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.factor)) if isinstance(value,sp.MatrixBase) else sp.simplify(sp.factor(value))


def high(value):
    return _high(jets.kind(value))


@cache
def _high(value):
    weights=canonical.data(value)["weights"]
    sub={jets.D:D,jets.z:1,jets.alpha[0]:4/(9*h),jets.beta[0]:28/(81*h)}
    A=sp.ImmutableMatrix([weights[key][0].subs(sub,simultaneous=True) for key in ("N","Z")])
    B=sp.ImmutableMatrix([weights[key][1].subs(sub,simultaneous=True) for key in ("N","Z")])
    return {"A":A,"B":B,"pair":clean(A-B)}


@cache
def polynomial():
    return sp.Poly(sp.expand(
        At*As*((-sp.I*k-bt)*(sp.I*k-bs))**2
        +Bt*Bs*(k*k+Ut)*(k*k+Us)
        +At*Bs*(k*k+Us)*(-sp.I*k-bt)**2
        +Bt*As*(k*k+Ut)*(sp.I*k-bs)**2),k)


@cache
def data():
    left,right=sp.symbols("positive_h_t positive_h_s",positive=True)
    pair_t=high("L")["pair"].subs({D:3,h:left})
    pair_s=high("L")["pair"].subs({D:3,h:right})
    return {"physical_dimensional_transverse":high("T"),
            "physical_dimensional_longitudinal":high("L"),
            "generic_acoustic_vertex_polynomial_coefficients":tuple(polynomial().nth(j) for j in range(5)),
            "longitudinal_leading_pair_bimatrix":pair_t*pair_s.T,
            "normalized_four_primitive_inverse_lag_bimatrix":pair_t*pair_s.T/4,
            "flat_fourth_order_pole":principal.pole(),
            "normalization":"Q=64*pi^2 times the unscaled Gaussian physical current response. For D=3 the normalized radial leading retarded kernel is 8*b(t)b(s)^T*k^4*sin(2*k*Delta_sigma)/(a(t)^4*a(s))."}


@cache
def checks():
    poly=polynomial()
    pair=high("L")["pair"]
    transverse=high("T")["pair"]
    out={"generic_leading_bilocal_pair":clean(poly.nth(4)-(At-Bt)*(As-Bs)),
         "generic_subleading_connection_not_dropped":clean(poly.nth(3)-2*sp.I*(bs*(At-Bt)*As-bt*At*(As-Bs))),
         "transverse_high_pair_full_dimension":clean(transverse-sp.Matrix([0,2*(D-3)])),
         "longitudinal_high_pair_full_dimension":clean(pair-sp.Matrix([sp.Rational(64,81)/h,2*(D-1)])),
         "transverse_no_physical_highest_log":clean(transverse.subs(D,3)),
         "transverse_no_first_dimensional_highest_finite_jet":clean(
             sp.diff((D-1)*transverse*transverse.T,D).subs(D,3)),
         "longitudinal_first_dimensional_pair_jet":clean(sp.diff(pair,D)-sp.Matrix([0,2])),
         "physical_rank_one_pole_matches_existing_full_local_pole":clean(
             pair.subs(D,3)*pair.subs(D,3).T/8-principal.pole()),
         "proper_momentum_scale_cancels_for_every_dimension":sp.simplify(
             a**(-D-2)*(a*p)**(D+1)*a-p**(D+1))}
    tau=sp.Symbol("lag",real=True)
    # a(t)=a0+a1*tau+..., sigma(t)-sigma(s)=tau/a0-a1*tau^2/(2a0^2)+...
    a0=sp.Symbol("a0",positive=True)
    a1=sp.Symbol("a1",real=True)
    geometry=(1+a1*tau/a0)**(-4)*(1-a1*tau/(2*a0))**(-5)
    out["leading_geometry_diagonal"] = geometry.subs(tau,0)-1
    out["leading_geometry_first_lag_jet"] = clean(sp.diff(geometry,tau).subs(tau,0)+3*a1/(2*a0))
    value=pair.subs(D,3)
    Avec=high("L")["A"].subs(D,3)
    connection=clean(2*(value*Avec.T-Avec*value.T))
    out["connection_diagonal_is_antisymmetric"] = clean(connection+connection.T)
    out["connection_diagonal_explicit"] = clean(connection-sp.Matrix([
        [0,8+32/(27*h)],[-8-32/(27*h),0]]))
    lag=sp.Symbol("positive_lag",positive=True)
    out["radial_leading_Abel_limit"] = clean(
        sp.im(sp.factorial(4)/(sp.Symbol("epsilon",positive=True)-2*sp.I*lag)**5)).subs(
            sp.Symbol("epsilon",positive=True),0)-sp.Rational(3,4)/lag**5
    out["radial_subleading_Abel_limit"] = clean(
        sp.re(sp.factorial(3)/(sp.Symbol("epsilon",positive=True)-2*sp.I*lag)**4)).subs(
            sp.Symbol("epsilon",positive=True),0)-sp.Rational(3,8)/lag**4
    out["normalized_highest_kernel_four_derivative_bridge"] = clean(
        8*sp.Rational(3,4)/lag**5-sp.diff(1/(4*lag),lag,4))
    out["normalized_log_matches_exact_massive_symbol"] = clean(
        value*value.T/4-2*principal.pole())
    return out


@cache
def adiabatic_principal_checks():
    out={}
    for sector in ("T","L"):
        weights=canonical.data(sector)["weights"]
        pairs=sp.ImmutableMatrix([weights[key][0]-weights[key][1] for key in ("N","Z")])
        for i,output in enumerate(("N","Z")):
            coefficient=variation.data(sector)["coefficients"][output][2]
            for j,row in enumerate((jets.n,jets.v)):
                actual=sp.expand(coefficient).coeff(row[4])
                out[sector+output+str(j)+"_full_dimensional_all_momentum_fourth_Taylor_contact"]=sp.factor(
                    actual-pairs[i]*pairs[j]/32)
    return out
