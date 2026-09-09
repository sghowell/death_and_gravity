"""Actual minimal-Proca flat bubble, local coefficient and scalar analytic block."""
from functools import cache

import sympy as sp
from p8_proca_local_response import local
from p8_vector_metric_dispersion import spectral as generic
from p8_vector_metric_local import canonical, jets, variation

z,y=sp.symbols("ordinary_momentum_fraction ordinary_radial_fraction",positive=True)
p=sp.Symbol("squared_Laplace_frequency")
m=sp.Symbol("positive_fixed_mass",positive=True)
d=sp.Symbol("radial_denominator",positive=True)
ZERO={field:sp.Integer(0) for row in (jets.alpha,jets.beta,jets.alpha2,jets.beta2) for field in row}


def pair(sector):
    sector=jets.kind(sector)
    return sp.ImmutableMatrix([0,2*(1-z)] if sector=="T" else [0,2*(1+z)])


def positive_axis(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact nonnegative squared-frequency over mass-squared ratio")
    value=sp.Rational(value)
    if value<0:
        raise ValueError("Require p/m^2>=0 on the positive axis")
    H=sp.Integer(4) if value==0 else data()["closed_scalar_H_in_d"].subs(d,1+4/value)
    return {"p_over_mass_squared":value,"H":H,"scalar_block":-H,
            "range_inverse":-1/H,"H_continuous_lower":sp.Integer(4),
            "range_inverse_absolute_upper":sp.Rational(1,4)}


@cache
def data():
    matrix=sp.ImmutableMatrix(2*pair("T")*pair("T").T+pair("L")*pair("L").T)
    polynomial=3-2*z+3*z*z
    weight=y*y*polynomial.subs(z,y*y)
    H=4+3*generic.moment(1,d)-2*generic.moment(2,d)+3*generic.moment(3,d)
    closed=sp.Rational(16,15)+d-3*d*d+sp.sqrt(d)*(3-2*d+3*d*d)*sp.atanh(1/sp.sqrt(d))
    coeffs={j:sp.factor((-1)**(j-1)*sp.integrate(weight*(1-y*y)**(j-1),(y,0,1))/4**j)
            for j in range(1,5)}
    return {"actual_physical_pair_Gram_matrix":matrix,"positive_pair_polynomial":polynomial,
            "positive_radial_weight":weight,"positive_sign_scalar_H_in_d":H,
            "closed_scalar_H_in_d":closed,"d_of_squared_Laplace_frequency":1+4*m*m/p,
            "analytic_scalar_H_integrand":p*weight/(4*m*m+p*(1-y*y)),
            "normalized_local_plus_thrice_subtracted_scalar_block":-closed.subs(d,1+4*m*m/p),
            "physical_two_source_block":sp.ImmutableMatrix(sp.diag(0,-closed.subs(d,1+4*m*m/p))),
            "first_sheet_cut":"p not in (-infinity,-4m^2]; p=0 is a removable point of the closed formula",
            "low_frequency_H_coefficients_in_p_over_mass_squared":coeffs,
            "H_at_zero":sp.Integer(4),"H_at_threshold":sp.Rational(16,15),
            "H_slope_at_zero":sp.Rational(9,35)/m**2,
            "high_frequency_H":"2 log(p/m^2)-14/15+O((m^2/p) log(p/m^2)) on the first sheet",
            "scope":"Isolated flat-vacuum reference bubble plus the new finite fourth-order local coefficient. Not the full actual curved selected-state or tree-plus-loop inverse."}


def contact(sector):
    return _contact(jets.kind(sector))


@cache
def _contact(sector):
    N,V,w,mass=sp.symbols("minimal_lapse minimal_logscale frequency fixed_mass",positive=True)
    q=w*w-mass*mass
    qp=q*sp.exp(-2*V)
    g2=sp.exp(V)/N if sector=="T" else sp.exp(3*V)*mass*mass*qp/(N*(qp+mass*mass))
    frequency2=N*N*(qp+mass*mass)
    point={N:1,V:0}
    g0=g2.subs(point)
    expectation=w*(g0/g2+g2*frequency2/(g0*w*w))/4
    seagull=sp.Matrix([[-sp.factor(sp.diff(expectation,left,right).subs(point)) for right in (N,V)] for left in (N,V)])
    static=sp.Matrix([[sp.factor(sp.diff(-sp.sqrt(frequency2)/2,left,right).subs(point)) for right in (N,V)] for left in (N,V)])
    pairs=pair(sector).subs(z,1-mass*mass/w**2)
    return {"actual_fixed_canonical_contact":sp.ImmutableMatrix(seagull),
            "static_vacuum_Hessian":sp.ImmutableMatrix(static),
            "zero_frequency_bubble":sp.ImmutableMatrix(w*pairs*pairs.T/8),
            "static_contact_identity":sp.ImmutableMatrix((seagull+w*pairs*pairs.T/8-static).applyfunc(sp.factor))}


@cache
def checks():
    out={}
    for sector in ("T","L"):
        weights=canonical.data(sector)["weights"]
        actual=sp.Matrix([sp.factor((weights[label][0]-weights[label][1]).xreplace(ZERO).subs({jets.D:3,jets.z:z}))
                          for label in ("N","Z")])
        out[sector+"_new_pair_from_actual_minimal_Hamiltonian"]=sp.ImmutableMatrix((actual-pair(sector)).applyfunc(sp.factor))
        out[sector+"_new_static_contact_not_dropped"]=contact(sector)["static_contact_identity"]
        frozen={field:0 for field in jets.H}
        frozen.update(ZERO)
        frozen.update({jets.D:3,jets.z:z})
        for order in (1,2):
            for i,output in enumerate(("N","Z")):
                value=variation.data(sector)["coefficients"][output][order]
                for j,row in enumerate((jets.n,jets.v)):
                    actual=sp.expand(value).coeff(row[2*order]).subs(frozen,simultaneous=True)
                    expected=(-1)**order*pair(sector)[i]*pair(sector)[j]/(2*4**order)
                    out[sector+f"_new_adiabatic_frequency_Taylor_{output}_{j}_{2*order}"]=sp.factor(actual-expected)
    item=data()
    matrix=item["actual_physical_pair_Gram_matrix"]
    out["new_rank_one_pair_Gram"]=sp.ImmutableMatrix((matrix-sp.diag(0,4*item["positive_pair_polynomial"])).applyfunc(sp.factor))
    out["pair_weight_strict_positive_square"]=sp.factor(item["positive_pair_polynomial"]-3*(z-sp.Rational(1,3))**2-sp.Rational(8,3))
    out["new_three_moment_closed_scalar_block"]=sp.simplify(item["positive_sign_scalar_H_in_d"]-item["closed_scalar_H_in_d"])
    actual=sp.Matrix([[sp.factor(local.operator(output,2).coeff(row[4])) for row in (local.n,local.v)] for output in ("N","Z")])
    out["actual_new_curved_local_fourth_coefficient"]=sp.ImmutableMatrix(actual-sp.diag(0,-4))
    out["threshold_value_by_independent_radial_integral"]=sp.factor(4-sp.integrate(item["positive_radial_weight"]/y**2,(y,0,1))-sp.Rational(16,15))
    out["static_slope_by_independent_radial_integral"]=sp.factor(
        sp.integrate(item["positive_radial_weight"],(y,0,1))/(4*m*m)-item["H_slope_at_zero"])
    eps,L=sp.symbols("positive_inverse_frequency large_log_p_over_mass_squared",positive=True)
    substituted=item["closed_scalar_H_in_d"].subs(sp.atanh(1/sp.sqrt(d)),L/2+sp.log((1+sp.sqrt(1+eps))/2)).subs(d,1+eps)
    out["uniform_large_circle_leading_constant"]=sp.simplify(substituted.subs(eps,0)-(2*L-sp.Rational(14,15)))
    out["spurious_large_frequency_truncation_root"]=2*sp.Rational(7,15)-sp.Rational(14,15)
    out["same_three_subtracted_denominator_and_radial_measure"]=generic.checks()["radial_measure_and_Kubo_normalization"]
    out["same_exact_three_subtraction_identity"]=generic.checks()["three_subtracted_retarded_denominator"]
    a,b=sp.symbols("real_squared_frequency positive_imaginary_squared_frequency",real=True)
    integrand=item["analytic_scalar_H_integrand"]
    expected=4*m*m*b*item["positive_radial_weight"]/((4*m*m+a*(1-y*y))**2+b*b*(1-y*y)**2)
    out["strict_half_plane_imaginary_sign_integrand"]=sp.factor(sp.im(integrand.subs(p,a+sp.I*b))-expected)
    out["strict_real_gap_monotonicity_integrand"]=sp.factor(sp.diff(integrand,p)
        -4*m*m*item["positive_radial_weight"]/(4*m*m+p*(1-y*y))**2)
    B=sp.Symbol("nonzero_scalar_block",nonzero=True)
    projector=sp.diag(0,1)
    out["range_inverse_only_not_full_two_source_inverse"]=sp.diag(0,B)*sp.diag(0,1/B)-projector
    return out
