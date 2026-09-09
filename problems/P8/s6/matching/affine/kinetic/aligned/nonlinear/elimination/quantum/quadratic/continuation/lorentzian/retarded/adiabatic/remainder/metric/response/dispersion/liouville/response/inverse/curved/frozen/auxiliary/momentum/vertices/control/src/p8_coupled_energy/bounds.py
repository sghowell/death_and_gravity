"""Explicit rational component and energy bounds on the two scalar charts."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_vertices import majorant as vertices

from . import scalars as s

ZMAX=sp.Rational(1,10**4)


def exact_expression(value,variables):
    if type(value) is bool or not isinstance(value,(int,sp.Expr)):
        raise ValueError("Require an exact rational expression")
    value=sp.sympify(value)
    if value.has(sp.Float,sp.oo,-sp.oo,sp.zoo,sp.nan) or value.free_symbols-set(variables):
        raise ValueError("Require finite exact coefficients in the declared variables")
    numerator,denominator=sp.fraction(sp.cancel(value))
    try:
        sp.Poly(numerator,*variables,domain=sp.QQ)
        sp.Poly(denominator,*variables,domain=sp.QQ)
    except (sp.PolynomialError,sp.CoercionFailed) as exc:
        raise ValueError("Require a rational function over QQ") from exc
    return value


def time_upper(value):
    value=exact_expression(value,(model.u,))
    numerator,denominator=sp.fraction(sp.cancel(value))
    num,den=sp.Poly(numerator,model.u),sp.Poly(denominator,model.u)
    if den.nth(0)<=0 or any(coefficient<0 or degree[0]%2 for degree,coefficient in den.terms()):
        raise ValueError("No positive-even time-denominator certificate")
    upper=sum(abs(coefficient)*sp.Rational(1,2)**degree[0] for degree,coefficient in num.terms())/den.nth(0)
    return {"absolute_half_interval_upper":upper,"denominator_constant_lower":den.nth(0)}


@cache
def backgrounds():
    actual=s.actual_background_jets()
    proofs={symbol:time_upper(value) for symbol,value in actual.items()}
    bounds={symbol:proof["absolute_half_interval_upper"] for symbol,proof in proofs.items()}
    # Independent elementary H, Lambda and Theta enclosures on I;
    # Je also inherits the actual second lapse derivative bound.
    bounds[s.H]=min(bounds[s.H],sp.Integer(2))
    bounds[s.ell]=min(bounds[s.ell],sp.Rational(1,10))
    bounds[s.theta]=min(bounds[s.theta],sp.Integer(3))
    bounds[s.lam]=min(bounds[s.lam],sp.Rational(1,2))
    bounds[s.Je]=min(bounds[s.Je],vertices.coefficient_bounds()["common_N_derivative_bounds"][2]/2)
    if bounds[s.Je]>=50:
        raise ValueError("The chosen finite-q denominator estimate needs Je<50")
    bounds[s.z]=ZMAX
    return {"actual_jets":actual,"positive_even_denominator_proofs":proofs,"absolute_bounds":bounds}


def coefficient_bound(value,chart):
    s.chart(chart)
    variables=s.BASE+s.FIRST+s.SECOND+(s.z,)
    value=exact_expression(value,variables)
    numerator,denominator=sp.fraction(sp.cancel(value))
    bounds=backgrounds()["absolute_bounds"]
    num=sp.Poly(numerator,*variables)
    top=sum(abs(coefficient)*sp.prod(bounds[variable]**degree
            for variable,degree in zip(variables,powers)) for powers,coefficient in num.terms())
    Je,ell,lam,z=s.Je,s.ell,s.lam,s.z
    D=lam**2-(Je+ell**2*lam**2/2)*z
    R=lam**2-Je*z
    primitives=((Je,sp.Rational(1,40)),(ell,sp.Rational(1,50)),
                (s.theta if chart=="unitary" else lam,sp.Rational(1,4)))
    if chart=="gamma":
        primitives+=((D,sp.Rational(1,32)),(R,sp.Rational(1,32)))
    constant,factors=sp.factor_list(denominator)
    lower=abs(constant)
    proof=[]
    for factor,power in factors:
        for primitive,minimum in primitives:
            ratio=sp.cancel(factor/primitive)
            if not ratio.free_symbols and ratio!=0:
                lower*=(abs(ratio)*minimum)**power
                proof.append({"factor":factor,"power":power,"absolute_lower":abs(ratio)*minimum})
                break
        else:
            raise ValueError("Uncertified chart denominator: "+str(factor))
    if lower<=0:
        raise ValueError("Nonpositive rational denominator bound")
    return {"absolute_upper":sp.factor(top/lower),"numerator_box_bound":top,
            "denominator_absolute_lower":lower,"denominator_factor_proofs":tuple(proof)}


def matrix_bound(matrix,chart):
    s.chart(chart)
    if not isinstance(matrix,sp.MatrixBase) or matrix.shape!=(2,2):
        raise ValueError("Require a two-component scalar matrix")
    rows=tuple(coefficient_bound(value,chart) for value in matrix)
    return {"entry_proofs":rows,"operator_upper":sum(row["absolute_upper"] for row in rows)}


def chart_bounds(chart):
    return _chart_bounds(s.chart(chart))


@cache
def _chart_bounds(chart):
    data=s.data(chart)
    names=("potential_remainder","antisymmetric_mixing","beta_leading","beta_remainder",
           "energy_velocity_coefficient","energy_principal_coefficient",
           "energy_remainder_coefficient","energy_mixing_coefficient")
    result={name:matrix_bound(data[name],chart) for name in names}
    upper=lambda name:result[name]["operator_upper"]
    qmin=max(sp.Integer(10**8),2000*upper("potential_remainder"))
    growth=(1000*upper("energy_velocity_coefficient")
            +2000*(upper("energy_principal_coefficient")+upper("energy_remainder_coefficient")/qmin)
            +sp.Rational(1,5)*upper("energy_mixing_coefficient"))
    return {"coefficient_bounds":result,"sufficient_q_lower":qmin,
            "energy_logarithmic_growth_upper":growth,
            "kinetic_lower":sp.Rational(1,1000),"kinetic_upper":sp.Integer(10**4),
            "principal_gradient_lower":sp.Rational(1,1000),"principal_gradient_upper":sp.Integer(10**4),
            "full_potential_lower_over_q":sp.Rational(1,2000)}


@cache
def elementary_margins():
    return {"gamma_Lambda_absolute_lower_margin":sp.Rational(1231,4913)-sp.Rational(1,4),
            "unitary_Theta_absolute_lower_margin":sp.Rational(864,3125)-sp.Rational(1,4),
            "gamma_D_denominator_margin":sp.Rational(1,16)-(50+sp.Rational(1,800))*ZMAX-sp.Rational(1,32),
            "gamma_R_denominator_margin":sp.Rational(1,16)-50*ZMAX-sp.Rational(1,32),
            "matter_density_lower_margin":sp.Rational(1,10)*(sp.Rational(4,5))**6-sp.Rational(1,50),
            "completed_square_eigenvalue_margin":sp.Rational(1,540)-sp.Rational(1,1000),
            "unitary_kinetic_trace_upper_margin":10**4-(1600+sp.Rational(1,25)+1),
            "gamma_kinetic_trace_upper_margin":10**4-32*(100+sp.Rational(1,400))-2,
            "gamma_window_chart_margin":sp.Rational(1,4)-sp.Rational(19,80)-sp.Rational(1,160),
            "unitary_window_chart_margin":sp.Rational(19,80)-sp.Rational(1,160)-sp.Rational(9,40)}
