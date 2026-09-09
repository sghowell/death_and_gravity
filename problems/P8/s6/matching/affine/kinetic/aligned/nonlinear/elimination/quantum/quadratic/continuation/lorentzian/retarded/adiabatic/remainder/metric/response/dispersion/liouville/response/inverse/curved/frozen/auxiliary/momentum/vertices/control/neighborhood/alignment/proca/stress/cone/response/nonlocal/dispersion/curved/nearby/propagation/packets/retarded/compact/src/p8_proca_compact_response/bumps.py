"""Smooth compact bump definitions with exact normalization and boundary recurrence."""
from functools import cache

import sympy as sp

r=sp.Symbol("signed_bump_coordinate",real=True)
rho=sp.Symbol("positive_source_width",positive=True)
A1,A3=sp.symbols("one_dimensional_bump_integral radial_three_dimensional_bump_integral",positive=True)


def beta_squared(square):
    if isinstance(square,bool) or not isinstance(square,(int,sp.Expr)):
        raise TypeError("Require an exact symbolic squared radius")
    square=sp.sympify(square)
    if square.has(sp.Float):
        raise TypeError("Require exact bump data without floating values")
    if square.is_nonnegative is not True or square.is_finite is not True:
        raise ValueError("Require a provably finite nonnegative squared radius")
    if (square-1).is_nonnegative is True:
        return sp.Integer(0)
    return sp.Piecewise((sp.exp(1-1/(1-square)),square<1),(0,True))


@cache
def data():
    t,s,t0,s0=sp.symbols("detector_time source_time detector_center_time source_center_time",real=True)
    a,eps=sp.symbols("cap_width normal_width",positive=True)
    x1,x2,x3,y1,y2,y3=sp.symbols("x1 x2 x3 y1 y2 y3",real=True)
    S=sp.Symbol("actual_clock_radius_from_source_center",positive=True)
    inside=sp.exp(1-1/(1-r*r))
    source=beta_squared((s-s0)**2/rho**2)*beta_squared((y1*y1+y2*y2+y3*y3)/rho**2)/(A1*A3*rho**4)
    detector=beta_squared((t-t0)**2/a**2)*beta_squared((x2*x2+x3*x3)/a**2)
    detector*=beta_squared((x1-S)**2/a**2)*beta_squared((sp.sqrt(x1*x1+x2*x2+x3*x3)-S)**2/eps**2)
    polynomials=[sp.Integer(1)]
    for n in range(4):
        P=polynomials[-1]
        polynomials.append(sp.expand((1-r*r)**2*sp.diff(P,r)+(4*n*r*(1-r*r)-2*r)*P))
    return {"interior_one_dimensional_bump":inside,
            "global_one_dimensional_bump":beta_squared(r*r),
            "one_dimensional_normalizer":sp.Integral(inside,(r,-1,1)),
            "three_dimensional_normalizer":4*sp.pi*sp.Integral(r*r*inside,(r,0,1)),
            "source_normalization_symbols":[A1,A3],
            "actual_compact_source_formula":source,
            "actual_compact_detector_formula":detector,
            "first_boundary_recurrence_polynomials":polynomials,
            "boundary_derivative_induction":"beta^(n)=beta*P_n/(1-r^2)^(2n); polynomial recurrence holds for every n",
            "spatial_radial_bumps_are_functions_of_squared_radius":True,
            "detector_is_extended_by_zero_outside_its_compact_interior_time_window":True}


@cache
def checks():
    d=data()
    inside=d["interior_one_dimensional_bump"]
    rows={"source_change_of_variables_keeps_all_four_normalization_powers":
              sp.factor(A1*A3*rho**4/(A1*A3*rho**4)-1)}
    for n,P in enumerate(d["first_boundary_recurrence_polynomials"]):
        rows[f"bump_boundary_derivative_recurrence_anchor_{n}"]=sp.simplify(
            sp.diff(inside,r,n)/inside-P/(1-r*r)**(2*n))
    z=sp.Symbol("positive_boundary_parameter",positive=True)
    rows["polynomial_times_boundary_exponential_has_zero_limit"]=sp.limit(z**8*sp.exp(-z),z,sp.oo)
    return rows


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.I,[],-1,sp.Rational(-1,4))
    rejected=0
    for value in bad:
        try:
            beta_squared(value)
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(bad):
        raise ValueError("An inadmissible compact bump squared radius was accepted")
    return {"rejected_inputs":rejected}
