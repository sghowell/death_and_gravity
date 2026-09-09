"""Actual finite ordinary-Proca metric response in the fixed covariant scheme."""
from functools import cache

import sympy as sp
from p8_constant_proca import quantum
from p8_vector_metric_local import counterterms, jets, radial

u,n,v=jets.u,jets.n,jets.v
H=4*u/(1+u*u)
ZERO={field:sp.Integer(0) for row in (jets.alpha,jets.beta,jets.alpha2,jets.beta2) for field in row}


def actual(value,extra=()):
    mapping={field:sp.diff(H,u,j) for j,field in enumerate(jets.H)}
    fields=n+v+tuple(field for row in extra for field in row)
    poly=sp.Poly(sp.expand(value.xreplace(ZERO)),*fields)
    return sp.Add(*(sp.factor(c.subs(mapping,simultaneous=True))*sp.prod(x**k for x,k in zip(fields,powers))
                    for powers,c in poly.terms()))


def clean(value,extra=()):
    fields=n+v+tuple(field for row in extra for field in row)
    poly=sp.Poly(sp.expand(value),*fields)
    return sp.Add(*(sp.factor(c)*sp.prod(x**k for x,k in zip(fields,powers)) for powers,c in poly.terms()))


def time(value,extra=()):
    return clean(sp.diff(value,u)+sum(sp.diff(value,row[j])*row[j+1]
                 for row in (n,v)+tuple(extra) for j in range(len(row)-1)),extra)


def weighted(value,extra=()):
    return clean(time(value,extra)+3*H*value,extra)


def euler(density,row,extra=()):
    value=sp.Integer(0)
    for j,field in enumerate(row):
        term=sp.diff(density,field)
        for _ in range(j):
            term=-weighted(term,extra)
        value+=term
    return clean(value,extra)


def density(order):
    return _density(jets.order(order))


@cache
def _density(order):
    e,D=counterterms.e,jets.D
    geometry=counterterms.geometry()
    invariant={0:sp.Rational(5,2),1:5*geometry["R"]/3,2:-4*geometry["a4sc"]}[order]
    volume=1+e*(n[0]+D*v[0])+e**2*(D*n[0]*v[0]+D**2*v[0]**2/2)
    return actual(counterterms.keep(volume*invariant).coeff(e,2).subs(D,3))


def operator(output,order):
    return _operator(jets.output(output),jets.order(order))


@cache
def _operator(output,order):
    return euler(density(order),n if output=="N" else v)


def coefficients(output,order,source):
    output,order,source=jets.output(output),jets.order(order),jets.output(source)
    value=operator(output,order)
    row=n if source=="N" else v
    return tuple(sp.factor(value.coeff(field)) for field in row[:5])


def physical_operator(output,order):
    if type(output) is not str or output not in ("energy","pressure"):
        raise ValueError("Require physical energy or pressure output")
    order=jets.order(order)
    baseline=quantum.local_coefficients()[order][output]
    result=(-operator("N",order)-3*v[0]*baseline if output=="energy"
            else operator("Z",order)/3-(n[0]+3*v[0])*baseline)
    return clean(result)


@cache
def checks():
    out={}
    for order in range(3):
        for output in ("N","Z"):
            raw=radial.laurent(output,order)
            ct=counterterms.generic_operator(output,order).xreplace(ZERO)
            out[output+f"_actual_constant_mass_pole_{order}"]=sp.expand(
                raw["pole"].xreplace(ZERO)-ct.subs(jets.D,3))
            matched=raw["finite_MSbar_mu_m"].xreplace(ZERO)+2*sp.diff(ct,jets.D).subs(jets.D,3)
            out[output+f"_independent_ordinary_finite_heat_action_{order}"]=sp.expand(
                actual(matched)-operator(output,order))
            out[output+f"_no_source_derivative_above_four_{order}"]=sp.expand(
                operator(output,order)-sum(sp.diff(operator(output,order),f)*f for f in n[:5]+v[:5]))
        for component in ("energy","pressure"):
            out[component+f"_physical_normalization_constant_spatial_rescaling_{order}"]=sp.factor(
                physical_operator(component,order).coeff(v[0]))
    return out
