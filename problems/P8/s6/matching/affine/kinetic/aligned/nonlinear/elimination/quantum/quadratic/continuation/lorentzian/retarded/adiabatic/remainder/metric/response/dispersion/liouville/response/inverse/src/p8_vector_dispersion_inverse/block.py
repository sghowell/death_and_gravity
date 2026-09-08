"""First-sheet matrix data; not the full tree-plus-loop response symbol."""
from functools import cache

import sympy as sp
from p8_vector_metric_dispersion import principal, spectral

h, z, y, m = spectral.h, spectral.z, spectral.y, spectral.m
p_real, p_imag = sp.symbols("p_real p_imag", real=True)
x = sp.Symbol("subthreshold_fraction", real=True)


def clean(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.factor))


def positive_exact(value):
    if isinstance(value, (bool, float, str)) or value in (sp.true, sp.false):
        raise ValueError("Require a finite positive exact number")
    value = sp.sympify(value)
    if not isinstance(value, sp.Expr) or value.is_number is not True or value.is_real is not True or value.is_positive is not True:
        raise ValueError("Require a finite positive exact number")
    if value.has(sp.Float) or value in (sp.oo, sp.nan, sp.zoo):
        raise ValueError("Require a finite positive exact number")
    return value


@cache
def data():
    M = spectral.actual(spectral.matrix())
    coefficients = tuple(spectral.actual(value) for value in spectral.coefficients())
    F, Q = principal.finite(), principal.change()
    threshold = clean(F+sum((value/sp.Integer(2*j+1) for j, value in enumerate(coefficients)), sp.zeros(2))/4)
    E_Q = clean(Q.T*principal.asymptotic()*Q)
    infinity = clean(Q[:, 0]*Q[:, 0].T/E_Q[0, 0])
    moment = clean(infinity-F.inv())
    return {"spectral_matrix": M, "local_matrix": F, "chart": Q,
            "threshold_matrix": threshold, "threshold_chart_matrix": clean(Q.T*threshold*Q),
            "threshold_determinant": sp.factor(threshold.det()),
            "infinity_chart_constant": E_Q, "inverse_infinity": infinity,
            "static_spectral_moment": moment, "static_moment_determinant": sp.factor(moment.det())}


def scale(value):
    value = positive_exact(value)
    return sp.ImmutableMatrix([[1/value, 0], [0, 1]])


@cache
def checks():
    item = data()
    M, Q, F, threshold = (item[key] for key in ("spectral_matrix", "chart", "local_matrix", "threshold_matrix"))
    bT, bL = (spectral.actual(spectral.pair(sector)) for sector in ("T", "L"))
    out = {"strict_interior_spectral_determinant": sp.factor(M.det()-512*z**2*(1-z)**2/(6561*h**2)),
           "physical_spectral_Gram": clean(M-2*bT*bT.T-bL*bL.T)}
    denominator = (4*m**2+(1-y**2)*p_real)**2+(1-y**2)**2*p_imag**2
    p = p_real+sp.I*p_imag
    scalar = -p/(4*(4*m**2+(1-y**2)*p))
    out["complex_dispersion_strict_sign_scalar"] = sp.factor(sp.im(scalar)+m**2*p_imag/denominator)
    weight = x*y**2/(1-x+x*y**2)
    out["subthreshold_positive_weight"] = sp.factor(weight-x*y**2/(1-x*(1-y**2)))
    out["subthreshold_weight_below_one"] = sp.factor(1-weight-(1-x)/(1-x+x*y**2))
    independent = clean(F+M.subs(z,y**2).applyfunc(lambda entry: sp.integrate(entry,(y,0,1)))/4)
    out["finite_threshold_by_polynomial_integral"] = clean(threshold-independent)
    out["threshold_negative_scale_diagonal"] = sp.factor(threshold[1,1]+sp.Rational(16,15))
    out["threshold_negative_determinant"] = sp.factor(threshold.det()+sp.Rational(526352,1476225)/h**2)
    out["positive_subthreshold_first_chart_diagonal"] = sp.factor(
        (Q.T*F*Q)[0,0]-sp.Rational(3128,19683)/h**2)
    E = item["infinity_chart_constant"]
    out["positive_high_frequency_complement"] = sp.factor(E[0,0]-sp.Rational(15616,98415)/h**2)
    out["inverse_infinity_explicit"] = clean(item["inverse_infinity"]-sp.Matrix([
        [sp.Rational(98415,15616)*h**2,-sp.Rational(1215,976)*h],
        [-sp.Rational(1215,976)*h,sp.Rational(15,61)]]))
    out["inverse_infinity_has_rank_one"] = sp.factor(item["inverse_infinity"].det())
    out["static_moment_sum_rule_matrix"] = clean(item["static_spectral_moment"]-(item["inverse_infinity"]-F.inv()))
    out["static_moment_strict_determinant"] = sp.factor(
        item["static_spectral_moment"].det()-sp.Rational(177147,99536384)*h**2)
    S = sp.diag(1/h,1)
    for name in ("spectral_matrix","local_matrix","threshold_matrix"):
        out[name+"_exact_congruence"] = clean(item[name]-S*item[name].subs(h,1)*S)
    for name in ("inverse_infinity","static_spectral_moment"):
        out[name+"_inverse_congruence"] = clean(item[name]-S.inv()*item[name].subs(h,1)*S.inv())
    return out
