"""Exact fourth-order subtraction block, not the full coupled propagator."""
from functools import cache

import sympy as sp
from p8_vector_metric_local import consistency, counterterms, jets

from . import spectral

h, z = spectral.h, spectral.z


def finite():
    return sp.ImmutableMatrix([[sp.Rational(1256, 6561)/h**2, -sp.Rational(76, 243)/h],
                               [-sp.Rational(76, 243)/h, -4]])


def pole():
    v = sp.ImmutableMatrix([sp.Rational(16, 81)/h, 1])
    return 2*v*v.T


def change():
    return sp.ImmutableMatrix([[1, 0], [-sp.Rational(16, 81)/h, 1]])


def asymptotic():
    return sp.ImmutableMatrix((finite()+spectral.actual(spectral.finite_asymptotic())).applyfunc(sp.factor))


@cache
def real_axis():
    Q = change()
    transformed = sp.ImmutableMatrix((Q.T*spectral.actual(spectral.matrix())*Q).applyfunc(sp.factor))
    F = sp.ImmutableMatrix((Q.T*finite()*Q).applyfunc(sp.factor))
    E = sp.ImmutableMatrix((Q.T*asymptotic()*Q).applyfunc(sp.factor))
    lower = E[0, 0]
    upper_cross = E[0, 1]
    inverse = sp.ImmutableMatrix([[1/lower, upper_cross/(4*lower)],
                                  [upper_cross/(4*lower), sp.Rational(1, 4)]])
    return {"transformed_spectral_matrix": transformed,
            "transformed_local_matrix": F, "transformed_asymptotic_constant": E,
            "positive_diagonal_lower": lower, "negative_diagonal_upper": -sp.Integer(4),
            "positive_cross_lower": F[0, 1], "positive_cross_upper": upper_cross,
            "determinant_upper": -4*lower,
            "absolute_inverse_entry_upper": inverse,
            "uniform_interval_inverse_row_sum_upper": max(
                sum(inverse.row(j)).subs(h, sp.Rational(125, 64)) for j in range(2)),
            "scope": "Real s>0, frozen fourth-order local term plus exact thrice-subtracted flat bubble, divided by s^4/(64*pi^2). This is neither the full tree-plus-loop symbol nor a causal operator norm."}


@cache
def local_checks():
    out = {}
    actual_h = (1+jets.u**2)**3
    for output_index, output in enumerate(("N", "Z")):
        for source_index, source in enumerate(("N", "Z")):
            row = jets.n if source == "N" else jets.v
            actual_pole = sp.expand(counterterms.operator(output, 2)).coeff(row[4])
            actual_finite = consistency.coefficients(output, 2, source)[4]
            out[output+source+"_frozen_curved_pole_principal"] = sp.factor(
                actual_pole-pole()[output_index, source_index].subs(h, actual_h))
            out[output+source+"_frozen_curved_finite_principal"] = sp.factor(
                actual_finite-finite()[output_index, source_index].subs(h, actual_h))
    return out


@cache
def checks():
    out = {}
    Q = change()
    data = real_axis()
    M = data["transformed_spectral_matrix"]
    out["triangular_spectral_chart_has_unit_determinant"] = sp.factor(Q.det()-1)
    out["one_uv_pole_channel"] = sp.factor(pole().det())
    out["finite_matrix_is_not_rank_one"] = sp.factor(finite().det()+sp.Rational(50992, 59049)/h**2)
    out["real_axis_A_spectral_positive_square"] = sp.factor(M[0, 0]-16*(1-z)**2/(2187*h**2))
    out["real_axis_B_spectral_negative"] = sp.factor(M[0, 1]+8*(3-z)*(1-z)/(81*h))
    out["real_axis_C_spectral_strictly_positive"] = sp.factor(
        M[1, 1]-12*(z-sp.Rational(1, 3))**2-sp.Rational(32, 3))
    y = spectral.y
    for i, j, name in ((0, 0, "A"), (0, 1, "B")):
        integrand = sp.cancel(y**2*M[i, j].subs(z, y**2)/(1-y**2))/4
        integral = sp.integrate(integrand, (y, 0, 1))
        out["finite_limit_"+name+"_by_independent_polynomial_integral"] = sp.factor(
            data["transformed_local_matrix"][i, j]-integral
            -data["transformed_asymptotic_constant"][i, j])
    ell = sp.Symbol("positive_axis_log_s_over_m", real=True)
    leading = asymptotic()-2*pole()*ell
    out["leading_log_determinant"] = sp.factor(leading.det()
        +16*(58560*ell+9137)/(1476225*h**2))
    out["spurious_asymptotic_root_log"] = sp.factor(
        leading.det().subs(ell, -sp.Rational(9137, 58560)))
    out["positive_real_exact_determinant_gap"] = sp.factor(
        data["determinant_upper"]+sp.Rational(62464, 98415)/h**2)
    return out
