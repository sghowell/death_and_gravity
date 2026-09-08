"""Uniform frequency-jet bounds with the comoving momentum held fixed."""
from functools import cache

import sympy as sp
from p8_affine_aligned import modes

u, r = modes.u, modes.R


@cache
def identities():
    data = modes.canonical()
    H = data["H"]
    alpha, beta = data["alpha"], data["beta"]
    A1 = sp.factor(r*sp.diff(alpha, r))
    B1 = sp.factor(r*sp.diff(beta, r))
    return {"A1": A1, "B1": B1,
            "H_second": sp.factor(sp.diff(H, u, 2)-8*u*(u**2-3)/(1+u**2)**3),
            "H_third": sp.factor(sp.diff(H, u, 3)+24*(u**4-6*u**2+1)/(1+u**2)**4),
            "alpha_r_weight": sp.factor(A1-r/(1+r)**2),
            "alpha_r_weight_derivative": sp.factor(r*sp.diff(A1, r)-r*(1-r)/(1+r)**3),
            "beta_r_weight": sp.factor(B1-r*(5*r-1)/(1+r)**3),
            "beta_r_weight_derivative": sp.factor(r*sp.diff(B1, r)-r*(-5*r**2+12*r-1)/(1+r)**4),
            "A1_quarter_bound": sp.expand((1+r)**2-4*r-(r-1)**2),
            "B1_absolute_bound": sp.expand((1+r)**3-(5*r**2+r)-(r*(r-1)**2+r+1)),
            "rB1_derivative_absolute_bound": sp.expand(3*(1+r)**4-(5*r**3+12*r**2+r)
                                                       -(3*r**4+7*r**3+6*r**2+11*r+3)),
            "H_third_numerator_bound": sp.expand(3*(1+u**2)**2-(u**4+6*u**2+1)-2*(u**4+1))}


@cache
def time_jets():
    H, Hp, Hpp, Hppp = sp.symbols("H H_prime H_second H_third", real=True)
    q = sp.Symbol("q", nonnegative=True)
    alpha, beta = modes.canonical()["alpha"], modes.canonical()["beta"]

    def derivative(value):
        return sp.expand(sp.diff(value, H)*Hp+sp.diff(value, Hp)*Hpp
                         +sp.diff(value, Hpp)*Hppp-2*H*r*sp.diff(value, r)-2*H*q*sp.diff(value, q))

    U = alpha*Hp+beta*H**2
    Up = derivative(U)
    Upp = derivative(Up)
    ap, app, bp, bpp = derivative(alpha), derivative(derivative(alpha)), derivative(beta), derivative(derivative(beta))
    target_Up = ap*Hp+alpha*Hpp+bp*H**2+2*beta*H*Hp
    target_Upp = app*Hp+2*ap*Hpp+alpha*Hppp+bpp*H**2+4*bp*H*Hp+2*beta*(Hp**2+H*Hpp)
    return {"H": H, "Hp": Hp, "Hpp": Hpp, "Hppp": Hppp, "q": q,
            "U": U, "U_prime": Up, "U_second": Upp,
            "alpha_prime": ap, "alpha_second": app, "beta_prime": bp, "beta_second": bpp,
            "longitudinal_curvature_first_jet": sp.factor(Up-target_Up),
            "longitudinal_curvature_second_jet": sp.factor(Upp-target_Upp),
            "comoving_q_second_jet": sp.factor(derivative(derivative(q))-(4*H**2-2*Hp)*q)}


@cache
def proof_checks():
    H, Hp, Hpp, Hppp = 2, 4, 12, 72
    ap = 2*H*sp.Rational(1, 4)
    app = 2*Hp*sp.Rational(1, 4)+4*H**2*sp.Rational(1, 4)
    bp, bpp = 2*H, 2*Hp+4*H**2*3
    Up = ap*Hp+sp.Rational(3, 2)*Hpp+bp*H**2+2*sp.Rational(9, 4)*H*Hp
    Upp = (app*Hp+2*ap*Hpp+sp.Rational(3, 2)*Hppp+bpp*H**2
           +4*bp*H*Hp+2*sp.Rational(9, 4)*(Hp**2+H*Hpp))
    return {"H_second_absolute_bound": 8*sp.Rational(1, 2)*3 == Hpp,
            "H_third_absolute_bound": 24*3 == Hppp,
            "alpha_first_time_bound": ap == 1,
            "alpha_second_time_bound": app == 6,
            "beta_first_time_bound": bp == 4,
            "beta_second_time_bound": bpp == 56,
            "longitudinal_U_first_bound": Up == 74,
            "longitudinal_U_second_bound": Upp == 688,
            "frequency_first_bound_for_q_le_one": 2*H+Up == 78,
            "frequency_second_bound_for_q_le_one": 4*H**2+2*Hp+Upp == 712,
            "normalization_log_rate_bound": H*sp.Rational(3, 2) == 3,
            "normalization_third_ratio_bound": Up+15*3 == 119}


@cache
def checks():
    data = identities()
    out = {name: value for name, value in data.items() if name not in ("A1", "B1")}
    jets = time_jets()
    for name in ("longitudinal_curvature_first_jet", "longitudinal_curvature_second_jet", "comoving_q_second_jet"):
        out[name] = jets[name]
    bg = modes.canonical()
    substitute = {jets["H"]: bg["H"], jets["Hp"]: sp.diff(bg["H"], u),
                  jets["Hpp"]: sp.diff(bg["H"], u, 2), jets["Hppp"]: sp.diff(bg["H"], u, 3), r: bg["r"]}
    out["actual_longitudinal_first_derivative"] = sp.factor(
        sp.diff(bg["longitudinal"]["correction"], u)-jets["U_prime"].subs(substitute))
    out["actual_longitudinal_second_derivative"] = sp.factor(
        sp.diff(bg["longitudinal"]["correction"], u, 2)-jets["U_second"].subs(substitute))
    return out
