"""Exact anisotropic tensor restriction of the formal own-f S_(6) term.

This is an off-shell unit-volume physical-g probe, not the full CD tensor
background and not a DHOST coefficient extraction.  All time derivatives
are taken before the center restriction.  M**2 is factored from densities.
"""

from functools import cache

import sympy as sp

M, TAU, C = sp.symbols("M tau c", positive=True)
R2 = sp.Symbol("r_squared", positive=True)
INV = sp.symbols("inv0:10", real=True)
LOG = sp.symbols("logr0:10", real=True)
Q = sp.symbols("q0:10", real=True)
EPS = sp.Symbol("epsilon", real=True)
AMPLITUDE = sp.Symbol("amplitude", real=True)


def total_derivative(expression):
    """Time jet derivation: inv=A=M²r/beta1, logr=(log r)', q_j=q^(j)."""
    expression = sp.sympify(expression)
    if any(symbol in expression.free_symbols for symbol in (INV[-1], LOG[-1], Q[-1])):
        raise ValueError("The finite jet array does not contain the required next derivative")
    result = 2*LOG[0]*R2*sp.diff(expression, R2)
    for jets in (INV, LOG, Q):
        result += sum(jets[j+1]*sp.diff(expression, jets[j]) for j in range(len(jets)-1))
    return sp.expand(result)


def _quadratic(expression):
    return sp.expand(expression.subs({q: AMPLITUDE*q for q in Q})).coeff(AMPLITUDE, 2)


def _euler(expression, highest=3):
    answer = sp.diff(expression, Q[0])
    for n in range(1, highest+1):
        term = sp.diff(expression, Q[n])
        for _ in range(n):
            term = total_derivative(term)
        answer += (-1)**n*term
    return sp.expand(answer)


@cache
def derive():
    """Literal diagonal coframes, formal hidden-metric path, and boundaries.

    r(t), A(t), and hence an arbitrary smooth composed clock theta(t), are
    retained.  No H=0 or clock-center substitution occurs here.
    """
    a, ell = INV[0], LOG[0]
    qp, qpp = Q[1], Q[2]
    P = (-5*qp**2/sp.Integer(12)-2*LOG[1]+ell**2,
         -qpp/2+qp**2/12-ell*qp-ell**2,
         qpp/2+qp**2/12+ell*qp-ell**2,
         qp**2/12-ell**2)
    h = tuple(a*p for p in P)
    dh = tuple(total_derivative(component) for component in h)
    rates = (ell+qp/2, ell-qp/2, ell)
    V1 = (sum(h[1:])-h[0])/2
    V2 = (sum(h[1:])-h[0])**2/8-(sum(component**2 for component in h[1:])-h[0]**2)/4
    S0 = 3*ell**2-qp**2/4
    S1 = sum(dh[i+1]*(3*ell-rates[i]) for i in range(3))/2
    S2 = (sum(dh[i+1]*dh[j+1] for i in range(3) for j in range(i+1, 3))/4
          -sum(h[i+1]*dh[i+1]*(3*ell-rates[i]) for i in range(3))/2)
    EH = -R2*(S2+V1*S1+V2*S0)
    # b=r*beta1=M²*r²/A.  The epsilon³ potential coefficient is S0'''/6.
    potential = R2*(sum(h)**3-6*sum(h)*sum(component**2 for component in h)
                    +5*sum(component**3 for component in h))/(24*a)
    raw = _quadratic(EH+potential)
    coefficients = {(i, j): sp.factor(sp.diff(raw, Q[i], Q[j])/(2 if i == j else 1))
                    for i in range(1, 4) for j in range(i, 4)}
    B3 = coefficients[3, 3]
    B2 = sp.factor(coefficients[2, 2]-total_derivative(coefficients[2, 3])/2-coefficients[1, 3])
    B1 = sp.factor(coefficients[1, 1]+total_derivative(total_derivative(coefficients[1, 3]))/2
                   -total_derivative(coefficients[1, 2])/2)
    normal = B3*Q[3]**2+B2*Q[2]**2+B1*Q[1]**2
    boundary = (coefficients[2, 3]*Q[2]**2/2+coefficients[1, 3]*Q[1]*Q[2]
                +(coefficients[1, 2]-total_derivative(coefficients[1, 3]))*Q[1]**2/2)
    # Covariant EH differs from ADM by +M² d_t[V_f/N_f sum beta_i'].
    EH_boundary = R2*(3*ell*V2+V1*sum(dh[1:])/2-sum(h[i]*dh[i] for i in range(1, 4))/2)
    EH_boundary_quadratic = _quadratic(EH_boundary)
    euler = _euler(normal)
    euler_coefficients = {n: sp.factor(euler.coeff(Q[n])) for n in range(7)}
    kappa_over_M2 = R2*a/4
    fourth_raw = kappa_over_M2*_quadratic(sum(p**2 for p in P)-sum(P)**2)
    fourth_normal = kappa_over_M2*Q[2]**2/2-total_derivative(kappa_over_M2)*ell*Q[1]**2
    fourth_boundary = kappa_over_M2*ell*Q[1]**2
    fourth_euler = _euler(fourth_normal, highest=2)
    return {"P_mixed_diagonal": P, "h_mixed_diagonal": h, "h_time_derivative": dh,
            "background_log_rates": rates, "volume_order1": V1, "volume_order2": V2,
            "rate_pair_order0": S0, "rate_pair_order1": S1, "rate_pair_order2": S2,
            "EH_second_coefficient": EH, "potential_third_coefficient": potential,
            "raw_quadratic": raw, "raw_coefficients": coefficients,
            "B1": B1, "B2": B2, "B3": B3, "normal_quadratic": normal,
            "time_IBP_boundary": boundary, "covariant_EH_boundary_second": EH_boundary,
            "covariant_EH_boundary_second_quadratic": EH_boundary_quadratic,
            "covariant_to_normal_boundary": EH_boundary_quadratic+boundary,
            "euler": euler, "euler_coefficients": euler_coefficients,
            "kappa_over_M2": kappa_over_M2, "fourth_raw": fourth_raw,
            "fourth_normal": fourth_normal, "fourth_boundary": fourth_boundary,
            "fourth_euler": fourth_euler,
            "fourth_euler_coefficients": {n: sp.factor(fourth_euler.coeff(Q[n])) for n in range(5)},
            "metric_path": sp.diag(R2*(1+EPS*h[0]), -R2*sp.exp(Q[0])*(1+EPS*h[1]),
                                    -R2*sp.exp(-Q[0])*(1+EPS*h[2]), -R2*(1+EPS*h[3]))}


@cache
def center_jets():
    """Independent rational-profile derivation; differentiation precedes u=0.

    These substitutions use u-time and A_bar=A/tau², ell_bar=tau*ell.
    They pull back the unchanged canonical clock via theta=T=tau*u.
    """
    v, c = sp.Symbol("v", real=True), C
    d = 1+v
    y = 2/d**4
    hp, h_squared = 4*(1-v)/d**2, 16*v/d**2
    b1 = y**3*hp/(c*(c-y))
    b4 = 3*h_squared/(2*c**2)-b1/y**3
    r3 = sp.factor(-b1/b4)
    kappa_bar = sp.factor(-1/(4*b4))
    ell_over_u = sp.factor(2*sp.diff(r3, v)/(3*r3))
    r3v, r3vv = (sp.factor(sp.diff(r3, v, n).subs(v, 0)) for n in (1, 2))
    r2v, r2vv = r3v/3, r3vv/3-r3v**2/72
    k = tuple(sp.factor(sp.diff(kappa_bar, v, n).subs(v, 0)) for n in range(3))
    # A_bar=4*kappa_bar/r²; the algebraic root has r²(0)=4.
    A0 = k[0]
    Av = k[1]-k[0]*r2v/4
    Avv = k[2]-k[1]*r2v/2-k[0]*r2vv/4+k[0]*r2v**2/8
    point = {R2: sp.Integer(4), INV[0]: A0, INV[1]: sp.S.Zero,
             INV[2]: sp.factor(2*Av), INV[3]: sp.S.Zero, INV[4]: sp.factor(12*Avv),
             INV[5]: sp.S.Zero, LOG[0]: sp.S.Zero, LOG[1]: sp.factor(ell_over_u.subs(v, 0)),
             LOG[2]: sp.S.Zero, LOG[3]: sp.factor(6*sp.diff(ell_over_u, v).subs(v, 0)),
             LOG[4]: sp.S.Zero, LOG[5]: sp.factor(60*sp.diff(ell_over_u, v, 2).subs(v, 0))}
    return {"v": v, "b1": b1, "b4": b4, "r_cubed": r3,
            "kappa_bar": kappa_bar, "log_r_u_over_u": ell_over_u,
            "substitution": point, "A_bar_u0_u2_u4": (point[INV[0]], point[INV[2]], point[INV[4]]),
            "log_r_u1_u3": (point[LOG[1]], point[LOG[3]])}


@cache
def center():
    d, j = derive(), center_jets()
    point = j["substitution"]
    B = {n: sp.factor(d[f"B{n}"].subs(point)) for n in (1, 2, 3)}
    E = {n: sp.factor(value.subs(point)) for n, value in d["euler_coefficients"].items()}
    E4 = {n: sp.factor(value.subs(point)) for n, value in d["fourth_euler_coefficients"].items()}
    combined = {n: sp.factor(E[n]+E4.get(n, 0)) for n in range(7)}
    delta = sp.Symbol("delta", positive=True)
    low_delta = sp.Poly(sp.expand(E[2].subs(C, 2+delta)), delta)
    combined_delta = sp.Poly(sp.expand(combined[2].subs(C, 2+delta)), delta)
    return {"B_dimensionless": B, "sixth_E_dimensionless": E, "fourth_E_dimensionless": E4,
            "combined_E_dimensionless": combined,
            "B_physical_including_M2": {n: M**2*TAU**(2*n-2)*value for n, value in B.items()},
            "sixth_E_physical": {n: M**2*TAU**(n-2)*value for n, value in E.items()},
            "fourth_E_physical": {n: M**2*TAU**(n-2)*value for n, value in E4.items()},
            "sixth_low_delta_polynomial": low_delta, "combined_low_delta_polynomial": combined_delta,
            "small_delta_domain": "0<delta<=1/100, delta=c-2",
            "sixth_low_absolute_lower_bound": sp.Rational(793, 400),
            "domain": "M,tau>0; 2<c<=4; coefficients from the frozen |u|<=1/10 family; theta=T at evaluation"}


@cache
def checks():
    d, j, c = derive(), center_jets(), center()
    A, Ap, App, Appp = INV[:4]
    ell, ellp, ellpp = LOG[:3]
    expected_B3 = R2*A**2/16
    expected_B2 = -R2*A*(7*A*ell**2+2*A*ellp+6*Ap*ell+App)/16
    expected_B1 = R2*(23*A**2*ell**4+4*A**2*ell**2*ellp+8*A**2*ell*ellpp+8*A**2*ellp**2
                     +24*A*Ap*ell**3+24*A*Ap*ell*ellp+2*A*Ap*ellpp+6*A*App*ell**2
                     +4*A*App*ellp+2*A*Appp*ell+8*Ap**2*ell**2+2*Ap**2*ellp+2*Ap*App*ell)/16
    raw_euler = _euler(d["raw_quadratic"])
    point = j["substitution"]
    values = {"exact_quadratic_time_boundary": d["raw_quadratic"]-d["normal_quadratic"]-total_derivative(d["time_IBP_boundary"]),
              "boundary_Euler_zero": _euler(total_derivative(d["time_IBP_boundary"])),
              "raw_and_normal_Euler": raw_euler-d["euler"],
              "B3_formula": d["B3"]-expected_B3, "B2_formula": d["B2"]-expected_B2,
              "B1_formula": d["B1"]-expected_B1,
              "constant_anisotropy_zero_symbol": d["euler_coefficients"][0],
              "highest_sixth_symbol": d["euler_coefficients"][6]+R2*A**2/8,
              "fourth_stationary_square": d["fourth_raw"]-d["kappa_over_M2"]*(Q[2]**2/2+2*ell*Q[1]*Q[2]+ellp*Q[1]**2),
              "fourth_exact_time_boundary": d["fourth_raw"]-d["fourth_normal"]-total_derivative(d["fourth_boundary"]),
              "fourth_raw_and_normal_Euler": _euler(d["fourth_raw"], highest=2)-d["fourth_euler"],
              "fourth_lower_symbol_before_center": d["fourth_euler_coefficients"][2]
                  -total_derivative(total_derivative(d["kappa_over_M2"]))-2*ell*total_derivative(d["kappa_over_M2"]),
              "inverse_coefficient_center": point[INV[0]]-C*(C-2)/16,
              "inverse_second_jet_center": point[INV[2]]-(13*C**2-22*C+8)/8,
              "log_ratio_first_jet_center": point[LOG[1]]+4*(C+2)/C}
    expected_B = {3: C**2*(C-2)**2/1024,
                  2: -C*(C-2)*(9*C**2-22*C+24)/512,
                  1: -(C-2)*(C+2)*(9*C**2-22*C+24)/32}
    expected_E = {6: -C**2*(C-2)**2/512,
                  4: -C*(C-2)*(75*C**2-154*C+120)/256,
                  2: -(349*C**4-2220*C**3+5780*C**2-7872*C+5056)/128}
    values.update({f"center_B{n}": c["B_dimensionless"][n]-value for n, value in expected_B.items()})
    values.update({f"center_E{n}": c["sixth_E_dimensionless"][n]-value for n, value in expected_E.items()})
    values.update({f"center_odd_E{n}": c["sixth_E_dimensionless"][n] for n in (1, 3, 5)})
    values.update({"center_fourth_lower_symbol": c["fourth_E_dimensionless"][2]-(9*C**2-22*C+24)/8,
                   "sixth_nonzero_low_limit": c["sixth_E_dimensionless"][2].subs(C, 2)+2,
                   "fourth_nonzero_low_limit": c["fourth_E_dimensionless"][2].subs(C, 2)-2,
                   "combined_low_limit_cancels": c["combined_E_dimensionless"][2].subs(C, 2),
                   "nonzero_low_limit_from_differentiated_B2": 2*total_derivative(total_derivative(d["B2"])).subs(point).subs(C, 2)+2,
                   "B1_itself_tends_to_zero": c["B_dimensionless"][1].subs(C, 2)})
    values = {name: sp.factor(sp.expand(value)) for name, value in values.items()}
    if any(value != 0 for value in values.values()):
        raise ValueError("An anisotropic sixth-order identity failed")
    return values


def calibration():
    c = center()
    if c["sixth_low_delta_polynomial"].all_coeffs() != [
            -sp.Rational(349, 128), -sp.Rational(143, 32), -sp.Rational(209, 32), sp.Rational(7, 4), -2]:
        raise ValueError("The continuous low-symbol polynomial failed")
    if 2-sp.Rational(7, 4)/100 != c["sixth_low_absolute_lower_bound"]:
        raise ValueError("The small-delta lower bound failed")
    return {"claim": "P8-S6.32.SIXTH", "sixth_highest_symbol": -M**2*R2*INV[0]**2/8,
            "sixth_center_low_limit": -2*M**2, "fourth_center_low_limit": 2*M**2,
            "combined_center_low_limit": 0,
            "sixth_low_over_Mstar_squared_limit": -sp.Rational(2, 5),
            "sixth_small_delta_low_absolute_floor_over_M2": c["sixth_low_absolute_lower_bound"],
            "scope": "off-shell unit-volume physical-g Euler contribution; not Xi, full CD propagation, an all-order estimate or UV exclusion",
            "scale_dictionary": "density B_n=M² tau^(2n-2)*Bbar_n; Euler q^(n) coefficient=M² tau^(n-2)*Ebar_n"}
