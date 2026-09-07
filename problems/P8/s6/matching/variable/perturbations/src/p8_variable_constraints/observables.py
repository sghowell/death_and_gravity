"""Actual Bardeen observable jets and conserved scalar symplectic form.

Review-stage algebra. The physical chart is q=(xi,delta_chi_B,Psi_f_B),
where xi=delta_phi_B/phi'=a²(B_g-E_g'), and E_g=delta_phi=0 in the
regular reduction chart. Explicit common time transformation:
xi -> xi+time_shift, delta_phi -> delta_phi-phi'*time_shift,
delta_chi -> delta_chi-chi'*time_shift, q_f -> q_f+h*time_shift.
Thus delta_phi_B=phi'*xi, delta_chi_B=z+chi'*xi and Psi_f_B=q_f-h*xi.

The frozen second-order operator is obtained from q=O Z with BOTH O'
and O'' retained, not by freezing the canonical Hamiltonian. The conserved
bilinear is pulled back from the canonical first-order system; its
q-velocity block is not the original Hamiltonian momentum block.
"""

from functools import cache, lru_cache
from math import comb

import sympy as sp

from . import action, reduction


@cache
def coefficient_rows():
    d, r = action.derive(), reduction.derive()
    a, w, K = d["a"], d["chi_speed"], action.K
    xi = sp.Matrix([[0, 0, 3*a**3*w/(2*a*K), 0, -3/(2*a*K), 0]])
    z = sp.Matrix([[0, 0, 1, 0, 0, 0]])
    R = sp.Matrix([[1, 0, 0, 0, 0, 0]])
    return {"xi": xi, "z": z, "R": R, "alpha": r["alpha"], "h": d["h"]}


@lru_cache(maxsize=None, typed=True)
def observable_jets(time, lapse):
    d, r = action.derive(), coefficient_rows()
    u = d["u"]
    at = {u: sp.sympify(time), d["c"]: sp.sympify(lapse)}
    # At the bounce alpha=0, so the otherwise needed L2 term vanishes.
    at_bounce = sp.sympify(time) == 0
    data = reduction.jets(time, lapse, 1 if at_bounce else 2)
    L = data["p0"]
    result = []
    for j in range(3):
        first = sp.diff(r["xi"], u, j).subs(at)
        second = sp.diff(r["z"]+d["chi_speed"]*r["xi"], u, j).subs(at)
        third = sp.diff(r["R"]-d["h"]*r["xi"], u, j).subs(at)
        for i in range(j+1):
            coefficient = sp.diff(r["alpha"], u, i).subs(at)
            if coefficient != 0:
                third += sp.binomial(j, i)*coefficient*L[j-i]
        result.append(first.col_join(second).col_join(third).applyfunc(sp.cancel))
    return {**data, "O": tuple(result)}


@lru_cache(maxsize=None, typed=True)
def equation(time, lapse):
    data = observable_jets(time, lapse)
    H0, H1 = data["H"][:2]
    O0, O1, O2 = data["O"]
    A0, A1 = reduction.J*H0, reduction.J*H1
    velocity = (O1+O0*A0).applyfunc(sp.cancel)
    phase_map = O0.col_join(velocity)
    determinant = sp.factor(phase_map.det(method="domain-ge"))
    if determinant == 0:
        raise ValueError("The chosen physical scalar Cauchy chart is singular")
    inverse = phase_map.inv().applyfunc(sp.cancel)
    acceleration = (O2+2*O1*A0+O0*(A1+A0*A0)).applyfunc(sp.cancel)
    operator = (acceleration*inverse).applyfunc(sp.cancel)
    symplectic = (inverse.T*reduction.J*inverse).applyfunc(sp.factor)
    return {**data, "A0": A0, "A1": A1, "phase_map": phase_map,
            "det_phase_map": determinant, "inverse_phase_map": inverse,
            "Mq": operator[:, :3], "Mv": operator[:, 3:],
            "symplectic": symplectic,
            "q_q_bracket": (O0*reduction.J*O0.T).applyfunc(sp.cancel),
            "q_velocity_bracket": (O0*reduction.J*velocity.T).applyfunc(sp.factor)}


@lru_cache(maxsize=None, typed=True)
def center(lapse=4):
    """Exact finite-K outputs and a NONUNIFORM frozen-symbol diagnostic.

The returned polynomial is not certified as a physical scalar cone. The
time/momentum controls below explicitly prevent that interpretation.
"""
    d = equation(0, lapse)
    K = action.K
    mq_limit = d["Mq"].applyfunc(lambda item: sp.simplify(sp.limit(item/K, K, sp.oo)))
    mv_limit = d["Mv"].applyfunc(lambda item: sp.simplify(sp.limit(item, K, sp.oo)))
    weight_limit = d["symplectic"][:3, 3:].applyfunc(
        lambda item: sp.simplify(sp.limit(item, K, sp.oo)))
    qq_limit = d["symplectic"][:3, :3].applyfunc(
        lambda item: sp.simplify(sp.limit(item, K, sp.oo)))
    v = sp.Symbol("v", real=True)
    principal = sp.factor((-v**2*sp.eye(3)-mq_limit).det())
    return {**d, "Mq_over_K_limit": mq_limit, "Mv_limit": mv_limit,
            "symplectic_qvelocity_limit": weight_limit, "symplectic_qq_limit": qq_limit,
            "principal": principal, "v": v}


@cache
def center_time_map():
    """Literal Cauchy-map second time jet exhibiting the high-K layer."""
    d, r = action.derive(), coefficient_rows()
    u, K = d["u"], action.K
    at = {u: 0, d["c"]: 4}
    data = reduction.jets(0, 4, 2)
    O = []
    for j in range(4):
        first = sp.diff(r["xi"], u, j).subs(at)
        second = sp.diff(r["z"]+d["chi_speed"]*r["xi"], u, j).subs(at)
        third = sp.diff(r["R"]-d["h"]*r["xi"], u, j).subs(at)
        for i in range(j+1):
            coefficient = sp.diff(r["alpha"], u, i).subs(at)
            if coefficient != 0:
                third += comb(j, i)*coefficient*data["p0"][j-i]
        O.append(first.col_join(second).col_join(third).applyfunc(sp.cancel))
    A = tuple(reduction.J*value for value in data["H"])
    matrices = []
    for j in range(3):
        velocity = O[j+1].copy()
        for i in range(j+1):
            velocity += comb(j, i)*O[i]*A[j-i]
        matrices.append(O[j].col_join(velocity).applyfunc(sp.cancel))
    inverse = matrices[0].inv()
    determinant = sp.factor(matrices[0].det(method="domain-ge"))
    X = inverse*matrices[1]
    first = sp.factor(determinant*sp.trace(X))
    second = sp.factor(determinant*(sp.trace(X)**2-sp.trace(X*X)+sp.trace(inverse*matrices[2])))
    return {"phase_map_jets": tuple(matrices), "O_jets": tuple(O),
            "det0": determinant, "det1": first, "det2": second,
            "u2_K_coefficient": sp.limit(second/(2*K), K, sp.oo),
            "relative_u2_K_coefficient": sp.factor(sp.limit(second/(2*K*determinant), K, sp.oo))}


def nonuniform_controls():
    """Exact failure of the inference 'frozen center polynomial = cone'."""
    u, K = sp.symbols("u K", real=True)
    x = sp.Function("x")(u)
    f = 1+K*u**2
    q = f*x
    equation = sp.diff(q, u, 2)-2*sp.diff(f, u)/f*sp.diff(q, u)
    coefficient = K+2*(sp.diff(f, u)/f)**2-sp.diff(f, u, 2)/f
    equation += coefficient*q
    identity = sp.simplify(equation-f*(sp.diff(x, u, 2)+K*x))
    if identity != 0:
        raise ValueError("The nonuniform derivative-map countercontrol failed")
    return {"exact_transformed_free_oscillator": identity,
            "original_frozen_frequency_squared": K,
            "transformed_center_frozen_frequency_squared": sp.factor(coefficient.subs(u, 0)),
            "map_second_derivative_at_center": sp.diff(f, u, 2).subs(u, 0),
            "meaning": "Same exact solutions under q=(1+K*u²)x; the opposite frozen sign is not an instability proof"}


@cache
def center_checks():
    """Finite-K center algebra; the frozen polynomial is only a control."""
    c, K = action.derive()["c"], action.K
    d = center(c)
    F = c**2*(c-2)*K**2+4*c*(c-2)*(c-8)*K+1536*(8-c)
    bracket = d["q_velocity_bracket"]
    expected_det = 5*F/(12*K**2*(c-2)*(6400-801*c))
    expected_two = 100*F/(K**2*c*(c-2)*(6400-801*c))
    expected_22 = (200*c*(c-2)*K**2-c*(c-2)*K+384)/(200*K**2*c*(c-2))
    result = {
        "physical_Cauchy_determinant": sp.factor(d["det_phase_map"]-expected_det),
        "physical_q_velocity_two_minor": sp.factor(bracket[:2, :2].det()-expected_two),
        "physical_q_velocity_second_diagonal": sp.factor(bracket[1, 1]-expected_22),
        "physical_q_velocity_third_diagonal": sp.factor(bracket[2, 2]-c/240),
        "finite_K_positive_completion": sp.expand(F-c**2*(c-2)*(K+2*(c-8)/c)**2
                                                   -4*(8-c)*(384-(c-2)*(8-c))),
        "second_diagonal_positive_completion": sp.expand(200*c*(c-2)*K**2-c*(c-2)*K+384
                                                            -200*c*(c-2)*(K-sp.Rational(1, 400))**2
                                                            -384+c*(c-2)/800),
    }
    for i, value in enumerate(d["q_q_bracket"]):
        result[f"center_commuting_observables_{i}"] = value
    expected_weight = sp.diag((6400-801*c)/(100*c), 1, 240/c)
    for i, value in enumerate(d["symplectic_qvelocity_limit"]-expected_weight):
        result[f"center_formal_symplectic_weight_{i}"] = sp.factor(value)
    # This is explicitly not a cone proof; replay it to detect omitted map jets.
    sigma = (c-1)*(9*c**2-22*c+144)/120
    result["center_frozen_polynomial_NOT_a_physical_cone"] = sp.factor(
        d["principal"]+(d["v"]**2-1)**2*(d["v"]**2-sigma))
    jet = center_time_map()
    result["actual_map_first_determinant_jet"] = jet["det1"]
    result["actual_map_second_determinant_jet"] = sp.factor(
        jet["det2"]-(19975*K**3-48440*K**2+1644744*K-313661952)/(5745609*K**2))
    result["actual_nonuniform_u2_K_coefficient"] = jet["relative_u2_K_coefficient"]-sp.Rational(5, 6)
    return result
