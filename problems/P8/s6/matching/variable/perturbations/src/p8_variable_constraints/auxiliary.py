"""Recover both lapses, both shifts and the actual metric perturbations.

The f lapse follows the secondary-momentum preservation equation, not a
division by H. These rows act on the six regular canonical phase entries.
All background coefficients are differentiated before time evaluation.
"""

from functools import lru_cache

import sympy as sp

from . import action, observables, reduction


@lru_cache(maxsize=None, typed=True)
def derive(time, lapse):
    d, r = action.derive(), reduction.derive()
    u = d["u"]
    at = {u: sp.sympify(time), d["c"]: sp.sympify(lapse)}
    data = reduction.jets(time, lapse, 1)
    L0, L1 = data["p0"]
    A = reduction.J*data["H"][0]
    T = r["old_state"].subs(reduction.Q, 0).jacobian(reduction.EXTENDED)
    T0 = T.subs(at).applyfunc(sp.cancel)
    T1 = sp.diff(T, u).subs(at).applyfunc(sp.cancel)
    elimination = sp.eye(6).col_join(L0)
    elimination_prime = sp.zeros(6).col_join(L1)
    old = (T0*elimination).applyfunc(sp.cancel)
    old_prime = (T1*elimination+T0*elimination_prime+old*A).applyfunc(sp.cancel)
    zero = dict.fromkeys(d["state"], 0)

    def linear_row(expression):
        row = sp.Matrix([[sp.diff(expression, value).subs(zero) for value in d["state"]]])
        return (row.subs(at)*old).applyfunc(sp.cancel)

    # The Q derivative of the transformed Hamiltonian is the old q_g
    # derivative; its time canonical correction is independent of Q.
    q_force = linear_row(sp.diff(d["H"], action.QG))
    kappa = (3*d["c"]*d["a"]**3*d["P"]).subs(at)
    nf = (-(q_force+L1+L0*A)/kappa).applyfunc(sp.cancel)
    ng = linear_row(d["S_g"]-(d["clock_constraint_numerator"]-d["chi_speed"]*action.PX)
                    /(d["a"]**3*d["kbar"]))
    bg = linear_row(d["B_g"])
    bf = (bg-old[6:7, :]/(2*d["C_shift"].subs(at))).applyfunc(sp.cancel)
    a, b, c, h = (d[key].subs(at) for key in ("a", "b", "c", "h"))
    xi = a**2*bg
    oxi = observables.coefficient_rows()["xi"]
    xi_prime = (sp.diff(oxi, u).subs(at)+xi*A).applyfunc(sp.cancel)
    return {**data, "old_phase": old, "old_phase_prime": old_prime,
            "n_g": ng, "n_f": nf, "B_g": bg, "B_f": bf,
            "xi": xi, "xi_prime": xi_prime,
            "Psi_g_B": (old[0:1, :]+h*xi).applyfunc(sp.cancel),
            "Psi_f_B": (old[1:2, :]-h*xi).applyfunc(sp.cancel),
            "n_g_B": (ng+xi_prime).applyfunc(sp.cancel),
            "n_f_B": (nf+xi_prime).applyfunc(sp.cancel),
            "B_f_B": (bf-c**2*xi/b**2).applyfunc(sp.cancel)}


def checks(time, lapse):
    """Replay every lapse/shift, momentum definition and phase equation."""
    d, lag, data = action.derive(), action.lagrangian(), derive(time, lapse)
    at = {d["u"]: sp.sympify(time), d["c"]: sp.sympify(lapse)}
    variables = (*d["coordinates"], *lag["velocities"], *lag["auxiliaries"])
    rows = data["old_phase"][:4, :].col_join(data["old_phase_prime"][:4, :])
    rows = rows.col_join(data["n_g"]).col_join(data["n_f"]).col_join(data["B_g"]).col_join(data["B_f"])
    residuals = {}
    for name, equation in zip(("g_lapse", "f_lapse", "g_shift", "f_shift"), lag["auxiliary_equations"], strict=True):
        coefficient = sp.Matrix([[sp.diff(equation, variable) for variable in variables]])
        row = (coefficient.subs(at)*rows).applyfunc(sp.cancel)
        for i, value in enumerate(row):
            residuals[f"reconstructed_{name}_{i}"] = value
    for j, velocity in enumerate(lag["velocities"]):
        equation = sp.diff(lag["L"], velocity)
        coefficient = sp.Matrix([[sp.diff(equation, variable) for variable in variables]])
        row = (coefficient.subs(at)*rows-data["old_phase"][j+4:j+5, :]).applyfunc(sp.cancel)
        for i, value in enumerate(row):
            residuals[f"reconstructed_momentum_{j}_{i}"] = value
    jold = sp.zeros(8)
    jold[:4, 4:] = sp.eye(4)
    jold[4:, :4] = -sp.eye(4)
    hessian = sp.hessian(d["H"], d["state"]).subs(at)
    constraint = sp.Matrix([sp.diff(d["C_f"], value) for value in d["state"]]).subs(at)
    phase_equation = (data["old_phase_prime"]-jold*(hessian*data["old_phase"]
                                                   +constraint*data["n_f"])).applyfunc(sp.cancel)
    for i, value in enumerate(phase_equation):
        residuals[f"full_unreduced_phase_equation_{i}"] = value
    return residuals
