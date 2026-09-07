"""Actual canonical matter vierbein and its first-order constraint dictionary.

epsilon is the source paper's matter-strength expansion parameter. It is
not the earlier composite family parameter 1/Y. The physical metric of
the full parent is h=u.T*eta*u; h0 is only its leading auxiliary solution.
"""

from functools import cache

import sympy as sp

from . import model


def canonical_source(u, gradient, potential=0):
    u, k = model.matrix(u), sp.Matrix(gradient)
    if k.shape != (4, 1):
        raise ValueError("A four-component covariant scalar gradient is required")
    h = u.T*model.ETA*u
    raised = h.inv()*k
    kinetic = (k.T*raised)[0]
    ju = u.det()*((kinetic/2-potential)*u.inv().T-model.ETA*u*raised*raised.T)
    source = model.trace_reduced_source(u, ju)
    return {"h": h, "Y": kinetic, "L_m": u.det()*(kinetic/2-potential),
            "J_u": ju, **source}


def first_order(e, v, gradient, potential=0, *, b=model.B, pg=model.PG, pf=model.PF):
    """Exact formal epsilon coefficients; require b!=0 and invertible u0.

    No finite-energy remainder or uniform convergence is inferred here.
    """
    b = sp.sympify(b)
    if b == 0:
        raise ValueError("The formal auxiliary expansion requires B!=0")
    e, v = model.matrix(e), model.matrix(v)
    u0 = -3*(pg*e+pf*v)/b
    source = canonical_source(u0, gradient, potential)
    u1 = 3*source["calT"]/b
    h1 = u0.T*model.ETA*u1+u1.T*model.ETA*u0
    return {"u0": u0, "u1": u1, "h0": source["h"], "h1": h1, **source}


@cache
def checks():
    us = sp.symbols("u0:4", positive=True)
    ks = sp.symbols("k0:4", real=True)
    potential = sp.Symbol("V", real=True)
    u, k = sp.diag(*us), sp.Matrix(ks)
    out = canonical_source(u, k, potential)
    result = {f"literal_canonical_u_Euler_{i}": sp.cancel(sp.diff(out["L_m"], us[i])-out["J_u"][i, i]) for i in range(4)}
    expected = u*((-out["Y"]/12-potential/6)*sp.eye(4)+out["h"].inv()*k*k.T/2)
    result["canonical_trace_reversal_00"] = sp.cancel(out["calT"][0, 0]-expected[0, 0])
    result["canonical_trace_reversal_01"] = sp.cancel(out["calT"][0, 1]-expected[0, 1])
    result["canonical_source_trace"] = sp.cancel(out["trace_T"]-out["Y"]/2+2*potential)
    u1 = 3*out["calT"]/model.B
    h1 = u.T*model.ETA*u1+u1.T*model.ETA*u
    predicted = (3*k*k.T-(out["Y"]/2+potential)*out["h"])/model.B
    result["actual_matter_metric_correction_00"] = sp.cancel(h1[0, 0]-predicted[0, 0])
    result["actual_matter_metric_correction_01"] = sp.cancel(h1[0, 1]-predicted[0, 1])
    result["actual_matter_metric_correction_11"] = sp.cancel(h1[1, 1]-predicted[1, 1])
    # A source with only a constant potential merely shifts B in the
    # FULL uneliminated action, not only in a truncated effective action.
    zero_source = canonical_source(u, sp.zeros(4, 1), potential)
    e, v = sp.eye(4), sp.diag(1, 2, 3, 4)
    with_source = model.euler_maps(e, v, u, matter_gradient=zero_source["J_u"])
    shifted = model.euler_maps(e, v, u, b=model.B+model.EPS*potential/2, epsilon=0)
    result["constant_potential_full_B_shift"] = sp.cancel(with_source["E_u"][0, 0]-shifted["E_u"][0, 0])
    return result


def controls():
    # Leading e-v symmetry holds, but retaining it at order epsilon would
    # violate the sourced Lorentz constraint. This is a local off-shell
    # dictionary fixture, not a complete cosmological solution.
    e, v = sp.diag(1, 2, 3, 4), sp.eye(4)
    k = sp.Matrix([1, 1, 0, 0])
    out = first_order(e, v, k, b=-6, pg=1, pf=1)
    anti = e.T*model.ETA*out["calT"]-out["calT"].T*model.ETA*e
    return {"leading_e_v_antisymmetric_01": (e.T*model.ETA*v-v.T*model.ETA*e)[0, 1],
            "positive_canonical_timelike_Y": sp.cancel(out["Y"]),
            "required_e_v_antisymmetric_epsilon_coefficient_01": sp.cancel(anti[0, 1]),
            "actual_h_correction_01": sp.cancel(out["h1"][0, 1]),
            "constant_V_h_correction_over_V_h0": -1/model.B}
