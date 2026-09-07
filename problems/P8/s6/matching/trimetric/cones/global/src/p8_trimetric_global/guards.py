"""Full-solution controls versus merely algebraic/spectral diagnostics."""

from functools import cache

import sympy as sp
from p8_trimetric import model


@cache
def checks():
    ident = sp.eye(4)
    # A genuine flat point with mixed nonzero links but singular TT inverse.
    # It is covered by the primary BACKGROUND theorem, not by the separate
    # regular positive-FP inverse corollary.
    singular = model.euler_maps(ident, ident, ident, b=0, pg=1, pf=-1, bg=-1, bf=1, epsilon=0)
    result = {f"mixed_singular_full_flat_{key}": singular[key][0, 0] for key in ("E_e", "E_v", "E_u")}
    pe, pv, ge, gf, gu = sp.symbols("P_g P_f gamma_g gamma_f gamma_u", real=True)
    lag = -(pe*(ge-gu)**2+pv*(gf-gu)**2)/4
    result["general_TT_auxiliary_Hessian"] = sp.expand(sp.diff(lag, gu, 2)+(pe+pv)/2)
    result["singular_TT_is_a_constraint_not_an_inverse"] = sp.expand(sp.diff(lag.subs(pv, -pe), gu)-pe*(ge-gf)/2)
    # All-zero links: physical u has no background or tensor equation when
    # the canonical scalar is constant and B=V=all endpoints=0. A completely
    # arbitrary metric is not an invertible auxiliary matching construction.
    t = sp.Symbol("t", real=True)
    au = 1+t**2
    u = au*ident
    disconnected = model.euler_maps(ident, ident, u, b=0, pg=0, pf=0, bg=0, bf=0, epsilon=0)
    for key in ("E_e", "E_v", "E_u"):
        result[f"all_zero_links_full_arbitrary_u_{key}"] = disconnected[key][0, 0]
    hu = sp.diff(au, t)/au**2
    proper_time = t+t**3/3
    result["all_zero_arbitrary_u_proper_clock"] = sp.diff(proper_time, t)-au
    result["all_zero_arbitrary_u_positive_bounce"] = sp.diff(hu, t).subs(t, 0)-2
    # Nonzero links and NEC=0 need not be static: actual common de Sitter.
    # H_u=1/2, G=F=q=1,R_i=c_i=1,B=-6,b_i=-5/8,V=0.
    einstein = sp.diag(sp.Rational(3, 4), *([-sp.Rational(3, 4)]*3))
    de_sitter = model.euler_maps(ident, ident, ident, b=-6, pg=1, pf=1,
                                bg=-sp.Rational(5, 8), bf=-sp.Rational(5, 8), g=1, f=1,
                                einstein_g=einstein, einstein_f=einstein, epsilon=0)
    for key in ("E_e", "E_v", "E_u"):
        result[f"actual_deSitter_nonzero_link_{key}"] = de_sitter[key][0, 0]
    return result


def controls():
    ident = sp.eye(4)
    # These calibrated e/v equations do not establish a full vacuum: B is
    # deliberately incompatible and the (free, V=0) u equation fails.
    no_vacuum = model.euler_maps(ident, ident, ident, b=4, pg=-2, pf=1, bg=2, bf=-1, epsilon=0)
    return {"positive_cap_does_not_supply_full_vacuum_E_u00": no_vacuum["E_u"][0, 0],
            "same_example_calibrated_E_e00": no_vacuum["E_e"][0, 0],
            "same_example_calibrated_E_v00": no_vacuum["E_v"][0, 0],
            "wrong_other_action_cap_counterexample_R": sp.Rational(3, 2),
            "actual_endpoint_b8_pminus1_H_squared_at_R3over2": sp.Rational(37, 12),
            "singular_TT_hessian": sp.Integer(0),
            "all_zero_links_arbitrary_metric_H_prime_at_zero": sp.Integer(2),
            "actual_deSitter_H_u": sp.Rational(1, 2),
            "kinematic_degenerate_crossing_needs_NEC_violation": -6*sp.Symbol("T", real=True)**2}
