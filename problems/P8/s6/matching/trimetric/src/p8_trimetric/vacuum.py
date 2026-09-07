"""Exact chosen-parent vacuum prerequisite and restricted effective potential.

The nonexistence theorem uses the full e/f equations, without eliminating
u, dividing by B, choosing an FP mass formula, or truncating matter.
An added beta4 pair is deliberately absent here and treated separately.
"""

from functools import cache

import sympy as sp

from . import model


@cache
def derive():
    b, pg, pf, eps = model.B, model.PG, model.PF, model.EPS
    vacuum_energy = sp.Symbol("V", real=True)
    effective_b = b+eps*vacuum_energy/2
    # Valid auxiliary elimination additionally requires effective_b!=0
    # and invertible Q=pg*e+pf*v. The undivided no-flat theorem does not.
    c = -54*pg**4/b**3
    ratio = pf/pg
    betas = tuple(-54*pg**(4-n)*pf**n/b**3 for n in range(5))
    # G0 is the leading h0=u0.T*eta*u0, not the full matter metric h.
    return {"B": b, "p_g": pg, "p_f": pf, "epsilon": eps, "V": vacuum_energy,
            "B_with_constant_matter_potential": effective_b,
            "geometric_C": c, "geometric_ratio": ratio,
            "effective_beta_values": betas,
            "effective_a": -3*pg/b, "effective_b": -3*pf/b,
            "source_paper_a": 3*pg/b, "source_paper_b": 3*pf/b,
            "vacuum_effective_potential_over_detQ": -54/b**3,
            "full_constant_V_action_over_detQ": 54/effective_b**3,
            "leading_physical_volume_over_detQ": 81/b**4,
            "leading_action_over_h0_volume": 2*b/3,
            "truncated_constant_V_action_over_h0_volume": 2*b/3-eps*vacuum_energy,
            "exact_constant_V_action_over_h0_volume": (2*b/3)*(b/effective_b)**3,
            "naive_cancel_epsilon_V_over_B": sp.Rational(2, 3),
            "exact_action_at_naive_cancel_over_B_h0_volume": sp.Rational(9, 32)}


@cache
def checks():
    d = derive()
    b, eps, potential = (d[key] for key in ("B", "epsilon", "V"))
    qdiag = sp.diag(*sp.symbols("Q0:4", nonzero=True))
    u = -3*qdiag/b
    literal = -2*u.det()*(b+sp.trace(u.inv()*qdiag))
    effective_b = d["B_with_constant_matter_potential"]
    shifted_u = -3*qdiag/effective_b
    exact_constant = -2*shifted_u.det()*(effective_b+sp.trace(shifted_u.inv()*qdiag))
    beta = d["effective_beta_values"]
    ratio = d["geometric_ratio"]
    cancel = {potential: 2*b/(3*eps)}
    y = sp.Symbol("positive_proportional_ratio", positive=True)
    return {
        "literal_vacuum_auxiliary_stationary_action": sp.cancel(literal-54*qdiag.det()/b**3),
        "literal_constant_V_exact_stationary_action": sp.cancel(exact_constant-54*qdiag.det()/effective_b**3),
        "effective_physical_volume_dictionary": sp.cancel(u.det()-81*qdiag.det()/b**4),
        "geometric_middle_minor": sp.cancel(beta[1]*beta[3]-beta[2]**2),
        "geometric_first_coefficient": sp.cancel(beta[0]-d["geometric_C"]),
        "geometric_last_coefficient": sp.cancel(beta[4]-d["geometric_C"]*ratio**4),
        "fixed_composite_weight_ratio": sp.cancel(d["effective_b"]/d["effective_a"]-ratio),
        "g_flat_lapse_geometric_polynomial": sp.cancel(beta[0]+3*y*beta[1]+3*y*y*beta[2]+y**3*beta[3]
                                                       -d["geometric_C"]*(1+ratio*y)**3),
        "f_flat_lapse_geometric_polynomial": sp.cancel(beta[1]+3*y*beta[2]+3*y*y*beta[3]+y**3*beta[4]
                                                       -d["geometric_C"]*ratio*(1+ratio*y)**3),
        "truncated_vacuum_cancellation": sp.cancel(d["truncated_constant_V_action_over_h0_volume"].subs(cancel)),
        "full_action_does_not_share_truncated_cancellation": sp.cancel(
            d["exact_constant_V_action_over_h0_volume"].subs(cancel)/b-sp.Rational(9, 32)),
        "constant_V_first_order_metric_rescaling": sp.diff((b/effective_b)**2, eps).subs(eps, 0)+potential/b,
    }


def controls():
    identity = sp.eye(4)
    original = model.euler_maps(identity, identity, identity, b=-6, pg=1, pf=1, epsilon=0)
    one_link = model.euler_maps(identity, identity, identity, b=-3, pg=1, pf=0, epsilon=0)
    bzero = model.euler_maps(identity, identity, identity, b=0, pg=1, pf=1, epsilon=0)
    disconnected = model.euler_maps(identity, identity, identity, b=0, pg=0, pf=0, epsilon=0)
    # Inserting a scalar vacuum potential can repair E_u at a chosen
    # point but leaves E_e,E_v unchanged. The original fixture already
    # has E_u=0, so dropping both dynamical metric equations would pass.
    return {"original_regular_g_flat_Euler_00": original["E_e"][0, 0],
            "original_regular_f_flat_Euler_00": original["E_v"][0, 0],
            "dropping_e_f_equations_false_vacuum_u_Euler_00": original["E_u"][0, 0],
            "only_one_nonzero_link_still_fails_g_Euler_00": one_link["E_e"][0, 0],
            "B_zero_still_fails_g_Euler_00": bzero["E_e"][0, 0],
            "both_links_zero_B_zero_flat_equations": sum(disconnected[key][0, 0] for key in ("E_e", "E_v", "E_u")),
            "pure_beta2_middle_minor": sp.Integer(-1),
            "old_beta2_with_vacuum_shift_middle_minor": sp.Rational(-1, 4),
            "formal_cancel_is_not_small_epsilon_V_over_B": sp.Rational(2, 3),
            "actual_action_at_formal_cancel_over_B_h0_volume": sp.Rational(9, 32)}
