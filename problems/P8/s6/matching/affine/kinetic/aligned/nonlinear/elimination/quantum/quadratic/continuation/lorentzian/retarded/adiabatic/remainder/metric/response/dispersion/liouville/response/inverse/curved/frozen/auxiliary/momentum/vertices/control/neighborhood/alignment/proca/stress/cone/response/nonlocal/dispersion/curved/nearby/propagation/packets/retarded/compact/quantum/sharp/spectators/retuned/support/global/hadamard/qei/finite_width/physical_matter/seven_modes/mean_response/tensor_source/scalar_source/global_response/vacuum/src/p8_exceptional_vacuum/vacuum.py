"""Independent labelled quartic contact and the actual massive vacuum limit."""

from functools import cache
from itertools import permutations

import sympy as sp
from p8_uv import vacuum as original

from . import analytic, family


@cache
def data():
    s, t, w = original.s, original.transfer, original.w
    m = original.mass
    ds, dt, dw = ((z - 2 * m * m) / 2 for z in (s, t, w))
    dots = sp.Matrix(
        [
            [m * m, ds, dt, dw],
            [ds, m * m, dw, dt],
            [dt, dw, m * m, ds],
            [dw, dt, ds, m * m],
        ]
    )
    # Fourier derivatives give i,-1 on first and second derivatives.
    # Each L3/L4 quartic thus carries a minus sign before subtraction.
    L3 = -sum(
        dots[a, b] * dots[b, c] * dots[d, d] for a, b, c, d in permutations(range(4))
    )
    L4 = -sum(
        dots[a, b] * dots[b, c] * dots[c, d] for a, b, c, d in permutations(range(4))
    )
    derivative = sp.factor(L3 - L4)
    g6, g0, lam = sp.symbols(
        "canonical_L3_minus_L4 canonical_potential_quartic canonical_derivative_quartic",
        real=True,
    )
    contact = original.labelled_contact(dots, lam) + g6 * derivative + 24 * g0
    contact = sp.factor(contact)
    b2 = original.forward_coefficient(contact)
    d = analytic.local_jets()
    kap = d["canonical_action_normalization_symbol"]
    coupling = d["canonical_coupling_symbol"]
    n = d["canonical_switch_order_symbol"]
    analytic_contact = contact.subs(
        {
            lam: d["canonical_derivative_quartic_coefficient"],
            g6: d["leading_canonical_L3_minus_L4_coefficient"],
            g0: d["canonical_potential_quartic_coefficient"],
        },
        simultaneous=True,
    )
    analytic_b2 = original.forward_coefficient(analytic_contact)
    canonical = kap * family.vacuum_scalar(kap)
    psi, Y = sp.symbols("vacuum_canonical_scalar vacuum_gradient_square", real=True)
    canonical = canonical.subs(
        {
            family.u: family.u0 + psi / sp.sqrt(kap * family.Z),
            family.X: Y / (kap * family.Z),
        },
        simultaneous=True,
    )
    return {
        "on_shell_mass": m,
        "independent_labelled_L3_minus_L4_contact": derivative,
        "full_scalar_quartic_contact_before_gravity_exchange": contact,
        "scalar_contact_forward_coefficient": b2,
        "analytic_candidate_scalar_contact": sp.factor(analytic_contact),
        "analytic_candidate_scalar_contact_b2": sp.factor(analytic_b2),
        "actual_analytic_calibration_contact_b2": sp.factor(
            analytic_b2.subs(
                {
                    kap: analytic.KAPPA,
                    coupling: analytic.VACUUM_LAMBDA_BAR,
                    n: analytic.ORDER,
                    m: 1,
                }
            )
        ),
        "coscaled_actual_mass_one_contact": sp.factor(
            analytic_contact.subs({kap: n / analytic.FIXED_GAMMA, m: 1})
        ),
        "smooth_candidate_canonical_vacuum": sp.expand(canonical),
        "smooth_candidate_kernel_with_restored_action_normalization": family.vacuum_scalar(
            kap
        ),
        "smooth_candidate_canonical_normalization": "Phi_bar=sqrt(kappa Z)*(u-u0); the bare derivative quartic is kappa*lambda*Z^2*X^2",
        "smooth_candidate_leading_decoupled_b2": 4 * family.lam,
        "finite_gravity_exchange_not_discarded_from_a_claimed_full_amplitude": True,
        "checks": {
            "smooth_actual_canonical_mass_and_interaction": sp.expand(
                canonical - (Y / 2 - family.m**2 * psi**2 / 2 + family.lam * Y**2)
            ),
            "all_incoming_contact_s_t_crossing": sp.expand(
                contact - contact.xreplace({s: t, t: s})
            ),
            "all_incoming_contact_s_w_crossing": sp.expand(
                contact - contact.xreplace({s: w, w: s})
            ),
            "labelled_DHOST_contact_is_minus_three_halves_stu_on_shell": sp.factor(
                (derivative + sp.Rational(3, 2) * s * t * w).subs(w, 4 * m * m - s - t)
            ),
            "decoupling_of_actual_analytic_contact_is_fixed_kessence": sp.factor(
                sp.limit(analytic_contact, kap, sp.oo)
                - original.contact().subs(original.coupling, coupling)
            ),
            "actual_analytic_vacuum_limit_forward_b2": sp.factor(
                sp.limit(analytic_b2, kap, sp.oo) - 4 * coupling
            ),
        },
        "bounds": {
            "actual_analytic_scalar_contact_b2_positive": analytic_b2.subs(
                {
                    kap: analytic.KAPPA,
                    coupling: analytic.VACUUM_LAMBDA_BAR,
                    n: analytic.ORDER,
                    m: 1,
                }
            )
            > 0,
        },
    }
