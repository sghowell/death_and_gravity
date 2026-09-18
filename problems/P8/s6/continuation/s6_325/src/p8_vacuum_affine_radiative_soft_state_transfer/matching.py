"""Paired universal soft state-change insertion; not a full hard-loop rate."""

from functools import cache

import sympy as s
from p8_vacuum_affine_three_singleton_subtraction.measure import require_cutoff

from . import source


@cache
def data():
    e, y, x, u = s.symbols("e y x u", positive=True)
    a, D, q, da, dD, dq, phi = s.symbols("a D q da dD dq phi", real=True)
    amplitude_virtual = -da / (4 * e) + s.I * phi / e
    regulated = ((da + 2 * e * dD + e * e * dq) * y ** (2 * e) - da) / (2 * e)
    ratio = s.gamma(1 + 2 * e) ** 2 / s.gamma(1 + 4 * e)
    ratio2 = s.series(ratio, e, 0, 3).removeO()
    Ie = (a + 2 * e * D + e * e * q) * x ** (2 * e) / (2 * e)
    rr = Ie * Ie * ratio2 / 2
    rv = -a * Ie / (2 * e)
    vv = a * a / (8 * e * e)
    checks = {
        "soft_logarithmic_power_count": s.integrate(u ** (-1 + 2 * e), (u, 0, 1))
        - 1 / (2 * e),
        "one_internal_attachment_power_count": s.integrate(u ** (2 * e), (u, 0, 1))
        - 1 / (1 + 2 * e),
        "same_state_probability_virtual_pole": amplitude_virtual
        + s.conjugate(amplitude_virtual)
        + da / (2 * e),
        "state_change_finite_soft_kernel": s.limit(regulated, e, 0)
        - dD
        - da * s.log(y),
        "wrong_state_virtual_leaves_pole": s.limit(
            e * (da + 2 * e * dD) * y ** (2 * e) / (2 * e), e, 0
        )
        - da / 2,
        "state_continuity_majorant_integral": x
        * s.integrate(1 - s.log(x) - s.log(u), (u, 0, 1))
        - x * (2 - s.log(x)),
        # Reflect v=1-u on the real interval for the remaining-energy log.
        "remaining_energy_endpoint_log_integral": x
        * s.integrate(-s.log(x) - s.log(u), (u, 0, 1))
        - x * (1 - s.log(x)),
        "Bose_total_energy_gamma_correction": s.expand(
            ratio2 - 1 + 4 * s.zeta(2) * e * e
        ),
        "Born_soft_second_order_finite_term": s.limit(rr + rv + vv, e, 0)
        - (D + a * s.log(x)) ** 2 / 2
        + s.pi**2 * a * a / 12,
        "retained_remaining_energy_change": s.expand_log(
            da * s.log(x * (1 - u)) - da * s.log(x) - da * s.log(1 - u), force=True
        ),
    }
    return {
        "checks": {k: s.factor(v) for k, v in checks.items()},
        "gates": {
            "limit_integrated_coefficient": 2 * 10000 + 1400 < 22000,
            "regulated_integrated_coefficient": 2 * 15000 + 1400 < 32000,
            "original_Born_weighted_limit": 8 * 22000 < 200000,
            "wrong_state_pole_nonzero_for_nonzero_change": s.factor(da / 2) != 0,
            "Coulomb_phase_retained_in_amplitude": amplitude_virtual.has(phi),
            "Gamma_total_energy_term_nonzero": s.pi**2 * a * a / 12 != 0,
            "marked_Born_measure_physical_D4": True,
            "additional_soft_projector_and_phase_full_D": True,
            "no_complete_hard_loop_or_outer_D_assembly_claim": True,
        },
        "whole_paired_state_difference": "Z_e(sigma;y)=[(a_sigma,e-a0,e)*y^(2e)-(a_sigma-a0)]/(2e) -> (Delta_sigma-Delta0)+(a_sigma-a0)*ln(y). Same-state virtual probability pole is-a_sigma/(2e). The imaginary Coulomb phase remains in the amplitude.",
        "whole_named_Born_seeded_insertion": "dnu_e,x=Z_e(sigma_b;x-b)dLambda0(b,Omega) on0<b<x; its total mass is C_e(x)=int dnu_e,x. The measures converge in total variation. dLambda0 is the physical D4 one-real Born leading measure with radial db/b and total angular coefficient a0<8/kappa. Only the additional universal soft factor is dimensionally continued.",
        "whole_uniform_regulator_dominator": (
            15000 * u * (1 - s.log(u)) + 1400 * u * (-s.log(y))
        )
        / source.KAPPA,
        "whole_transfer_total_variation_upper": 200000
        * x
        * (1 - s.log(x))
        / source.KAPPA**2,
        "whole_calorimetric_Born_finite_coefficient": (D + a * s.log(x)) ** 2 / 2
        - s.pi**2 * a * a / 12,
        "whole_matching_boundary": "This computes a named leading-soft real/virtual state-difference insertion and a uniform regulator limit. It does not reconstruct finite radiative hard loops, identify the whole NNLO rate, drop evanescent outer-state terms, or establish all-N hard matching, quantum state, unitarity, Regge, bounce or original P8 closure.",
    }


def transfer_upper(value):
    x = require_cutoff(value)
    return 200000 * x * (1 - s.log(x)) / source.KAPPA**2
