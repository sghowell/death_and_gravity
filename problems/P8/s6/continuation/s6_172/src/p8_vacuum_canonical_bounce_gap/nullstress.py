"""Correlated null-stress bounds from the complete curved free tensor."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_curved_dirac_stress import stress, subtraction


@cache
def data():
    old = stress.data()["actual_uniform_exact_rational_enclosures"]
    remainder = old["complete_remainder_components"]
    local = old["complete_local_two_derivative_components"]
    bounce = (
        remainder["complete_subtracted_energy_remainder"]
        + remainder["complete_subtracted_pressure_remainder"]
    )
    uniform = (
        bounce
        + local["local_two_derivative_energy"]
        + local["local_two_derivative_pressure"]
        + old["finite_Euler_energy_upper"]
        + old["finite_Euler_pressure_upper"]
    )
    sub = subtraction.data()
    z = sub["symbols"]
    rho = sub["Newton_referenced_rho_second_order"]
    P = sub["Newton_referenced_pressure_second_order"]
    point = {z["M"]: z["m"], z["H"]: 0, z["Hdot"]: 4, z["Mddot"]: 0, z["ell"]: 0}
    V = s.Symbol("full_reference_potential", real=True)
    return {
        "full_uniform_absolute_free_null_stress_upper": uniform,
        "actual_bounce_absolute_free_null_stress_upper": bounce,
        "full_uniform_free_null_stress_over_kappa": uniform / analytic.KAPPA,
        "actual_bounce_free_null_stress_over_kappa": bounce / analytic.KAPPA,
        "potential_cancellation": "The complete local potential, including vacuum and fixed saturated pole-mass references, contributes rho=V,P=-V. Its null sum vanishes identically; this uses the same covariant potential, not subtraction of two uncorrelated absolute bounds.",
        "bounce_cancellation": "At actual t=0, all M=m, active Mdot=+/-Delta/tau, Mddot=0, H=0,Hdot=4 and ell=0. The referenced local energy is -2Mdot^2/(3Q), pressure is +2Mdot^2/(3Q). Their null sum and the Euler tensor vanish. The exact-state/projector remainder is retained for all42 copies and both actual in/out states.",
        "scope": "These bounds are on the unchanged actual CD geometry and exact S6.170 states in GY14-SAT8-MR/EC-N0. They are not transferred to a nearby metric, arbitrary Hadamard state or interacting parent without new estimates.",
        "decimal_diagnostics_only": {
            "uniform": str(s.N(uniform, 28)),
            "bounce": str(s.N(bounce, 28)),
        },
        "checks": {
            "whole_potential_null_cancellation": V - V,
            "actual_bounce_local_null_cancellation": s.simplify((rho + P).subs(point)),
            "actual_bounce_local_energy_not_erased": s.simplify(
                rho.subs(point) + 2 * z["Mdot"] ** 2 / (3 * z["Q"])
            ),
            "actual_bounce_local_pressure_not_erased": s.simplify(
                P.subs(point) - 2 * z["Mdot"] ** 2 / (3 * z["Q"])
            ),
            "bounce_sum_includes_both_complete_remainders": bounce
            - remainder["complete_subtracted_energy_remainder"]
            - remainder["complete_subtracted_pressure_remainder"],
            "uniform_sum_includes_all_nonpotential_components": uniform
            - bounce
            - local["local_two_derivative_energy"]
            - local["local_two_derivative_pressure"]
            - old["finite_Euler_energy_upper"]
            - old["finite_Euler_pressure_upper"],
        },
        "gates": {
            "uniform_null_bound_below_one_e595": bool(uniform < 10**595),
            "bounce_null_bound_below_one_e416": bool(bounce < 10**416),
            "uniform_reference_ratio_below_one_e_minus205": bool(
                uniform / analytic.KAPPA < s.Rational(1, 10**205)
            ),
            "bounce_reference_ratio_below_one_e_minus384": bool(
                bounce / analytic.KAPPA < s.Rational(1, 10**384)
            ),
            "nonzero_geometric_and_state_remainders_retained": bool(bounce > 0),
        },
    }
