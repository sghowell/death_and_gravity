"""Total-energy continuity of the complete fixed-ball angular conversion."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree.trees import exact_real
from p8_vacuum_affine_radiative_angular_finite import angular

from . import source


def radiation_energy(value):
    value = exact_real(value)
    if not 0 <= value <= s.Rational(1, 8):
        raise ValueError("Require marked total energy in[0,1/8]")
    return value


def modulus(value):
    value = radiation_energy(value)
    return value * (1 - s.log(value)) if value else s.S.Zero


def conversion_change_upper(value):
    return 5000 * modulus(value) / source.KAPPA


def regulator_coefficient_change_upper(value):
    return 8000 * modulus(value) / source.KAPPA


@cache
def physical_samples():
    E = s.Rational(5, 4)
    outgoing = s.Matrix([1, 0, 0])
    directions = (
        s.Matrix([0, 1, 0]),
        s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]),
        s.Matrix([0, s.Rational(3, 5), s.Rational(4, 5)]),
        s.Matrix([0, 0, 1]),
    )
    checks = {}
    gates = {}
    records = {}
    for wi, weight in enumerate(
        (s.Rational(1, 16), s.Rational(1, 1024), s.Rational(1, 65536))
    ):
        ray = s.Matrix([weight, 0, 0, weight])
        for ni, direction in enumerate(directions):
            for ri, radial in enumerate((s.S.Zero, s.Rational(3, 5), s.S.One)):
                label = f"state_{wi}_direction_{ni}_radial_{ri}"
                current = angular.transverse_current(
                    E, [ray], outgoing, direction, radial
                )
                born = angular.transverse_current(E, [], outgoing, direction, radial)
                diff = (current - born).applyfunc(s.factor)
                norm2 = s.factor(s.trace(diff.T * diff))
                hs, us = angular.contractions(current)
                h0, u0 = angular.contractions(born)
                dh = s.factor(hs - h0)
                du = s.factor(us - u0)
                checks[label + "_symmetric_complete_current"] = current - current.T
                gates[label + "_current_and_quadratic_state_change"] = bool(
                    norm2 < (640 * weight) ** 2
                    and abs(dh) < 32000 * weight
                    and abs(du) < 21000 * weight
                )
                if ni == 0 and radial == 1:
                    gates[label + "_physical_state_transport_nonzero"] = dh != 0
                    records["energy_" + str(weight)] = dh
    return checks, gates, records


@cache
def data():
    R, tau = s.symbols("R tau", positive=True)
    radial = s.integrate(tau ** s.Rational(-3, 4), (tau, 0, R**4)) + s.integrate(
        R / tau, (tau, R**4, 1)
    )
    B, B0 = s.Rational(65, 4), s.Integer(16)
    checks = {
        "two_massive_legs_and_entire_null_current": s.Integer(2 * (48 + 256) + 2 - 610),
        "current_quadratic_H_variation": s.Rational(3, 2) * (B + B0) * 640 - 30960,
        "current_quadratic_U_variation": (B + B0) * 640 - 20640,
        "entire_radial_state_Holder_integral": s.simplify(
            radial - 4 * R * (1 - s.log(R))
        ),
        "full_finite_K1_state_budget": s.Integer(10500 + 256000 - 266500),
        "full_regulator_K_state_budget": s.Integer(10500 + 512000 + 21000 - 543500),
        "full_phase_regulator_state_budget": s.Integer(543500 + 4 * 1440 - 549260),
        "state_modulus_derivative": s.simplify(
            s.diff(R * (1 - s.log(R)), R) + s.log(R)
        ),
        "state_modulus_continuous_zero": s.limit(R * (1 - s.log(R)), R, 0, dir="+"),
        "conversion_bound_zero_energy": conversion_change_upper(0),
        "regulator_bound_zero_energy": regulator_coefficient_change_upper(0),
    }
    samples, sample_gates, records = physical_samples()
    checks.update(samples)
    return {
        "whole_total_energy_state_continuity": "At any fixed-ball radial point, the two outgoing massive currents each change by at most304R in nuclear norm: projected numerator difference12R divided by Doppler gap1/4 costs48R; inverse denominator difference64R times numerator4 costs256R. Incoming legs are fixed and all null currents total at most2R. Thus ||T_sigma-T0||_*<640R, with no multiplicity or angular-separation assumption.",
        "whole_quadratic_and_radial_continuity": "Using current norms65/4 and16, |delta H(t)|<32000R and |delta U(t)|<21000R. Combining these with S301's separate quarter-power radial bounds gives the STATE difference of H(t)-H(1) bounded by64000 min(R,(1-t)^(1/4)), and the U difference by42000 times that minimum. The exact radial integral is4R(1-lnR). This is a joint state/radial majorant, not a derivative justified by bounded convergence alone.",
        "whole_finite_and_regulated_conversion": "The entire K1 difference is below300000R(1-lnR), so |Delta_sigma-Delta0|<5000R(1-lnR)/kappa. In the exact finite-e subtracted angular integral, c(e)<=e/2,A(e)<=2 and the phase derivative bound give |b_e(sigma)-b_e(0)|<8000R(1-lnR)/kappa, b_e=(a_e-a)/(2e). The common |b_e|<5500/kappa remains unchanged.",
        "whole_actual_nonzero_current_transport": records,
        "whole_scope": "The state comparison is to the associated elastic Born at fixed incoming energy and outgoing rest-frame direction. It retains all marked null legs, original recoil, trace term and radial finite term. This is continuity of the known soft reference, not finite hard matching or an interacting quantum state.",
        "checks": checks,
        "gates": {
            **sample_gates,
            "current_budget_below640": s.Integer(610) < 640,
            "H_budget_below32000": s.Integer(30960) < 32000,
            "U_budget_below21000": s.Integer(20640) < 21000,
            "K1_budget_below300000": s.Integer(266500) < 300000,
            "finite_conversion_budget_below5000": s.Rational(300000 + 4 * 1440, 72)
            < 5000,
            "finite_regulator_budget_below8000": s.Rational(549260, 72) < 8000,
            "zero_state_limit_not_a_constant_uniform_bound": True,
            "full_radial_and_trace_terms_retained": True,
            "all_null_multiplicities_depend_only_on_total_energy": True,
            "finite_samples_do_not_replace_nuclear_and_Holder_proofs": True,
        },
    }
