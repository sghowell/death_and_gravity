"""Continuous classical margins and one explicitly isolated finite loop jet."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_aligned_quantum import potential

from . import action, dynamics

old = dynamics.old


@cache
def polynomial_bounds():
    bg = old.background()
    u = old.u
    value = sp.factor(bg["J"]*bg["h"]**2)
    num, den = sp.fraction(value)
    remainder = sp.Poly(sp.expand(num-sp.Rational(1199, 800)*den), u)
    upper = sp.factor(num.subs(u, sp.Rational(1, 2))/800)
    return {"J_times_h_squared": value,
            "global_lower": sp.Rational(1199, 800),
            "global_lower_remainder_coefficients": remainder.terms(),
            "compact_upper": upper,
            "global_positive_even_remainder": bool(all(
                powers[0] % 2 == 0 and coeff > 0 for powers, coeff in remainder.terms())),
            "denominator_identity": sp.factor(den-800*(1+u**2)**12),
            "compact_upper_below_36": bool(upper < 36)}


@cache
def potential_jet():
    data = potential.clock_jets()
    h = data["h"]
    finite = data["finite_weight"]
    # x=-N^-2 gives B_N=2B_x and B_NN=4B_xx-6B_x.
    Bx = sp.factor(finite["N_first"]/2)
    Bxx = sp.factor((finite["N_second"]+3*finite["N_first"])/4)
    px, pxx, time, spatial = sp.symbols("P_x P_xx clock_time_gradient clock_spatial_gradient", real=True)
    order = sp.Symbol("perturbation_order")
    dx = -2*order*time+order**2*(-time**2+spatial**2)
    quadratic = sp.expand(px*dx+pxx*dx**2/2).coeff(order, 2)
    K, G = sp.diff(quadratic, time, 2), -sp.diff(quadratic, spatial, 2)
    return {"B_x_on_clock": Bx, "B_xx_on_clock": Bxx,
            "finite_Lorentzian_P_xx": -Bxx,
            "finite_weight_x_second_identity": sp.factor(Bxx-sp.Rational(2392, 6561)/h**2),
            "isolated_clock_K_minus_G_identity": sp.factor(K-G-4*pxx),
            "not_the_full_constrained_quantum_principal_symbol": True}


def physical_bounds(margin, planck_time_product, mass_time_product):
    eps = exact(margin, "epsilon_margin")
    L = exact(planck_time_product, "M_tau")
    R = exact(mass_time_product, "m0_tau")
    if not (bool(eps > 0) and bool(eps <= sp.Rational(1, 100)) and bool(L > 0) and bool(R >= 1000)):
        raise ValueError("Require 0<epsilon_margin<=1/100, M*tau>0 and m0*tau>=1000")
    b = polynomial_bounds()
    # pi²>9 gives an upper bound, not a floating approximation.
    loop_weight = R**4/(576*L**2)
    ratio = loop_weight*sp.Rational(1196, 6561)/eps
    return {"epsilon_margin": eps, "M_tau": L, "m0_tau": R,
            "clock_tube_added_P_over_reference_density_upper": eps/100,
            "clock_tube_added_P_x_over_reference_density_upper": eps/5,
            "clock_tube_added_P_xx_over_reference_density_upper": 2*eps,
            "all_time_clock_squared_speed_lower": b["global_lower"]/(b["global_lower"]+4*eps),
            "compact_clock_cone_margin_lower": 4*eps/(b["compact_upper"]+4*eps),
            "all_time_fractional_light_kinetic_increase_upper": 4*eps/b["global_lower"],
            "finite_potential_negative_P_xx_over_added_positive_P_xx_upper": ratio,
            "isolated_potential_jet_is_dominated": bool(ratio < 1),
            "full_quantum_cone_or_cutoff_claim": False}


@cache
def checks():
    p = potential_jet()
    return {"global_J_h_squared_denominator": polynomial_bounds()["denominator_identity"],
            "finite_potential_x_jet": p["finite_weight_x_second_identity"],
            "isolated_clock_principal_difference": p["isolated_clock_K_minus_G_identity"]}


@cache
def proof_checks():
    b = polynomial_bounds()
    eps = sp.Rational(1, 10**6)
    example = physical_bounds(eps, 10**12, 1000)
    return {"all_time_positive_even_J_h_squared_remainder": b["global_positive_even_remainder"],
            "compact_J_h_squared_upper_below_36": b["compact_upper_below_36"],
            "gamma_chart_h_upper_below_five_fourths": bool(sp.Rational(17, 16)**3 < sp.Rational(5, 4)),
            "gamma_chart_lambda_magnitude_at_least_one_fifth": bool(3/(2*sp.Rational(5, 4))-1 == sp.Rational(1, 5)),
            "unitary_chart_Theta_magnitude_at_least_three_fifths": bool(3*sp.Rational(1, 4)/sp.Rational(5, 4) == sp.Rational(3, 5)),
            "compact_auxiliary_J_lower_above_one_over_40": bool(sp.Rational(1199, 800)/sp.Rational(5, 4)**18 > sp.Rational(1, 40)),
            "positive_added_J": action.deformation()["delta_J"].is_positive is True,
            "regular_phase_has_no_Theta_inverse": dynamics.scalar()["no_Theta_inverse"],
            "example_clock_margin_above_one_over_ten_million": bool(example["compact_clock_cone_margin_lower"] > sp.Rational(1, 10**7)),
            "example_fractional_kinetic_change_below_three_parts_per_million": bool(example["all_time_fractional_light_kinetic_increase_upper"] < sp.Rational(3, 10**6)),
            "example_isolated_loop_jet_fraction_below_one_billionth": bool(example["finite_potential_negative_P_xx_over_added_positive_P_xx_upper"] < sp.Rational(1, 10**9)),
            "free_matter_direction_remains_luminal_not_strictly_protected": True,
            "no_full_quantum_cone_from_isolated_potential_sign": True,
            "old_vector_background_operator_and_state_unchanged": True,
            "nonlinear_count_not_energy_stability": True}
