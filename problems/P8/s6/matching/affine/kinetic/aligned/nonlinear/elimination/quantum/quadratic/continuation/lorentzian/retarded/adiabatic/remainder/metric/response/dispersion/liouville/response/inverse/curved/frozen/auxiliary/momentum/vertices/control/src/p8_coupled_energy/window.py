"""A concrete high-frequency, short-window scalar column bound."""
from functools import cache

import sympy as sp

from . import bounds

LOWER_K=10**12
UPPER_K=10**13
HALF_WINDOW=sp.Rational(1,4*10**9)
GROWTH=10**9
RAW_SEED=10**16*UPPER_K**2


@cache
def data():
    charts={name:bounds.chart_bounds(name) for name in ("unitary","gamma")}
    qmin=max(row["sufficient_q_lower"] for row in charts.values())
    growth=max(row["energy_logarithmic_growth_upper"] for row in charts.values())
    beta=max(row["coefficient_bounds"][name]["operator_upper"] for row in charts.values()
             for name in ("beta_leading","beta_remainder","antisymmetric_mixing"))
    # Initial proper frequency lies in [LOWER_K/2,4*UPPER_K].
    # The dimensionless scalar free columns satisfy E_initial<1e9*kappa,
    # then E<2e9*kappa on this window, with the canonical data of energy.py.
    margins={
        "both_chart_energy_growth_below_named_constant":GROWTH-growth,
        "all_beta_parts_below_common_bound":10**5-beta,
        "whole_window_kinetic_and_potential_high_q":sp.Rational(LOWER_K**2,4)-qmin,
        "initial_antisymmetric_velocity_bound":sp.Rational(LOWER_K,2)-1000*beta,
        "positive_width_inside_both_chart_overlap_margins":sp.Rational(1,160)-HALF_WINDOW,
        "many_oscillatory_radians_in_half_window":LOWER_K*HALF_WINDOW-100,
        "scale_drift_exponent_below_half":sp.Rational(1,2)-8*HALF_WINDOW,
        "energy_growth_exponent_at_most_half":sp.Rational(1,2)-2*GROWTH*HALF_WINDOW,
        "coordinate_norm_squared_below_64":64-sp.Rational(32*10**12,LOWER_K),
        "velocity_squared_below_declared_linear_upper":(10**7*UPPER_K)**2-16*10**12*UPPER_K,
        "phase_momentum_density_below_declared_upper":
            10**12*UPPER_K**2-(10**11*UPPER_K+8*(16*UPPER_K**2+1)*10**5),
        "gamma_swap_nonzero_q_lower":sp.Rational(LOWER_K**2,4)-1,
        "new_matter_boundary_phase_bounds_fit_before_Fourier_conversion":
            2*10**12*UPPER_K**2-(sp.Rational(23,20)*10**12*UPPER_K**2+3),
        "constant_center_Fourier_and_TT_factors_fit_raw_seed":
            RAW_SEED-8*10**12*UPPER_K**2}
    if any(value<0 for value in margins.values()) or any(value==0 for key,value in margins.items()
            if key!="energy_growth_exponent_at_most_half"):
        raise ValueError("The concrete scalar free-column domain failed")
    return {"charts":charts,"external_fixed_center_momentum_band":(LOWER_K,UPPER_K),
            "nonempty_proper_subset_fixed_center_transfer_lower":LOWER_K,
            "largest_internal_tree_fixed_center_momentum":2*UPPER_K,
            "all_time_proper_mode_frequency_bounds":(sp.Rational(LOWER_K,2),4*UPPER_K),
            "Hamiltonian_derivative_bound":8*UPPER_K,
            "York_inverse_transfer_bound":sp.Rational(2,LOWER_K),
            "half_window":HALF_WINDOW,"full_window":2*HALF_WINDOW,
            "scalar_energy_growth_upper":GROWTH,"scalar_energy_inflation_upper":2,
            "initial_energy_upper_over_kappa":10**9,
            "evolved_coordinate_norm_upper":8,
            "evolved_velocity_norm_upper":10**7*UPPER_K,
            "original_chart_momentum_density_norm_upper":10**12*UPPER_K**2,
            "new_matter_boundary_phase_upper":2*10**12*UPPER_K**2,
            "common_raw_phase_seed_after_center_Fourier_conversion":RAW_SEED,
            "margins":margins,
            "mode_choice":"Scalar canonical initial columns use an arbitrary real Cholesky root of the actual positive kinetic matrix and the symmetric momentum boundary. The antisymmetric initial-velocity contribution is retained. This does not reselect the vector Gaussian state.",
            "window_domain":"Any interval of this full length contained in I. Choose gamma if its center has |u|<=19/80, unitary otherwise. The half-width is below 1/160, keeping the whole interval in the corresponding certified chart."}
