"""Finite-time, hard-transfer, connected seven-mode tree operator bounds."""
from functools import cache

import sympy as sp
from p8_coupled_vertices import majorant

from . import other_modes, window

SPECIES=7
WICK_FACTOR=720
EXCHANGE_FACTOR=3*SPECIES*2*2
TARGET=sp.Rational(1,1000)


@cache
def data():
    domain=window.data()
    if not all(bool(value) for value in other_modes.gates().values()):
        raise ValueError("The tensor or unchanged-state vector phase bound failed")
    raw=majorant.bound(domain["common_raw_phase_seed_after_center_Fourier_conversion"],
                       domain["Hamiltonian_derivative_bound"],domain["York_inverse_transfer_bound"])
    h3=sp.Rational(raw["cubic_labelled_kernel_upper"])
    h4=sp.Rational(raw["quartic_labelled_kernel_upper"])
    # S6.77 is per instantaneous local spatial volume. The fixed-center
    # Fourier Hamiltonian has the additional (a/a_center)^3<=8.
    B3,B4=8*WICK_FACTOR*h3,8*WICK_FACTOR*h4
    phase_measure=SPECIES*(4*window.UPPER_K+1)**3
    schur=8*phase_measure**3
    T=domain["full_window"]
    C3=schur*T*B3
    C4=schur*(T*B4+EXCHANGE_FACTOR*T*T*B3*B3)
    scale=sp.Integer(1)
    while C3/scale>TARGET or C4/scale**2>TARGET:
        scale*=10
    return {"domain":domain,"propagating_mode_columns":SPECIES,
            "raw_Hamiltonian_bound":majorant.serialize_bound(raw),
            "Wick_majorant":WICK_FACTOR,"fixed_center_background_volume_upper":8,
            "connected_exchange_contraction_majorant":EXCHANGE_FACTOR,
            "momentum_species_measure_upper":phase_measure,"fixed_total_momentum_Schur_upper":schur,
            "cubic_mode_kernel_upper":B3,"quartic_mode_kernel_upper":B4,
            "cubic_transition_numerator":C3,"quartic_connected_tree_numerator":C4,
            "sufficient_named_M_tau":scale,"sufficient_named_M_tau_power10":len(str(scale))-1,
            "cubic_block_bound_at_named_scale":C3/scale,
            "quartic_connected_tree_bound_at_named_scale":C4/scale**2,
            "target_block_norm":TARGET,
            "domain_description":"Finite-time 1-to-2 and 2-to-1 cubic blocks and the hard-masked connected 2-to-2 tree block through H4 and H3 squared, on fixed-total-momentum fibers of R^3. All seven modes, including on-shell massive Proca modes, are retained.",
            "not_a_light_only_heavy_elimination_domain":True,
            "not_a_Wilsonian_or_all_orders_cutoff":True,
            "not_the_previous_M_tau_10_to_24_example":True}


@cache
def scale_checks():
    L=sp.Symbol("positive_M_tau",positive=True)
    q,p=sp.symbols("raw_coordinate raw_momentum")
    out={}
    for n in (2,3,4):
        expression=sum(sp.Symbol(f"coefficient_{n}_{j}")*q**j*p**(n-j) for j in range(n+1))
        out[f"full_phase_homothety_degree_{n}"]=sp.expand(
            L**2*expression.subs({q:q/L,p:p/L},simultaneous=True)-L**(2-n)*expression)
    return out
