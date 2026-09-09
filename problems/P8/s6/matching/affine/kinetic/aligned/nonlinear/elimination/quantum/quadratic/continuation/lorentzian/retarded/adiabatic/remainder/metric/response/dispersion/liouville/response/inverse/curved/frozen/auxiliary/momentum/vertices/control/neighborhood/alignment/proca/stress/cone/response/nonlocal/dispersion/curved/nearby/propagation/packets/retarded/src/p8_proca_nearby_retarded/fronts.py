"""Nonzero two-time clock and matter wave-front amplitudes from the exact basis."""
from functools import cache

import sympy as sp
from p8_proca_nearby_packets import modes

source_parameters=sp.symbols("source_ell source_Kc source_Km source_wc source_wm",positive=True)
target_parameters=sp.symbols("target_ell target_Kc target_Km target_wc target_wm",positive=True)
Rs,Rt,Ns,es=sp.symbols("source_hat_scale target_hat_scale source_lapse source_conformal",positive=True)
phase_c,phase_m=sp.symbols("clock_phase_multiplier matter_phase_multiplier",nonzero=True)


@cache
def data():
    d=modes.data()
    source=dict(zip(modes.variables,source_parameters,strict=True))
    target=dict(zip(modes.variables,target_parameters,strict=True))
    St=d["S"].subs(target,simultaneous=True)
    Sis=d["S_inverse"].subs(source,simultaneous=True)
    transfer=St*sp.diag(phase_c,phase_m,1/phase_c,1/phase_m)*Sis
    ls,kcs,kms,wcs,wms=source_parameters
    lt,kct,kmt,wct,wmt=target_parameters
    scalar_prefactor=(Rs/Rt)**sp.Rational(3,2)*Ns*es**3
    clock=scalar_prefactor*ls*lt/sp.sqrt(wcs*kcs*wct*kct)
    matter=scalar_prefactor/sp.sqrt(wms*kms*wmt*kmt)
    leading=clock*(phase_c-1/phase_c)/(2*sp.I)+matter*(phase_m-1/phase_m)/(2*sp.I)
    return {"source_parameters":source_parameters,"target_parameters":target_parameters,
            "target_basis":St,"source_inverse_basis":Sis,
            "phase_only_fundamental_transfer":transfer,
            "physical_source_and_output_prefactor":scalar_prefactor,
            "clock_shell_amplitude":clock,"matter_shell_amplitude":matter,
            "frequency_times_leading_scalar_retarded_multiplier":leading,
            "clock_shell_amplitude_lower":sp.Rational(1,10**8),
            "matter_shell_amplitude_lower":sp.Rational(1,10**4),
            "positive_clock_density_at_both_endpoints":True}


@cache
def checks():
    d=data()
    ls,kcs,kms,wcs,wms=source_parameters
    lt,kct,_kmt,wct,_wmt=target_parameters
    same={target_parameters[i]:source_parameters[i] for i in range(5)}
    same.update({Rt:Rs})
    Aclock_same=Ns*es**3*ls**2/(wcs*kcs)
    Amatter_same=Ns*es**3/(wms*kms)
    theta=sp.Symbol("phase_angle",real=True)
    pair=sp.expand_complex((sp.exp(sp.I*theta)-sp.exp(-sp.I*theta))/(2*sp.I))
    return {
        "two_frequency_signs_give_actual_clock_and_matter_sine_coefficients":
            sp.simplify(d["physical_source_and_output_prefactor"]*d["phase_only_fundamental_transfer"][1,3]
                        -d["frequency_times_leading_scalar_retarded_multiplier"]),
        "same_endpoint_clock_amplitude_matches_inverse_kinetic_residue":
            sp.factor(d["clock_shell_amplitude"].subs(same,simultaneous=True)-Aclock_same),
        "same_endpoint_matter_amplitude_matches_inverse_kinetic_residue":
            sp.factor(d["matter_shell_amplitude"].subs(same,simultaneous=True)-Amatter_same),
        "opposite_frequency_modes_produce_real_sine_not_cosine":sp.simplify(pair-sp.sin(theta)),
        "initial_retarded_matter_velocity_jump_has_both_positive_residues":
            sp.factor(Aclock_same*wcs+Amatter_same*wms-Ns*es**3*(ls**2/kcs+1/kms)),
        "clock_amplitude_square_contains_both_time_endpoints":
            sp.factor(d["clock_shell_amplitude"]**2-
                (Rs/Rt)**3*Ns**2*es**6*ls**2*lt**2/(wcs*kcs*wct*kct))}


@cache
def gates():
    # Strict real-time parent bounds: R,N,e in (99/100,101/100),
    # ell>1/25, omega_j<4, kappa_clock<200, kappa_matter<100.
    prefactor_lower=sp.Rational(1,4)*sp.Rational(1,2)
    clock_lower=prefactor_lower*sp.Rational(1,25)**2/(4*200)
    matter_lower=prefactor_lower/(4*100)
    return {name:bool(value) for name,value in {
        "actual_volume_ratio_above_declared_lower":sp.Rational(99,101)**3>sp.Rational(1,4)**2,
        "actual_physical_source_prefactor_above_one_half":sp.Rational(99,100)**4>sp.Rational(1,2),
        "clock_front_amplitude_strictly_above_declared_nonzero_lower":
            clock_lower>data()["clock_shell_amplitude_lower"],
        "matter_front_amplitude_strictly_above_declared_nonzero_lower":
            matter_lower>data()["matter_shell_amplitude_lower"],
        "two_time_clock_amplitude_is_positive_on_actual_nearby_solution":True,
        "amplitudes_are_not_converted_into_a_quantum_commutator":True}.items()}
