"""Actual changed quantum readouts; no automatic old tadpole or response transfer."""
from functools import cache

import sympy as sp
from p8_vector_clock_matching import bounds as old_local
from p8_vector_clock_matching import counterterms
from p8_vector_dimensional import local
from p8_vector_state import wkb
from p8_vector_variation import modes, readouts


@cache
def matrices():
    zero={modes.alpha:0,modes.beta:0,readouts.alpha2:0,readouts.beta2:0}
    result={}
    for kind,prior in readouts.matrices().items():
        old=prior["Hamiltonian_quadratic_matrix"]
        new=old.subs(zero)
        g2=modes.canonical()[kind]["g_squared"].subs({modes.N:1,modes.am:1,modes.bm:1})
        C=sp.diag(1/sp.sqrt(g2),sp.sqrt(g2))
        row={"new_Hamiltonian_matrix":new,
             "old_Hamiltonian_matrix":old,
             "on_clock_Hamiltonian_difference":readouts.clean((new-old).subs(modes.N,1))}
        for observable in ("energy","pressure"):
            physical=prior[observable+"_matrix"]
            before=readouts.clean(C.T*physical*C)
            after=readouts.clean(before.subs(zero))
            row["old_canonical_"+observable]=before
            row["new_canonical_"+observable]=after
            row["canonical_"+observable+"_change"]=readouts.clean(after-before)
            row["new_"+observable+"_contact"]=readouts.clean(prior[observable+"_contact_matrix"].subs(zero))
        row["old_minus_new_first_lapse_operator_jet"]=readouts.clean(
            (old-new).diff(modes.N).subs(modes.N,1))
        row["old_minus_new_second_lapse_operator_jet"]=readouts.clean(
            (old-new).diff(modes.N,2).subs(modes.N,1))
        result[kind]=row
    return result


@cache
def local_coefficients():
    H,u=wkb.background()["H"],wkb.u
    mapping={local.alpha:0,local.beta:0,local.ell:0}
    mapping.update({symbol:sp.diff(H,u,n) for n,symbol in enumerate(
        (local.H,local.Hd,local.Hdd,local.Hddd,local.Hdddd))})
    return {n:{name:sp.factor(value[name].subs(mapping)) for name in ("energy","pressure")}
            for n,value in counterterms.finite_coefficients().items()}


@cache
def checks():
    a,q,m2=modes.scale,modes.q,modes.mass2
    z=q/(q+m2)
    target={"transverse":sp.diag(-m2*modes.beta,0)/a**3,
            "longitudinal":sp.diag(-(q+m2)*modes.beta,modes.alpha*z)/a**3}
    out={}
    for kind,row in matrices().items():
        out[kind+"_same_on_clock_vector_operator"]=row["on_clock_Hamiltonian_difference"]
        out[kind+"_actual_changed_quantum_energy_readout"]=readouts.clean(row["canonical_energy_change"]-target[kind])
        out[kind+"_new_quantum_energy_is_ordinary_Proca"]=readouts.clean(
            row["new_canonical_energy"]-sp.diag(q+m2,1)/a**3)
        out[kind+"_unchanged_on_clock_pressure_readout"]=row["canonical_pressure_change"]
    H,u=wkb.background()["H"],wkb.u
    ordinary={0:-sp.Rational(5,2),1:-10*H**2,
              2:2*(6*H**2*sp.diff(H,u)+2*H*sp.diff(H,u,2)-sp.diff(H,u)**2)}
    for n,row in local_coefficients().items():
        out[f"actual_constant_mass_finite_local_energy_{n}"]=sp.factor(row["energy"]-ordinary[n])
        out[f"unchanged_covariantly_matched_local_pressure_{n}"]=sp.factor(row["pressure"]-old_local.actual_coefficients()[n]["pressure"])
        out[f"ordinary_local_quantum_Ward_conservation_{n}"]=sp.factor(sp.diff(row["energy"],u)+3*H*(row["energy"]+row["pressure"]))
    return out


@cache
def boundary():
    h=(1+wkb.u**2)**3
    return {"actual_old_first_mass_lapse_jets":{"temporal":4/(9*h),"spatial":28/(81*h)},
            "actual_new_first_and_higher_mass_lapse_jets":sp.Integer(0),
            "new_constant_mass_local_finite_coefficients":local_coefficients(),
            "same_original_clock_modes_and_initial_covariance":True,
            "state_not_reselected":True,
            "old_fixed_scalar_profiles_in_this_comparison":"Retained as fixed old coefficients, not silently changed. Their previous cancellation is not asserted for the new energy readout. A replacement must be computed once from the new finite Proca stress and separately budgeted.",
            "required_profile_change_if_new_finite_stress_is_constructed":"The on-clock pressure mode readout and matched local pressure are unchanged. With Delta rho = rho_new-rho_old, the required new fixed scalar profile difference would be Delta rho(phi)*(x+1)/2 on the flat clock tube, with the same fixed cutoff. Its quantitative norm is not supplied here.",
            "quantum_response_boundary":"Agreement of the vector operator along the original clock does not imply agreement of its functional metric derivatives. The first and second lapse operator jets and energy contacts change. Old Gaussian light response, frozen pole brackets and exact background tadpole cancellation are not transferred. No newly quantum-corrected bounce is asserted."}


@cache
def gates():
    rows=matrices()
    return {"both_actual_energy_readouts_change":all(row["canonical_energy_change"]!=sp.zeros(2) for row in rows.values()),
            "both_first_lapse_operator_jets_change":all(row["old_minus_new_first_lapse_operator_jet"]!=sp.zeros(2) for row in rows.values()),
            "both_second_lapse_operator_jets_change":all(row["old_minus_new_second_lapse_operator_jet"]!=sp.zeros(2) for row in rows.values()),
            "new_fixed_background_operator_not_declared_equal_off_clock":True,
            "old_tadpole_and_frozen_quantum_diagnostic_not_transferred":True,
            "no_new_quantum_stress_profile_norm_claim":True}
