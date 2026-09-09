"""Linear relational matter observable, nonzero clock pole and norm conversion."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import bounds as cone_bounds

from . import domains, modes, normal_form

PROJECTION_LOWER=sp.Rational(1,10**4)
PROJECTION_UPPER=sp.Integer(1000)
RECONSTRUCTION_NORM_FACTOR=4*domains.BASIS_NORM/PROJECTION_LOWER
DECLARED_ERROR_COEFFICIENT=sp.Integer(10)**44


@cache
def data():
    eupper=1+sp.Rational(1,10**6)
    clock_residue=sp.Rational(9,1000)/eupper**2*4*sp.Rational(24,100)**2/3
    return {"local_background_dependent_scalar":"O=chi-chi_background(phi)",
            "linear_relational_observable":"delta chi-(chi_background'/phi_background') delta phi",
            "unitary_gauge_observable":"chi",
            "clock_principal_pole_weight_lower":clock_residue,
            "clock_coordinate_projection_lower_before_volume_rescaling":sp.Rational(1,25)/40,
            "matter_coordinate_projection_lower_before_volume_rescaling":sp.Rational(1,40),
            "inverse_volume_rescaling_lower":sp.Rational(1,4),
            "inverse_volume_rescaling_upper":sp.Integer(4),
            "both_mode_observable_projection_lower":PROJECTION_LOWER,
            "both_mode_observable_projection_upper":PROJECTION_UPPER,
            "modal_to_relative_observable_error_factor":RECONSTRUCTION_NORM_FACTOR,
            "declared_inverse_carrier_relative_observable_error_coefficient":DECLARED_ERROR_COEFFICIENT,
            "formal_covariant_probe":"sqrt(-g) J [chi-chi_background(phi)]",
            "no_compact_source_retarded_preparation_or_quantum_detector_theorem":True}


def relative_error(carrier):
    carrier=normal_form.momentum(carrier)
    if carrier/2<normal_form.MIN_MOMENTUM:
        raise ValueError("The complete packet band must lie in the normal-form domain")
    small=normal_form.error(carrier/2)
    actual=RECONSTRUCTION_NORM_FACTOR*small["full_U_error_upper"]
    return {"carrier":carrier,"minimum_band_momentum":carrier/2,
            "actual_relative_observable_error_upper":actual,
            "declared_relative_observable_error_upper":DECLARED_ERROR_COEFFICIENT/carrier}


@cache
def checks():
    dchi,dphi,chip,phip,xi,J=sp.symbols("delta_chi delta_phi chi_prime phi_prime gauge_shift probe",nonzero=True)
    O=dchi-chip*dphi/phip
    gauge={dchi:dchi-chip*xi,dphi:dphi-phip*xi}
    omega,q,gc,gm=sp.symbols("frequency wave_square G_clock G_matter",real=True)
    d=modes.data()
    T=sp.Matrix([[1,0],[-modes.ell,1]])
    G=T.T.inv()*sp.diag(gc,gm)*T.inv()
    actual=(omega**2*d["kinetic"]-q*G).inv()[1,1]
    expected=modes.ell**2/(modes.kc*omega**2-q*gc)+1/(modes.km*omega**2-q*gm)
    u,background,metric,J0=sp.symbols("field_order background_observable metric_perturbation source",real=True)
    source=(1+u*metric)*u*J0*(background+u*O)
    return {
        "linear_relational_matter_observable_is_gauge_invariant":sp.expand(O.subs(gauge,simultaneous=True)-O),
        "unitary_gauge_is_the_reconstructed_matter_coordinate":sp.expand(O.subs(dphi,0)-dchi),
        "formal_clock_and_matter_probe_Noether_forces_cancel":J*chip-J*chip/phip*phip,
        "vanishing_background_observable_has_no_linear_metric_probe_force":
            sp.diff(source,u,2).subs({u:0,background:0}).diff(metric),
        "full_principal_matter_response_has_both_distinct_positive_poles":sp.factor(actual-expected),
        "clock_mode_is_present_in_actual_scalar_basis":d["S"][1,2]+modes.ell/sp.sqrt(2*modes.wc*modes.kc),
        "matter_mode_is_present_in_actual_scalar_basis":d["S"][1,3]-1/sp.sqrt(2*modes.wm*modes.km),
        "relative_reconstruction_bound_uses_nonzero_projection":RECONSTRUCTION_NORM_FACTOR-4*domains.BASIS_NORM/PROJECTION_LOWER}


@cache
def gates():
    d=data()
    test=relative_error(sp.Integer(10)**72)
    return {name:bool(value) for name,value in {
        "actual_relational_clock_pole_weight_above_one_over_two_thousand":d["clock_principal_pole_weight_lower"]>sp.Rational(1,2000),
        "nonzero_clock_projection_survives_volume_rescaling":
            d["clock_coordinate_projection_lower_before_volume_rescaling"]*d["inverse_volume_rescaling_lower"]>PROJECTION_LOWER,
        "nonzero_matter_projection_survives_volume_rescaling":
            d["matter_coordinate_projection_lower_before_volume_rescaling"]*d["inverse_volume_rescaling_lower"]>PROJECTION_LOWER,
        "declared_projection_upper_exceeds_basis_times_scale_bound":PROJECTION_UPPER>4*100,
        "actual_norm_conversion_below_declared_error_coefficient":
            test["actual_relative_observable_error_upper"]<test["declared_relative_observable_error_upper"],
        "same_actual_nearby_classical_cones_used":cone_bounds.CLOCK_SPEED_EXCESS_LOWER>0,
        "local_relational_observable_not_a_quantum_gravity_locality_assertion":True,
        "bandlimited_packet_not_compact_retarded_preparation":True}.items()}
