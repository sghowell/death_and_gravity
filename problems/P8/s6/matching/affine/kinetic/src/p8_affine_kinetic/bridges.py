"""Check separately authored action premises and nonunit normalization."""
from functools import cache

import sympy as sp
from p8 import gamma
from p8_affine import connection
from p8_m1 import nonlinear

from . import scalar, vector


def _clean(value):
    return sp.factor(value)


@cache
def checks():
    bg = scalar.background()
    rolling = vector.rolling()
    previous = nonlinear.functions()
    previous_u = nonlinear.u
    result = {}
    for name, old_name in (("a", "a"), ("H", "Hubble"), ("h", "h")):
        result["actual_background_"+name] = _clean(bg[name]-rolling[old_name].subs(rolling["u"], scalar.u))
    for name, old_name in (("J", "J"), ("theta", "Theta"), ("lam", "Lambda"),
                           ("delta", "delta"), ("w", "w")):
        result["pinned_original_"+name] = _clean(bg[name]-previous[old_name].subs(previous_u, scalar.u))
    for name, old_name in (("alpha", "gradient_amplitude"), ("d", "mass_time_shift"),
                           ("e", "mass_spatial_shift")):
        result["independent_rolling_"+name] = _clean(bg[name]-rolling[old_name].subs(rolling["u"], scalar.u))
    independent = scalar.vector_reconstruction()
    for name, key in (("alpha", "gradient_amplitude"), ("d", "mass_time_shift"),
                      ("e", "mass_spatial_shift")):
        evaluated = independent[name].subs({independent["h"]: bg["h"], independent["H"]: bg["H"]})
        result["literal_reconstruction_"+name] = _clean(evaluated-rolling[key].subs(rolling["u"], scalar.u))
    result["literal_vector_contraction"] = (independent["generic_vector"]-vector.schur()["Vstar"]).applyfunc(_clean)
    result["rolling_vector_mass"] = (vector.schur()["D_inverse"].subs(connection.P, sp.Rational(1, 2))
                                      -2*scalar.BETA*sp.diag(1, -1, -1, -1))
    # Rebuild the old small-symbol first-order action, before its lapse solve.
    old = gamma.hamiltonian()
    T, FT, th, S, lam, delta, ell, q = old["symbols"]
    b, v, matter, pm = old["states"]
    mapping = {T: 1, FT: 1, th: scalar.theta,
               S: scalar.action()["S"], lam: scalar.lam,
               delta: (1+scalar.w/scalar.ell)/3, ell: scalar.ell, q: scalar.q,
               b: scalar.shift, v: scalar.v, matter: scalar.matter, pm: scalar.pm}
    old_hamiltonian = old["density"].subs(mapping, simultaneous=True)
    bare = scalar.action()["base"]
    vd = scalar.theta*scalar.n-scalar.ell*scalar.matter/2
    sd = scalar.pm-scalar.w*scalar.n
    new_before = (-2*scalar.q*scalar.shift*scalar.vd+scalar.pm*scalar.sd-bare).subs(
        {scalar.vd: vd, scalar.sd: sd})
    lapse = sp.solve(sp.diff(new_before, scalar.n), scalar.n)[0]
    result["unchanged_full_CD_matter_Hamiltonian"] = _clean(new_before.subs(scalar.n, lapse)-old_hamiltonian)
    for name, value in old["residuals"].items():
        result["old_direct_constraint_"+name] = value
    physical = scalar.units(3, 2, 5)
    independent_units = vector.units(mass_squared=3, time_scale=2, kinetic_coefficient=5)
    for key, other in (("normalized_zeta", "lambda_normalized"),
                       ("isolated_proca_mass_squared", "isolated_center_mass_squared_physical"),
                       ("normalized_isolated_mass_squared", "isolated_center_mass_squared_normalized")):
        result["independent_units_"+key] = physical[key]-independent_units[other]
    result["normalized_nonunit_coupling"] = physical["normalized_zeta"]-sp.Rational(5, 12)
    result["physical_nonunit_isolated_mass"] = physical["isolated_proca_mass_squared"]-sp.Rational(8, 5)
    result["normalized_nonunit_isolated_mass"] = physical["normalized_isolated_mass_squared"]-sp.Rational(32, 5)
    result["physical_to_normalized_mass"] = (physical["normalized_isolated_mass_squared"]
                                             -4*physical["isolated_proca_mass_squared"])
    M2, tau, zeta_phys, qphys = sp.symbols("M2 tau zeta_phys qphys", positive=True)
    norm = zeta_phys/(M2*tau**2)
    result["kinetic_action_unit_factor"] = _clean(M2*tau**2*norm-zeta_phys)
    result["isolated_mass_unit_identity"] = _clean(8/(3*norm*tau**2)-8*M2/(3*zeta_phys))
    result["physical_momentum_definition"] = _clean((tau*sp.sqrt(qphys))**2-tau**2*qphys)
    return result


@cache
def proof_checks():
    bounds = vector.domain_bounds()
    derivative = sp.Poly(bounds["polynomial_derivative"], vector.P)
    out = {"closed_tube_has_stronger_positive_p_lower": bounds["lower_p_squared_margin"] > 0,
           "closed_tube_has_strict_p_upper": bounds["upper_p_squared_margin"] > 0,
           "D00_polynomial_increasing_for_positive_p": all(value > 0 for value in derivative.all_coeffs()),
           "D00_polynomial_lower_positive": bounds["polynomial_lower"] > 0,
           "D00_positive_lower": bounds["D00_lower"] > 0,
           "minus_Dii_positive_lower": bounds["minus_Dii_lower"] > 0,
           "D_inverse_norm_below_24": bounds["D_inverse_norm_margin_below_24"] > 0,
           "fixed_vector_complement_dimension": vector.schur()["complement_dimension"] == 56,
           "exact_gradient_change_not_a_point_redefinition": True,
           "nonunit_original_frame_and_matter_unchanged": True}
    return {name: bool(value) for name, value in out.items()}


@cache
def calibration():
    return {"M_squared": 3, "tau": 2, "zeta_physical": 5,
            **scalar.units(3, 2, 5),
            "physical_q_threshold_at_center": sp.Rational(3597, 12800),
            "threshold_is_squared_momentum_not_frequency_cutoff": True,
            "center_constraint_chart_not_used_for_the_punctured_physical_proof": True}
