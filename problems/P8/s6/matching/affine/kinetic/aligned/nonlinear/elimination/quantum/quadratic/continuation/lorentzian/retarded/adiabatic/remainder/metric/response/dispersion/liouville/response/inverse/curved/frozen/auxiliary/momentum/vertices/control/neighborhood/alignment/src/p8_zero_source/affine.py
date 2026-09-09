"""Literal local affine lift of the normal retained-source counterterm.

All light fields, including their derivatives in S, are fixed in the
connection Hessian. D has the original (+---) mass-square convention.
The physical matter metric itself is not changed.
"""
from functools import cache

import sympy as sp
from p8_affine import connection
from p8_affine_retuned import degeneracy, geometry


@cache
def data():
    d=geometry.update()
    old=geometry.source_centering()
    split=geometry.complement()
    L,E,W,D=d["N_full"],d["embedding"],d["W"],d["D"]
    sigma=sp.Symbol("actual_normal_light_source",real=True)
    S=sp.ImmutableMatrix([sigma,0,0,0])
    added_source=geometry.clean(L.T*D.inv()*S)
    stationary=geometry.clean(old["stationary"]-E*W*D.inv()*S)
    # This is the existing shift in W = L C - B dphi, expressed via
    # the original retained center and its normal residual source.
    shift=old["Tstar"]-S
    Ctrace=sp.ImmutableMatrix(sp.symbols("affine_retained_trace0:4",real=True))
    change=geometry.expanded_rational(((Ctrace-shift).T*D.inv()*S)[0]
                                     -(S.T*D.inv()*S)[0]/2)
    return {"normal_source_vector":S,"retained_trace_variables":Ctrace,
            "literal_counterterm_at_fixed_lights":change,
            "added_connection_source":added_source,
            "new_stationary_connection":stationary,"old_one_form_shift":shift,
            "new_connection_source":old["source_new"]+added_source,
            "unchanged_full_connection_Hessian":old["M_full"],
            "unchanged_retained_mass_matrix":D.inv(),
            "unchanged_56_complement_projector":split["projector"]}


@cache
def checks():
    d=data()
    old=geometry.source_centering()
    g=geometry.update()
    S=d["normal_source_vector"]
    trace=d["retained_trace_variables"]
    change=d["literal_counterterm_at_fixed_lights"]
    gradient=sp.ImmutableMatrix([sp.diff(change,x) for x in trace])
    normal=trace[0]-d["old_one_form_shift"][0]
    expected=(normal*S[0]-S[0]**2/2)/g["gamma_t"]
    out={
        "literal_normal_counterterm_matches_retained_mass_convention":geometry.expanded_rational(change-expected),
        "retained_gradient_lifts_to_all64_connection_sources":geometry.clean(g["N_full"].T*gradient-d["added_connection_source"]),
        "connection_Hessian_change_is_identically_zero":sp.hessian(change,trace),
        "all64_new_stationary_Euler_equations":geometry.clean(
            old["M_full"]*d["new_stationary_connection"]+d["new_connection_source"]),
        "new_stationary_shifted_one_form_has_zero_normal_source":geometry.clean(
            g["N_full"]*d["new_stationary_connection"]-d["old_one_form_shift"]),
        "projective_source_annihilation":geometry.clean(
            connection.quadratic()["gauge"].T*d["added_connection_source"]),
        "complement_source_annihilation":geometry.clean(
            d["unchanged_56_complement_projector"].T*g["N"].T*g["D"].inv()*S),
        "old_stationary_recovered_when_source_zero":geometry.clean(
            d["new_stationary_connection"].subs(S[0],0)-old["stationary"]),
        "retained_mass_response_signs_unchanged":geometry.clean(g["D"]-g["expected_D"]),
        "actual_original_source_is_covariantly_normal":degeneracy.source()["covariant_reconstruction"],
        "actual_original_source_has_no_spatial_component":degeneracy.source()["generic_component_residual"],
    }
    return out
