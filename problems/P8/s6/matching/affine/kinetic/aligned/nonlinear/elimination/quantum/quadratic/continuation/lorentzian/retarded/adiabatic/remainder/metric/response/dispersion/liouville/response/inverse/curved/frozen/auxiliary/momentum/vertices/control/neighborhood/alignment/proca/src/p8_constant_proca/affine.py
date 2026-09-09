"""Literal retained-mass update of the full source-free affine action."""
from functools import cache

import sympy as sp
from p8_affine import connection
from p8_affine_retuned import geometry


@cache
def data():
    d=geometry.update()
    source=geometry.source_centering()
    split=geometry.complement()
    L,E=d["N_full"],d["embedding"]
    oldmass=d["D"].inv()
    newmass=sp.diag(1,-1,-1,-1)
    change=geometry.clean(newmass-oldmass)
    B=sp.Symbol("fixed_original_B_phi_x",real=True)
    shift=sp.ImmutableMatrix([connection.S*B,0,0,0])
    lift=geometry.clean(E*split["lift"])
    center=geometry.clean(source["stationary"]-lift*(source["Tstar"]-shift))
    old_source=source["source_new"]+L.T*oldmass*(source["Tstar"]-shift)
    added_Hessian=geometry.clean(L.T*change*L)
    added_source=geometry.clean(-L.T*change*shift)
    new_Hessian=geometry.clean(source["M_full"]+added_Hessian)
    new_source=geometry.clean(old_source+added_source)
    trace=sp.ImmutableMatrix(sp.symbols("unshifted_retained_trace0:4",real=True))
    W=trace-shift
    counterterm=geometry.clean((W.T*change*W)[0]/2)
    return {"original_one_form_shift":shift,"source_free_stationary_connection":center,
            "new_full_connection_Hessian":new_Hessian,"new_full_connection_source":new_source,
            "added_full_connection_Hessian":added_Hessian,"added_full_connection_source":added_source,
            "retained_mass_change":change,"new_retained_mass":newmass,
            "old_retained_mass":oldmass,"literal_added_affine_density":counterterm,
            "unshifted_retained_trace_variables":trace,"retained_right_inverse":lift,
            "new_full_inverse_trace_lift":geometry.clean(lift*newmass),
            "quotient_determinant_ratio":sp.factor((sp.eye(4)+change*d["D"]).det())}


@cache
def checks():
    new=data()
    old=geometry.update()
    split=geometry.complement()
    L,E=old["N_full"],old["embedding"]
    R=new["retained_right_inverse"]
    trace=new["unshifted_retained_trace_variables"]
    counterterm=new["literal_added_affine_density"]
    gradient=sp.Matrix([sp.diff(counterterm,x) for x in trace])
    return {"literal_mass_counterterm_Hessian":geometry.clean(sp.hessian(counterterm,trace)-new["retained_mass_change"]),
            "literal_mass_counterterm_source":geometry.clean(L.T*gradient.subs(dict.fromkeys(trace,0))-new["added_full_connection_source"]),
            "new_all64_stationary_Euler_equations":geometry.clean(
                new["new_full_connection_Hessian"]*new["source_free_stationary_connection"]+new["new_full_connection_source"]),
            "new_stationary_trace_is_same_original_B_dphi":geometry.clean(
                L*new["source_free_stationary_connection"]-new["original_one_form_shift"]),
            "new_full_inverse_trace_lift_Euler_equations":geometry.clean(
                new["new_full_connection_Hessian"]*new["new_full_inverse_trace_lift"]-L.T),
            "new_retained_response_exact_constant_Lorentz_matrix":geometry.clean(
                L*new["new_full_inverse_trace_lift"]-new["new_retained_mass"]),
            "new_Hessian_preserves_all_projective_columns":geometry.clean(
                new["new_full_connection_Hessian"]*connection.quadratic()["gauge"]),
            "new_source_preserves_projective_annihilation":geometry.clean(
                connection.quadratic()["gauge"].T*new["new_full_connection_source"]),
            "retained_map_still_annihilates_all_56_complement_directions":geometry.clean(L*E*split["projector"]),
            "new_mass_has_no_complement_retained_cross_block":geometry.clean(
                split["projector"].T*E.T*new["new_full_connection_Hessian"]*R),
            "new_retained_mass_from_full_connection_Hessian":geometry.clean(
                R.T*new["new_full_connection_Hessian"]*R-new["new_retained_mass"]),
            "exact_quotient_determinant_lemma_ratio":sp.factor(
                new["quotient_determinant_ratio"]-old["gamma_t"]*old["gamma_s"]**3)}
