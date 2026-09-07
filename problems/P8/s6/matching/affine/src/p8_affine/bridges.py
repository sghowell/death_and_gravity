"""Discharge literal-action, original-target and independent-domain interfaces."""
from functools import cache

import sympy as sp

from . import audit, connection, lower
from . import dictionary as d


@cache
def identities():
    values = d.lift()
    q, qphi = sp.symbols("q qphi", real=True)
    cubic = -(q+values["f_phi"])/(4*values["p"]*d.x)
    substitutions = {
        connection.P: values["p"], connection.S: sp.sqrt(-d.x),
        connection.PP: values["p_phi"], connection.PX: values["px"],
        connection.CP: sp.diff(values["c"], d.u), connection.CX: values["cx"],
        connection.F3: cubic,
    }
    coefficients = connection.coefficients()
    actual = {key: sp.simplify(coefficients[key].subs(substitutions))
              for key in ("P", "Q1", "Q2", "A1", "A2", "A3", "A4", "A5")}
    expected = d.target()
    result = {f"literal_original_alpha{i}": sp.simplify(actual[f"A{i}"]-expected[f"alpha{i}"])
              for i in range(1, 6)}
    curvature = connection.source_formula_comparison()["curvature_coefficient"].subs(substitutions)
    scalar = expected["scalar_F"]+d.x*qphi-3*d.x*(q+2*values["f_phi"])**2/(16*values["p"]**2)
    equation = lower.ode()
    result.update({
        "literal_original_curvature": sp.simplify(curvature-expected["f"]),
        "literal_lower_Q1": sp.simplify(actual["Q1"]-q),
        "literal_lower_Q2_ODE": sp.simplify(actual["Q2"]-2*(equation["forcing"]-equation["coefficient"]*q)),
        "literal_full_original_scalar_after_boundary": sp.simplify(actual["P"]+scalar-d.x*qphi-expected["scalar_F"]),
        "same_pointwise_c": sp.simplify(connection.cd_substitution()[connection.C].subs(substitutions)-values["c"]),
        "same_pointwise_quartic": sp.simplify(connection.cd_substitution()[connection.F4].subs(substitutions)-values["quartic"]),
    })
    independent = audit.parameters()
    independent_map = {audit.X: d.x, audit.U: d.u, audit.H: values["h"]}
    for key, other in (("p", "p"), ("f", "f"), ("c", "c"), ("quartic", "F4"), ("px", "p_x"), ("cx", "c_x")):
        result[f"independent_dictionary_{key}"] = sp.simplify(values[key]-independent[other].subs(independent_map))
    independent_lower = audit.lower_order()
    for key, other in (("coefficient", "ode_coefficient"), ("forcing", "ode_forcing")):
        result[f"independent_ODE_{key}"] = sp.simplify(equation[key]-independent_lower[other].subs(independent_map))
    for key in ("f", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5"):
        result[f"independent_original_target_{key}"] = sp.simplify(expected[key]-audit.target()["paper"][key].subs(independent_map))
    return result


@cache
def checks():
    bounds, independent = connection.calibration(), audit.calibration()
    embedding = connection.quotient()["embedding"]
    embed_norm = max(sum(abs(value) for value in embedding.row(i)) for i in range(64))
    transpose_norm = max(sum(abs(value) for value in embedding.col(i)) for i in range(60))
    result = {
        "closed_lower_p_squared": bounds["CD_p_squared_lower"] == independent["p_squared_min"],
        "closed_upper_p_squared": bounds["CD_p_squared_upper"] == independent["p_squared_max"],
        "full_rank_extra_factor_discharged": bounds["CD_exceptional_factor_lower"] == independent["quotient_factor_min"],
        "extra_factor_strictly_positive": bounds["CD_exceptional_factor_lower"] > 0,
        "inverse_parameter_lower_valid": sp.Rational(9, 20)**2 < bounds["CD_p_squared_lower"],
        "inverse_parameter_upper_valid": bounds["CD_p_squared_upper"] < sp.Rational(11, 20)**2,
        "inverse_row_sum_below_25": bounds["inverse_norm_upper"] < 25,
        "trace_gauge_embedding_norm": embed_norm == 3,
        "trace_gauge_source_projection_norm": transpose_norm == 2,
        "full_gauge_fixed_solution_bound": embed_norm*transpose_norm*bounds["inverse_norm_upper"] < 150,
        "source_palatini_correction_discharged": connection.palatini_control()["Q2_residual"] == 0,
        "printed_formula_is_not_silently_used": connection.source_formula_comparison()["printed_Q2_defect"] != 0,
        "excluded_locus_remains_a_rank_control": bounds["exceptional_quotient_rank"] == 57,
        "algebraic_not_a_propagating_health_claim": bounds["algebraic_only"] is True and bounds["UV_or_kinetic_health_claim"] is False,
    }
    return {name: bool(value) for name, value in result.items()}
