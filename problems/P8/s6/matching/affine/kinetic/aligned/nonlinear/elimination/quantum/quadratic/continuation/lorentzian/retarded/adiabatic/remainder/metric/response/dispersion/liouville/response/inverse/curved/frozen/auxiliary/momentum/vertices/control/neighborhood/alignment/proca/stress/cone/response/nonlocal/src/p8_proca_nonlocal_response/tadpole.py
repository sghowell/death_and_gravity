"""Response of the already-fixed new profile and complete clock-chart pullback."""
from functools import cache

import sympy as sp
from p8_proca_local_response import chart as local_chart
from p8_proca_stress import estimates, profiles
from p8_vector_state import wkb

from . import evolution


@cache
def physical_vertices():
    N,V=sp.symbols("ordinary_profile_physical_lapse ordinary_profile_logscale",real=True)
    rho,pressure=sp.symbols("fixed_new_ordinary_energy fixed_new_ordinary_pressure",real=True)
    chosen=profiles.action()
    P=chosen["new_clock_tube_profile"].subs({
        chosen["new_fixed_energy"]:rho,chosen["same_fixed_pressure"]:pressure,
        profiles.x:-N**-2},simultaneous=True)
    density=N*sp.exp(3*V)*P
    point={N:1,V:0}
    gradient=sp.Matrix([sp.diff(density,field).subs(point) for field in (N,V)])
    hessian=sp.Matrix([[sp.diff(density,left,right).subs(point) for right in (N,V)] for left in (N,V)])
    energy=-sp.exp(-3*V)*sp.diff(density,N)
    spatial=sp.exp(-3*V)*sp.diff(density,V)/(3*N)
    jacobian=sp.Matrix([[sp.diff(value,field).subs(point) for field in (N,V)] for value in (energy,spatial)])
    expected_hessian=sp.Matrix([[-rho-pressure,3*rho],[3*rho,-9*pressure]])
    expected_jacobian=sp.Matrix([[rho+pressure,0],[rho+pressure,0]])
    w1,w2=local_chart.mapping()["omega_N"],local_chart.mapping()["omega_NN"]
    J=sp.Matrix([[1,0],[w1,1]])
    actual_clock_output=J.T*sp.diag(-1,3)
    target_clock_output=sp.Matrix([[-1,3*w1],[0,3]])
    gaussian_gradient=sp.Matrix([-rho,3*pressure])
    return {"rho":rho,"pressure":pressure,"fixed_profile_density_gradient":sp.ImmutableMatrix(gradient),
            "fixed_profile_density_hessian":sp.ImmutableMatrix(hessian),
            "fixed_profile_physical_stress_jacobian":sp.ImmutableMatrix(jacobian),
            "linear_metric_chart_Jacobian":sp.ImmutableMatrix(J),
            "background_cancelled_clock_output_matrix":sp.ImmutableMatrix(actual_clock_output),
            "checks":{"new_fixed_profile_gradient":sp.ImmutableMatrix(gradient+gaussian_gradient),
                "new_fixed_profile_Hessian":sp.ImmutableMatrix(hessian-expected_hessian),
                "new_fixed_profile_physical_normalization":sp.ImmutableMatrix(jacobian-expected_jacobian),
                "total_second_metric_map_contact_cancels_only_after_fixed_profile":sp.factor(
                    (gaussian_gradient[1]+gradient[1])*w2),
                "actual_background_cancelled_clock_current_map":sp.ImmutableMatrix(actual_clock_output-target_clock_output)},
            "policy":"S6.82's new ordinary-Proca c-number profile is fixed, including its cutoff plateau. No state or profile is reevaluated. The nonlinear chart's extra Hessian contact cancels for the total background-cancelled Gaussian-plus-profile action, not for its metric piece alone."}


@cache
def chart_lift():
    w1=local_chart.mapping()["omega_N"]
    derivatives={j:wkb.box_bound(sp.factor(sp.diff(w1,wkb.u,j))) for j in range(11)}
    rows={j:1+sum(sp.binomial(j,k)*derivatives[j-k]["absolute_upper"] for k in range(j+1))
          for j in range(11)}
    return {"coefficient_derivative_envelopes":derivatives,"Leibniz_C10_row_bounds":rows,
            "joint_input_C10_upper":max(sp.Integer(1),*rows.values()),
            "joint_output_C0_upper":sp.Integer(3),
            "output_bound_uses_omega_N_at_most_one_half":bool(derivatives[0]["absolute_upper"]<=sp.Rational(1,2))}


def bound(scale):
    return _bound(estimates.exact_scale(scale))


@cache
def _bound(scale):
    vector=evolution.bound(scale)
    state=estimates.physical_bounds(scale)
    values={"energy":state["new_energy_derivative_bounds"][0]["total"],
            "pressure":state["identical_pressure_derivative_bounds"][0]["total"]}
    additional=sum(values.values())
    total={key:value+additional for key,value in vector["complete_metric_response_C10_to_C0_component_bounds"].items()}
    lift=chart_lift()
    clock=lift["joint_input_C10_upper"]*lift["joint_output_C0_upper"]*max(total.values())
    return {"M_tau":scale,"fixed_mass_time_product":estimates.MASS,
            "new_fixed_profile_C0_upper_bounds":values,
            "new_fixed_profile_physical_response_C0_to_C0_upper":additional,
            "background_cancelled_physical_C10_to_C0_component_bounds":total,
            "joint_background_cancelled_physical_C10_to_C0_upper":max(total.values()),
            "clock_chart_input_C10_upper":lift["joint_input_C10_upper"],
            "joint_background_cancelled_clock_C10_to_C0_upper":clock,
            "no_no_loss_or_coupled_inverse_inference":True}
