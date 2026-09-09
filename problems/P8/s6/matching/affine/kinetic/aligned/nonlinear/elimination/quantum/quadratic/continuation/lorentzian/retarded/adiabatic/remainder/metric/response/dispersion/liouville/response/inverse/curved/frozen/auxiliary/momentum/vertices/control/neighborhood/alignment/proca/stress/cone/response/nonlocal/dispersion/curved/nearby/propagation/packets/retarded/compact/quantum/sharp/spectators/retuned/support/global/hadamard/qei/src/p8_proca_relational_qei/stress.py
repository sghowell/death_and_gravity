"""Exact divergence of the fixed-background test tensor and its actual source defect."""
from functools import cache

import sympy as sp
from p8_proca_global_support import charts, model

from . import modes

clean=modes.clean


@cache
def data():
    metric=sp.diag(-1,1,1,1)
    grad=sp.Matrix(sp.symbols("scalar_gradient_0:4",real=True))
    hessian=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"scalar_Hessian_{min(i,j)}_{max(i,j)}",real=True))
    norm=(grad.T*metric*grad)[0]
    tensor=grad*grad.T-metric*norm/2
    divergence=sp.Matrix([sum(metric[mu,mu]*sum(
        sp.diff(tensor[mu,nu],grad[r])*hessian[r,mu] for r in range(4)) for mu in range(4))
        for nu in range(4)])
    box=sum(metric[i,i]*hessian[i,i] for i in range(4))
    bg=model.background()
    force=clean(charts.data()["gamma"]["complete_leading_spatial_force"].subs(model.substitution(),simultaneous=True))
    b,chi,q=sp.symbols("clock_gamma_coordinate relational_matter_coordinate physical_spatial_momentum_squared",real=True)
    scalar_box=q*((force-sp.eye(2))*sp.Matrix([b,chi]))[1]
    expected=q*bg["ell"]*(1-bg["clock_physical_squared_speed"])*b
    return {"local_orthonormal_metric":metric,"scalar_gradient":grad,"symmetric_scalar_Hessian":hessian,
        "fixed_background_test_stress_tensor":tensor,"test_energy_sum_of_squares":tensor[0,0],
        "actual_covariant_divergence_in_normal_coordinates":divergence,
        "scalar_dalembertian":box,
        "actual_complete_principal_scalar_force":force,
        "actual_relational_scalar_principal_dalembertian":clean(scalar_box),
        "actual_clock_mode_principal_test_tensor_source":expected,
        "bounce_principal_box_source_coefficient_per_q_b":sp.Rational(8,6075),
        "positive_clock_physical_cotangent_metric_norm_per_q":1-bg["clock_physical_squared_speed"],
        "test_tensor_is_not_the_full_covariant_matter_metric_clock_stress":True,
        "no_interacting_Ward_identity_or_semiclassical_source_is_assigned":True}


@cache
def checks():
    d=data()
    bg=model.background()
    c2=bg["clock_physical_squared_speed"]
    expected=sp.Matrix([[c2,0],[bg["ell"]*(1-c2),1]])
    return {
        "test_tensor_divergence_is_box_phi_times_gradient":clean(
            d["actual_covariant_divergence_in_normal_coordinates"]-d["scalar_dalembertian"]*d["scalar_gradient"]),
        "test_energy_density_is_positive_sum_of_four_derivative_squares":sp.factor(
            d["test_energy_sum_of_squares"]-sum(v*v for v in d["scalar_gradient"])/2),
        "actual_coupled_clock_matter_principal_force_rebuilt":clean(d["actual_complete_principal_scalar_force"]-expected),
        "actual_matter_dalembertian_has_nonzero_clock_principal_source":clean(
            d["actual_relational_scalar_principal_dalembertian"]-d["actual_clock_mode_principal_test_tensor_source"]),
        "bounce_matter_box_source_coefficient_is_eight_over_6075":sp.factor(
            (bg["ell"]*(1-c2)).subs(model.u,0)-d["bounce_principal_box_source_coefficient_per_q_b"]),
        "bounce_clock_covector_is_physical_metric_spacelike":sp.factor(
            d["positive_clock_physical_cotangent_metric_norm_per_q"].subs(model.u,0)-sp.Rational(16,1215)),
    }


@cache
def gates():
    d=data()
    return {"actual_clock_principal_stress_source_is_not_zero":bool(d["bounce_principal_box_source_coefficient_per_q_b"]>0),
        "actual_clock_wavefront_is_not_on_single_matter_metric_null_cone":bool(sp.Rational(16,1215)>0)}
