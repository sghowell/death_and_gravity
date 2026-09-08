"""Response of the already-fixed scalar tadpole action in the physical metric."""
from functools import cache

import sympy as sp
from p8_clock_tadpole import profiles
from p8_vector_hadamard import transfer

from . import evolution


@cache
def physical_vertices():
    N, zeta = sp.symbols("physical_lapse physical_logscale", real=True)
    rho, pressure = sp.symbols("fixed_selected_energy fixed_selected_pressure", real=True)
    original = profiles.covariant()
    P = original["P_on_flat_clock_tube"].subs({
        original["rho"]: rho, original["pressure"]: pressure, profiles.x: -N**-2}, simultaneous=True)
    density = N*sp.exp(3*zeta)*P
    point = {N: 1, zeta: 0}
    gradient = sp.ImmutableMatrix([sp.diff(density, field).subs(point) for field in (N, zeta)])
    hessian = sp.ImmutableMatrix([[sp.diff(density, left, right).subs(point) for right in (N, zeta)] for left in (N, zeta)])
    target_gradient = sp.ImmutableMatrix([rho, -3*pressure])
    target_hessian = sp.ImmutableMatrix([[-rho-pressure, 3*rho], [3*rho, -9*pressure]])
    physical_energy = -sp.exp(-3*zeta)*sp.diff(density, N)
    physical_pressure = sp.exp(-3*zeta)*sp.diff(density, zeta)/(3*N)
    physical_jacobian = sp.ImmutableMatrix([[sp.diff(readout, field).subs(point) for field in (N, zeta)]
                                           for readout in (physical_energy, physical_pressure)])
    target_jacobian = sp.ImmutableMatrix([[rho+pressure, 0], [rho+pressure, 0]])
    return {"rho": rho, "pressure": pressure, "density_gradient": gradient, "density_hessian": hessian,
            "physical_stress_jacobian": physical_jacobian,
            "checks": {"existing_tadpole_physical_gradient": sp.ImmutableMatrix(gradient-target_gradient),
                       "existing_tadpole_physical_hessian": sp.ImmutableMatrix(hessian-target_hessian),
                       "existing_tadpole_physical_output_jacobian": sp.ImmutableMatrix(physical_jacobian-target_jacobian)},
            "policy": "The S6.60 selected c-number profiles and its scalar-clock cutoff are held fixed. The cutoff is identically one near the clock; no state or profile is reselected under either metric variation."}


def bound(planck_time_product, mass_time_product):
    vector = evolution.bound(planck_time_product, mass_time_product)
    state = transfer.physical_bounds(planck_time_product, mass_time_product)
    profiles_upper = {component: state[component+"_new_state_total_over_reference_density"] for component in ("energy", "pressure")}
    additional = sum(profiles_upper.values())
    complete = {component: value+additional for component, value in vector["complete_metric_response_C10_to_C0_component_bounds"].items()}
    return {"fixed_profile_C0_upper_bounds": profiles_upper,
            "existing_tadpole_physical_C0_to_C0_upper_bound": additional,
            "background_cancelled_metric_C10_to_C0_component_bounds": complete,
            "joint_background_cancelled_metric_C10_to_C0_upper_bound": max(complete.values()),
            "scope": "The Gaussian vector plus the already-fixed S6.60 scalar tadpole, for prepared homogeneous metric sources. No other loops, interacting feedback or quantum health verdict."}
