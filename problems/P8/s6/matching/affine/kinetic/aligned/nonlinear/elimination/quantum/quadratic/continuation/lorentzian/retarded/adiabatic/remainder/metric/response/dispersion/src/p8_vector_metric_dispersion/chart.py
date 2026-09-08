"""Regular metric chart and cancellation of its nonlinear Hessian contact."""
from functools import cache

import sympy as sp
from p8_clock_tadpole import profiles
from p8_margin_response import system
from p8_vector_mass_adiabatic.bounds import rational_bound
from p8_vector_metric_response import tadpole

from . import principal, spectral

u = profiles.u
h = spectral.h


def linear():
    return sp.ImmutableMatrix([[1, 0], [1/(2*h), 1]])


@cache
def data():
    N = sp.Symbol("chart_lapse", real=True)
    omega = -sp.log((h-1+N**-2)/h)/4
    delta = sp.diff(omega, N).subs(N, 1)
    second = sp.diff(omega, N, 2).subs(N, 1)
    return {"omega": omega, "first_lapse_jet": delta, "second_lapse_jet": second,
            "linear": linear(),
            "pole_principal": sp.ImmutableMatrix((linear().T*principal.pole()*linear()).applyfunc(sp.factor)),
            "finite_principal": sp.ImmutableMatrix((linear().T*principal.finite()*linear()).applyfunc(sp.factor)),
            "source_order": ["regular_lapse", "regular_logscale"],
            "physical_source_order": ["physical_lapse", "physical_logscale"]}


@cache
def checks():
    out = {}
    N, v = sp.symbols("chart_lapse chart_scale", real=True)
    qn, qz, knn, knz, kzz = sp.symbols("gradient_N gradient_Z H_NN H_NZ H_ZZ", real=True)
    n, zeta = sp.symbols("physical_n physical_zeta", real=True)
    density = qn*n+qz*zeta+(knn*n**2+2*knz*n*zeta+kzz*zeta**2)/2
    changed = density.subs({n: N-1, zeta: v+data()["omega"]}, simultaneous=True)
    actual = sp.ImmutableMatrix([[sp.diff(changed, left, right).subs({N: 1, v: 0})
                                 for right in (N, v)] for left in (N, v)])
    expected = linear().T*sp.ImmutableMatrix([[knn, knz], [knz, kzz]])*linear()
    extra = sp.diag(qz*data()["second_lapse_jet"], 0)
    out["nonlinear_chart_Hessian_chain_rule"] = sp.ImmutableMatrix((actual-expected-extra).applyfunc(sp.factor))
    out["nonlinear_chart_first_jet"] = sp.factor(data()["first_lapse_jet"]-1/(2*h))
    out["nonlinear_chart_second_jet"] = sp.factor(data()["second_lapse_jet"]+3/(2*h)-1/h**2)
    vertices = tadpole.physical_vertices()
    rho, pressure = vertices["rho"], vertices["pressure"]
    vector_gradient = sp.ImmutableMatrix([-rho, 3*pressure])
    out["total_quantum_background_gradient_vanishes"] = sp.ImmutableMatrix(
        vector_gradient+vertices["density_gradient"])
    out["chart_contacts_cancel_only_after_fixed_tadpole"] = sp.factor(
        (vector_gradient[1]+vertices["density_gradient"][1])*data()["second_lapse_jet"])
    drho, dp = sp.symbols("varied_total_energy varied_total_pressure", real=True)
    force = linear().T*sp.ImmutableMatrix([-drho, 3*dp])
    old_force = system.source()
    mapping = {system.rho: drho, system.pressure: dp, system.delta: 1/(2*h)}
    out["regular_stress_force_reproduces_frozen_tree_source"] = sp.ImmutableMatrix(
        force-sp.ImmutableMatrix([old_force["lapse_force"], old_force["curvature_force"]]).subs(mapping))
    out["regular_finite_fourth_lapse_coefficient"] = sp.factor(
        data()["finite_principal"][0, 0]+sp.Rational(7357, 6561)/h**2)
    return out


@cache
def norm():
    delta = 1/(2*(1+u**2)**3)
    derivatives = tuple(rational_bound(sp.diff(delta, u, j)) for j in range(11))
    by_order = tuple(1+sum(sp.binomial(k, j)*derivatives[j] for j in range(k+1)) for k in range(11))
    coefficient = max(by_order)
    kappa = tadpole.bound(10**24, 1000)["joint_background_cancelled_metric_C10_to_C0_upper_bound"]
    return {"chart_derivative_upper": derivatives, "chart_C10_product_upper": coefficient,
            "source_stress_C10_to_C0_upper": coefficient*kappa,
            "one_retarded_tree_application_C10_to_C0_upper": {
                "regular_phase": 131*coefficient*kappa,
                "regular_lapse": 94*coefficient*kappa,
                "matter": 44*coefficient*kappa,
                "physical_logscale": 225*coefficient*kappa},
            "scope": "One zero-initial classical forced response to the background-cancelled prepared Gaussian stress. No same-norm feedback iteration, lapse derivative estimate, resummed solution or residual error bound."}
