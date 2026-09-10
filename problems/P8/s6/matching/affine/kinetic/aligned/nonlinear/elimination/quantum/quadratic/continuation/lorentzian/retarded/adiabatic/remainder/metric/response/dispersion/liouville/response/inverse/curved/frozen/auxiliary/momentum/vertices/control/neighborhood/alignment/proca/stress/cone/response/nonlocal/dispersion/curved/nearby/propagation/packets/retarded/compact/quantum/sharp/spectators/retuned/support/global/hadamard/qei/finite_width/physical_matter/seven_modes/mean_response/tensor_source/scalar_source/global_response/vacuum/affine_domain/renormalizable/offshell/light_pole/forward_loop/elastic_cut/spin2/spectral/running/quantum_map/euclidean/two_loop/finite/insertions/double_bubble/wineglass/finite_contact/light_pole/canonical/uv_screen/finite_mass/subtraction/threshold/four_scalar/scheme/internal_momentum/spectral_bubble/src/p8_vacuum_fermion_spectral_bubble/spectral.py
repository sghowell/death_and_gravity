"""Direct spectral identity from the complete on-shell-subtracted kernel."""

from functools import cache

import sympy as sp
from p8_vacuum_global_one_loop_insertions import halfplane


@cache
def data():
    s, u, tau, T = sp.symbols("s u tau T")
    A, m = sp.symbols("A m", positive=True)
    H = (s - 1) ** 2 / ((u - 1) ** 2 * (u - s))
    antiderivative = -1 / (u - s) - (s - T) / (2 * (u - s) ** 2)
    integrated_second = -2 * (1 / (tau - s) + (s - T) / (2 * (tau - s) ** 2))
    original = -(T - s) * sp.log(1 - s / tau)
    density_integrand = (u - T) / ((u - 1) ** 2 * (u - s))
    old = halfplane.data()
    old_second = old["second_derivative_resolvent_decomposition"]
    v = old["symbols"]
    mapped = old_second.subs({v["m"]: m, v["A"]: A, v["s"]: s}, simultaneous=True)
    symbolic_f = sp.Function("renormalized_first_fermion_inverse")
    fR = symbolic_f(s) - symbolic_f(1) - (s - 1) * sp.diff(symbolic_f(s), s).subs(s, 1)
    h = sp.symbols("formal_loop_order")
    D0 = 1 / (1 - s)
    first_correction = sp.diff(1 / (1 - s + h * fR), h).subs(h, 0)
    checks = {
        "twice_differentiated_subtracted_Cauchy_kernel": sp.factor(
            sp.diff(H, s, 2) - 2 / (u - s) ** 3
        ),
        "spectral_mass_anchor": H.subs(s, 1),
        "spectral_residue_anchor": sp.diff(H, s).subs(s, 1),
        "elementary_spectral_integral_primitive": sp.factor(
            sp.diff(antiderivative, u) - (u - T) / (u - s) ** 3
        ),
        "integrated_second_derivative_matches_log_kernel": sp.factor(
            integrated_second - sp.diff(original, s, 2)
        ),
        "same_frozen_resolvent_derivative": sp.factor(
            integrated_second.subs({tau: m * m / A, T: 4 * m * m}, simultaneous=True)
            - mapped
        ),
        "parameter_spectral_density_gap": sp.factor(
            (tau - T).subs({tau: m * m / A, T: 4 * m * m}, simultaneous=True)
            - m * m * (1 - 4 * A) / A
        ),
        "formal_propagator_correction_sign": sp.factor(first_correction + fR * D0**2),
        "no_exact_geometric_resummation": sp.factor(
            sp.diff(1 / (1 - s + h * fR), h, 2).subs(h, 0) - 2 * fR**2 / (1 - s) ** 3
        ),
    }
    return {
        "symbols": {"s": s, "u": u, "tau": tau, "T": T},
        "per_parameter_spectral_integrand": density_integrand,
        "subtracted_Cauchy_kernel": H,
        "spectral_integral_second_derivative": integrated_second,
        "per_parameter_identity": "-g_x,R(s)/(s-1)^2=integral_tau^infinity (u-T)/[(u-1)^2(u-s)] du; tau=mF^2/A>=T=4mF^2.",
        "complete_identity": "-f_R(s)/(s-1)^2=integral_T^infinity w(u)/(u-s) du, w(u)=(2NY/Q) u(1-T/u)^(3/2)/(u-1)^2.",
        "formal_propagator_first_correction": first_correction,
        "domain": "Re(s)<T; the mass-one point is below threshold. Endpoint A=0 contributes zero, obtained by tau going to infinity.",
        "scope": "Positive spectral representation of one formal propagator correction, not the exact normalized propagator, full reflection positivity or an all-orders spectral theorem.",
        "checks": checks,
    }
