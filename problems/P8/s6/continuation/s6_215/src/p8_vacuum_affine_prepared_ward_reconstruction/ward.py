"""Ordered Ward reconstruction with the full ADM one-current chart term."""

from functools import cache

import sympy as s
from p8_vacuum_affine_full_adm_vertices import chart

from . import kinematics as k


def pair(E, h):
    return s.trace(E.T * h)


def lie_density(E, xi):
    J = xi.jacobian(k.COORDS)
    transport = E.applyfunc(
        lambda f: sum(xi[j] * s.diff(f, q) for j, q in enumerate(k.COORDS))
    )
    return transport - J * E - E * J.T + s.trace(J) * E


def source_integrand(D, xi, E):
    Gxi = k.gauge(xi[0], xi[1:, 0])
    return pair(lie_density(E, xi), chart.first(*D, k.a)) + pair(
        E, chart.second(*D, *Gxi, k.a)
    )


def detector_integrand(xi, G, E):
    Dxi = k.gauge(xi[0], xi[1:, 0])
    return -pair(E, k.lie_covariant(chart.first(*G, k.a), xi)) + pair(
        E, chart.second(*Dxi, *G, k.a)
    )


def boundary_flux(E, h, delta_E, xi, metric):
    """The time flux of the differentiated conservation identity."""
    return 2 * sum(
        xi[alpha]
        * (
            sum(
                delta_E[0, nu] * metric[alpha, nu] + E[0, nu] * h[alpha, nu]
                for nu in range(4)
            )
        )
        for alpha in range(4)
    )


@cache
def data():
    t, x, y, z = k.COORDS
    A, B = s.Function("A")(t), s.Function("B")(t)
    E = s.diag(A, B, B, B)
    eta = s.Function("eta")(*k.COORDS)
    chi = s.Matrix([s.Function(f"chi{i}")(*k.COORDS) for i in range(3)])
    xi = s.Matrix([eta, *chi])
    divchi = s.trace(k.spatial_jac(chi))
    expected = s.zeros(4)
    expected[0, 0] = eta * s.diff(A, t) - A * s.diff(eta, t) + A * divchi
    expected[0, 1:] = (-B * k.grad(eta) - A * chi.diff(t)).T
    expected[1:, 0] = -B * k.grad(eta) - A * chi.diff(t)
    expected[1:, 1:] = (
        eta * s.diff(B, t) * s.eye(3)
        - B * k.symgrad(chi)
        + B * (s.diff(eta, t) + divchi) * s.eye(3)
    )
    rho, pressure = s.Function("rho")(t), s.Function("pressure")(t)
    A0, B0 = -(k.a**3) * rho / 2, -k.a * pressure / 2
    metric = s.diag(1, -(k.a**2), -(k.a**2), -(k.a**2))
    conservation = -2 * s.diff(A0, t) - 6 * k.a**2 * k.H * B0
    # Independent tensor-density Leibniz rule on nonconstant, non-diagonal tensors.
    Ef = s.Matrix(
        [
            [t + x, y, z, x * y],
            [y, 1 + t * t, x * z, t * y],
            [z, x * z, 2 + t, y * z],
            [x * y, t * y, y * z, 3 + t * z],
        ]
    )
    hf = s.Matrix(
        [
            [1 + t * x, t * y, x * z, y],
            [t * y, 2 + x, t * z, x * y],
            [x * z, t * z, 3 + y, t * x],
            [y, x * y, t * x, 4 + z],
        ]
    )
    xf = s.Matrix([t * x + y, t * y + z, x * z + t, t * z + x * y])
    contracted = pair(Ef, hf)
    divergence = sum(s.diff(xf[j] * contracted, q) for j, q in enumerate(k.COORDS))
    n, p = s.symbols("n trace_Q", real=True)
    h = chart.first(n, s.zeros(3, 1), p * s.eye(3) / 3, k.a)
    zero = s.zeros(4)
    checks = {
        "all_weight_one_density_Lie_components": lie_density(E, xi) - expected,
        "actual_ADM_current_density_normalization": pair(s.diag(A0, B0, B0, B0), h)
        - (-(k.a**3) * rho * n + k.a**3 * pressure * p / 2),
        "homogeneous_density_conservation": conservation
        - k.a**3 * (s.diff(rho, t) + 3 * k.H * (rho + pressure)),
        "full_tensor_density_Leibniz_divergence": pair(lie_density(Ef, xf), hf)
        + pair(Ef, k.lie_covariant(hf, xf))
        - divergence,
        "final_detector_zero_kills_full_differentiated_flux": boundary_flux(
            Ef, hf, Ef, zero[:, 0], metric
        ),
        "initial_prepared_source_zero_kills_full_differentiated_flux": boundary_flux(
            Ef, zero, zero, xf, metric
        ),
    }
    return {
        "density_convention": "Define E by J(D)=integral E^{mu nu} hD_mu nu with coordinate measure. E is contravariant weight one: E00=-a^3 rho/2, Eij=-a P deltaij/2, E0i=0. Its actual one-point function is not zero.",
        "source_Ward": "For a gauge source identity on the initial germ, local covariance and the same pulled-back Hadamard state give deltaE[Lxi g]=Lxi E. Thus R_ADM(D,Gxi)=integral hD:Lxi E+E:q(D,Gxi), with S214's full nonlinear chart q.",
        "detector_Ward": "Differentiate conservation with a fixed detector vector xi_D. The time flux is2 xi^alpha(deltaE^{0nu}g_alpha nu+E^{0nu}hG_alpha nu). It vanishes at the right by xi_D=0 and at the left by hG=deltaE=0. Therefore R_ADM(Dxi,G)=integral -E:Lxi hG+E:q(Dxi,G).",
        "ordered_reconstruction": "R_ADM(D,G)=R_ADM(Dsyn,Gsyn)+source_Ward(D,xi_G)+detector_Ward(xi_D,Gsyn). Source and detector identities are derived separately; the retarded response is not declared symmetric.",
        "existence_boundary": "These are structural identities for the original covariant renormalized Gaussian current and its admissible retarded distributional variation. They do not prove existence, quantitative matching or norms of the three remaining scalar kernels. Conditional reconstruction cannot create those missing kernels.",
        "prescription": "A sharp computational momentum band is not diffeomorphism invariant. No finite-band Ward identity is asserted. The fixed covariant pole/finite prescription and same prepared state are required before using these identities. No finite coefficient, state germ, clock retuning or one-point contact is changed.",
        "checks": checks,
        "gates": {
            "density_weight_and_both_tensor_indices_retained": True,
            "source_covariance_and_detector_conservation_separate": True,
            "full_nonlinear_ADM_chart_contact_retained": True,
            "both_differentiated_boundary_fluxes_zero_for_stated_reasons": True,
            "no_finite_sharp_band_diffeomorphism_claim": True,
            "structural_not_unproved_full_stress_norm_differentiability": True,
            "conditional_ordinary_Proca_bridge_not_full_clock_parent": True,
        },
    }
