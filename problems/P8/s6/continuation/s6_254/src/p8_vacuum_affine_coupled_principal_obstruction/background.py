"""Entire current off-clock parent and the heavy source in the four-mode block."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_source_filtration import source as heavy_source
from p8_vacuum_affine_physical_background_vertices import parent

from . import coupled as c

u, N = parent.u, parent.N
R = s.Function("entire_current_R", positive=True)(u, N)
F = s.Function("entire_current_F", real=True)(u, N)
I = s.Function("current_ADM_primitive_fixed_at_N1", real=True)(u, N)
JH = s.Function("entire_physical_heavy_source", real=True)(u, N)
kappa, mass2 = s.symbols(
    "actual_positive_kappa actual_positive_heavy_mass2", positive=True
)
Nd, md, Hd, m = s.symbols(
    "background_Ndot background_mdot background_Hdot background_m", real=True
)


def dtime(expression):
    return (
        s.diff(expression, u)
        + Nd * s.diff(expression, N)
        + c.a * c.H * s.diff(expression, c.a)
        + md * s.diff(expression, m)
        + Hd * s.diff(expression, c.H)
    )


@cache
def coefficients():
    U, M, C3, Cchi = (
        R ** -s.Rational(3, 4),
        R ** s.Rational(1, 4),
        R ** s.Rational(3, 4) / 2,
        R ** -s.Rational(1, 4),
    )
    Ru = s.diff(R, u)
    B = -U * Ru / (2 * N) - I
    Fhat = U * (F + 9 * Ru**2 / (16 * R * N**2)) - s.diff(I, u) / N
    D, Z = M / N, U / N
    theta = -c.H * s.diff(D, N) + s.diff(B, N) / 2
    charge, w = Z * m, s.diff(Z, N) * m
    base0 = -3 * D * c.H**2 + 3 * B * c.H + N * Fhat + Z * m**2 / 2
    cv = -18 * D * c.H + 9 * B
    Cnn = s.diff(base0, N, 2) / 2
    Qn = s.diff(base0, N)
    source_density = N * U * JH / s.sqrt(kappa)
    return {
        "whole_parameter_substitution": {
            c.D: D,
            c.Z: Z,
            c.J: Cnn + 3 * theta**2 / D - w**2 / (2 * Z),
            c.th: theta,
            c.r: R - 1,
            c.rN: s.diff(R, N),
            c.dH: c.H - parent.H,
            c.L0: 3 * Qn,
            c.L2: 4 * s.diff(N * C3, N),
            c.C: N * C3,
            c.Y: N * Cchi,
            c.K: c.zeta * Cchi / (N * c.a**2),
            c.Yv: N * Cchi / c.a**2,
            c.Vvv: 9 * base0 / 2 - (dtime(cv) + 3 * c.H * cv) / 2,
            c.cs[0]: charge,
            c.cs[1]: 0,
            c.ws[0]: w,
            c.ws[1]: 0,
            c.ds[0]: 0,
            c.ds[1]: s.diff(source_density, N),
            c.vs[0]: -3 * (dtime(charge) + 3 * c.H * charge),
            c.vs[1]: 3 * source_density,
            c.m00: 0,
            c.m01: 0,
            c.m11: -N * U * mass2,
        },
        "whole_heavy_source_density": source_density,
        "whole_ADM_families": {
            "M": M,
            "U": U,
            "C3": C3,
            "Cchi": Cchi,
            "B": B,
            "Fhat": Fhat,
        },
        "primitive_lapse_derivative": 3 * U * Ru * s.diff(R, N) / (4 * R * N),
        "primitive_value_on_clock": 0,
        "primitive_interpretation": "I(u,N) is the exact integral from1 toN of its listed lapse derivative. I(u,1)=0. All mixed derivatives in the coefficient map act on that same integral, not on an independently varied function or truncated clock jet.",
    }


@cache
def data():
    packet = coefficients()
    full = parent.fixed_functions()
    source = heavy_source.physical_source(u, N**-2)
    bindings = {
        R: full["R_full"].subs(parent.X, N**-2),
        F: full["F_full"].subs(parent.X, N**-2),
        JH: source,
        kappa: heavy_source.KAPPA,
        mass2: heavy_source.MASS2,
    }
    # Independent scalar density variation in an arbitrary lapse/volume chart.
    x, n, v, h, hd, P, scale, eps = s.symbols(
        "x lapse_fluct volume_fluct h hdot FourierP scale eps", real=True
    )
    fU, fY, fJ = (
        s.Function(name)(x)
        for name in ("generic_U", "generic_Cchi", "generic_physical_J")
    )
    shifted = x + eps * n
    kinetic = fU.subs(x, shifted) * s.exp(3 * eps * v) * (eps * hd) ** 2 / (2 * shifted)
    gradient = (
        -shifted
        * fY.subs(x, shifted)
        * s.exp(eps * v)
        * P**2
        * (eps * h) ** 2
        / (2 * scale**2)
    )
    potential = (
        -shifted * fU.subs(x, shifted) * s.exp(3 * eps * v) * mass2 * (eps * h) ** 2 / 2
    )
    source_term = (
        shifted
        * fU.subs(x, shifted)
        * s.exp(3 * eps * v)
        * fJ.subs(x, shifted)
        * eps
        * h
        / s.sqrt(kappa)
    )
    all_terms = kinetic + gradient + potential + source_term
    linear = s.diff(all_terms, eps).subs(eps, 0).doit()
    quadratic = s.diff(all_terms, eps, 2).subs(eps, 0).doit() / 2
    density = x * fU * fJ / s.sqrt(kappa)
    expected = (
        fU * hd**2 / (2 * x)
        - x * fY * P**2 * h**2 / (2 * scale**2)
        - x * fU * mass2 * h**2 / 2
        + (s.diff(density, x) * n + 3 * density * v) * h
    )
    RN = parent.lapse_jets()["R_lapse_jets_zero_through_four"][1]
    ref = parent.reference_data()["whole_reference_coefficients"]
    J0 = ref["Jc"].subs(u, 0)
    r1, L20 = RN.subs(u, 0), (2 * ref["E"]).subs(u, 0)
    eta = s.Symbol("nonzero_lapse_probe", real=True)
    mapped = packet["whole_parameter_substitution"]
    checks = {
        "whole_heavy_density_linear_source_not_erased": linear - density * h,
        "whole_heavy_density_quadratic_source_and_kinetic_terms": quadratic - expected,
        "complete_heavy_lapse_source_in_full_block": mapped[c.ds[1]]
        - s.diff(packet["whole_heavy_source_density"], N),
        "complete_heavy_volume_source_in_full_block": mapped[c.vs[1]]
        - 3 * packet["whole_heavy_source_density"],
        "complete_heavy_mass_in_full_block": mapped[c.m11]
        + N * R ** -s.Rational(3, 4) * mass2,
        "actual_Hbar0_heavy_rate_zero": mapped[c.cs[1]] + mapped[c.ws[1]],
        "complete_heavy_source_clock_value_zero": source.subs(N, 1),
        "full_current_reference_R_lapse_derivative": RN + 2 / parent.h,
        "actual_central_lapse_derivative": r1 + 2,
        "actual_central_lapse_spatial_coefficient": L20 + 1,
        "actual_central_full_current_pivot": J0
        - s.Rational(243, 160)
        - parent.RHO.subs(u, 0)
        + 7 * parent.PRESSURE.subs(u, 0) / 8,
        "actual_nonreference_fast_coefficient_second_jet": 2
        * r1**2
        * eta**2
        * L20**2
        / (8 * J0)
        - eta**2 / J0,
        "exact_literal_gamma_domain": c.GAMMA.subs(
            {
                c.D: R ** s.Rational(1, 4) / N,
                c.Z: R ** -s.Rational(3, 4) / N,
                c.r: R - 1,
            },
            simultaneous=True,
        )
        - 1
        + 3 * (R - 1) ** 2 / (2 * R),
        "literal_R_strip_gamma_margin": 1
        - 3 * (R - 1) ** 2 / (2 * R)
        - s.Rational(1, 4)
        - 3 * (R - s.Rational(1, 2)) * (2 - R) / (2 * R),
    }
    return {
        **packet,
        "entire_current_function_bindings": bindings,
        "whole_fixed_profile_bindings": {
            parent.RHO: full["rho"],
            parent.PRESSURE: full["pressure"],
        },
        "all_three_fixed_vacuum_constants": full["vacuum_pressure_total"],
        "whole_heavy_quadratic_density": expected,
        "nonzero_offshell_heavy_linear_density": linear,
        "actual_reference_central_pivot": J0,
        "actual_fast_rate_fourth_power_second_lapse_jet": eta**2 / J0,
        "actual_probe": "N=1+epsilon eta(u), with eta constant nonzero on a small interval aboutu0 and zero near the original initial slice; a_phys=(1+u^2)^2; a_hat=R(u,N)^(1/4)a_phys; m=ell(u); Hbar=0; W0=3(R-1)(Hhat-Hclock). This is a smooth homogeneous OFF-SHELL background, not a constructed unforced nonlinear solution. The full R,F,J bindings, both profiles and every constant remain. Atu0, R_N=-2,L2=-1,Jc=243/160+rho-7P/8>0 by the unchanged reference bounds. Continuity gives a compact interval with strict margins for each sufficiently small fixed nonzero epsilon.",
        "off_reference_source_scope": "The degree-eight/1024 clock filtration licenses S253's finite probe vertices but does not delete this full off-clock H source. AtHbar0 it adds the stated lapse-H and volume-H contacts and heavy mass, with no higher spatial derivative. These coefficients are retained in the entire four-mode Hamiltonian; only then are they shown lower order in its large-P symbol.",
        "checks": checks,
        "gates": {
            "entire_heavy_source_is_not_identically_zero": source != 0,
            "entire_fixed_finite_onepoint_switch_retained": source.has(s.pi),
            "full_current_R_F_not_clock_germs_off_reference": bindings[R]
            != full["R_clock"].subs(parent.X, N**-2),
            "all_heavy_quadratic_mixed_sources_retained": mapped[c.ds[1]] != 0
            and mapped[c.vs[1]] != 0,
            "actual_central_r_first_derivative_nonzero": r1 != 0,
            "actual_central_L2_nonzero": L20 != 0,
            "unchanged_reference_J_positivity_not_new_mean": True,
            "homogeneous_aligned_Hbar0_background_scope_only": True,
            "no_on_shell_bounce_or_finite_cutoff_refutation_inferred": True,
        },
    }
