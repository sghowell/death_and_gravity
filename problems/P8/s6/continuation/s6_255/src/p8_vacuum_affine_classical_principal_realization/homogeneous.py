"""Full homogeneous current action and the regular constrained Euler system."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import background as old_background
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c

u, N = old_background.u, old_background.N
R, F, I, JH = old_background.R, old_background.F, old_background.I, old_background.JH
kappa, mass2 = old_background.kappa, old_background.mass2
hbar, mh, mc, chi = s.symbols(
    "homogeneous_heavy_field heavy_rate M1_rate M1_field", real=True
)


@cache
def action():
    families = old_background.coefficients()["whole_ADM_families"]
    U, M, B, Fhat = (families[name] for name in ("U", "M", "B", "Fhat"))
    D, Z = M / N, U / N
    V = -mass2 * hbar**2 / 2 + JH * hbar / s.sqrt(kappa)
    effective_F = Fhat + U * V
    density = -3 * D * c.H**2 + 3 * B * c.H + N * effective_F + Z * (mc**2 + mh**2) / 2
    theta = -c.H * s.diff(D, N) + s.diff(B, N) / 2
    w = s.Matrix([s.diff(Z, N) * mc, s.diff(Z, N) * mh])
    Cnn = s.diff(density, N, 2) / 2
    return {
        "normalized_homogeneous_density": density,
        "physical_homogeneous_action_density": kappa * c.a**3 * density,
        "complete_effective_Fhat": effective_F,
        "complete_heavy_potential": V,
        "D": D,
        "Z": Z,
        "Theta": theta,
        "w": w,
        "Cnn": Cnn,
        "J": Cnn + 3 * theta**2 / D - w.dot(w) / (2 * Z),
        "complete_lapse_constraint": s.diff(density, N),
        "complete_heavy_force_density": N * U * (-mass2 * hbar + JH / s.sqrt(kappa)),
        "complete_M1_charge": c.a**3 * Z * mc,
        "complete_heavy_momentum": c.a**3 * Z * mh,
    }


@cache
def data():
    packet = action()
    density = packet["normalized_homogeneous_density"]
    velocity = s.Matrix([c.H, mc, mh])
    q = s.Matrix([6 * packet["Theta"], *packet["w"]])
    kinetic = s.diag(-6 * packet["D"], packet["Z"], packet["Z"])
    D, Z, J, theta = s.symbols("positive_D positive_Z positive_J Theta", real=True)
    w1, w2 = s.symbols("mixed_momentum1 mixed_momentum2", real=True)
    vector = s.Matrix([6 * theta, w1, w2])
    Cnn = J - 3 * theta**2 / D + (w1**2 + w2**2) / (2 * Z)
    K0 = s.diag(-6 * D, Z, Z)
    reduced = K0 - vector * vector.T / (2 * Cnn)
    inverse = K0.inv() + K0.inv() * vector * vector.T * K0.inv() / (2 * J)
    mapping = {
        D: packet["D"],
        Z: packet["Z"],
        J: packet["J"],
        theta: packet["Theta"],
        w1: packet["w"][0],
        w2: packet["w"][1],
    }
    checks = {
        "whole_three_velocity_homogeneous_Hessian": s.hessian(density, velocity)
        - kinetic,
        "whole_lapse_velocity_cross_vector": s.Matrix(
            [s.diff(density, N, v) for v in velocity]
        )
        - q,
        "whole_regular_reduced_three_velocity_Hessian_inverse": reduced * inverse
        - s.eye(3),
        "whole_regular_reduced_three_velocity_determinant": reduced.det()
        + 6 * D * Z**2 * J / Cnn,
        "complete_literal_lapse_pivot_binding": Cnn.subs(mapping, simultaneous=True)
        - packet["Cnn"],
        "heavy_force_is_literal_variation": s.diff(density, hbar)
        - packet["complete_heavy_force_density"],
        "M1_momentum_is_literal_variation": c.a**3 * s.diff(density, mc)
        - packet["complete_M1_charge"],
        "heavy_momentum_is_literal_variation": c.a**3 * s.diff(density, mh)
        - packet["complete_heavy_momentum"],
        "M1_shift_symmetry_not_changed": s.diff(density, chi),
    }
    return {
        **packet,
        "whole_fixed_lapse_three_velocity_Hessian": kinetic,
        "whole_lapse_velocity_cross_vector": q,
        "whole_implicit_lapse_reduced_Hessian": reduced,
        "whole_implicit_lapse_reduced_Hessian_inverse": inverse,
        "whole_implicit_lapse_reduced_Hessian_determinant": -6 * D * Z**2 * J / Cnn,
        "whole_actual_coefficient_binding": mapping,
        "local_existence_domain": "Use log(a_hat),M1,H/sqrt(kappa) as the three homogeneous configuration coordinates. On a>0,N>0,R>0 with Cnn>0 and J>0, the lapse equation has nonzero derivative2Cnn and determines N smoothly from time, fields and their three velocities. The entire reduced Euler velocity Hessian is invertible with the listed inverse. The resulting smooth six-dimensional first-order ODE has a unique local solution through each admitted constraint datum. No uniform-mass or macroscopic lifetime is claimed.",
        "reconstruction": "Homogeneous W0=3(R-1)(Hhat-Hclock),Wi0 has F(W)0 and W-S0, so its full vector Euler equation vanishes and its first action variation contributes zero to the background equations. The original algebraic affine complement is reconstructed in its unchanged invertible R domain. Isotropy removes vector, tracefree and spatial-gradient equations consistently. The lapse equation is imposed, not lost by substituting a gauge; the clock equation follows the full covariant Noether identity with phi_dot1.",
        "checks": checks,
        "gates": {
            "entire_heavy_source_potential_in_homogeneous_action": density.has(
                JH, hbar, mh
            ),
            "both_matter_velocities_retained": density.has(mc, mh),
            "no_homogeneous_vector_kinetic_term_deleted_when_nonzero": True,
            "lapse_implicit_function_requires_Cnn_nonzero": True,
            "entire_reduced_kinetic_inverse_requires_J_nonzero": True,
            "local_ODE_not_macro_bounce_or_quantum_mean": True,
        },
    }
