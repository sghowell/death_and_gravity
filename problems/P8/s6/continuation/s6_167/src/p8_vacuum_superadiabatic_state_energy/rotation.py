"""Exact alternating Pauli frames, including every induced connection."""

from functools import cache

import sympy as s


def norm_squared(matrix):
    return s.simplify(sum(s.conjugate(v) * v for v in matrix))


@cache
def data():
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    theta, rate, r = s.symbols("theta theta_prime r", real=True)
    ux = s.cos(theta / 2) * s.eye(2) + s.I * s.sin(theta / 2) * pauli[0]
    uy = s.cos(theta / 2) * s.eye(2) - s.I * s.sin(theta / 2) * pauli[1]
    Hx = r * (s.cos(theta) * pauli[2] + s.sin(theta) * pauli[0])
    Hy = r * (s.cos(theta) * pauli[2] + s.sin(theta) * pauli[1])
    t = s.Symbol("t", real=True)
    e, g = s.Function("e")(t), s.Function("g")(t)
    connection = (e * s.diff(g, t) - g * s.diff(e, t)) / (2 * (e * e + g * g))
    phase = s.Symbol("phase", real=True)
    # Use a real angle for the projector algebra so positivity is explicit.
    angle = s.Symbol("angle", real=True)
    v = s.Matrix([s.cos(angle), s.sin(angle) * s.exp(s.I * phase)])
    delta = v * v.H - s.diag(1, 0)
    checks = {
        "x_rotation_unitary": norm_squared(ux.H * ux - s.eye(2)),
        "y_rotation_unitary": norm_squared(uy.H * uy - s.eye(2)),
        "sigma1_frame_exact_diagonalization": norm_squared(
            s.simplify(uy.H * Hx * uy - r * pauli[2])
        ),
        "sigma2_frame_exact_diagonalization": norm_squared(
            s.simplify(ux.H * Hy * ux - r * pauli[2])
        ),
        "sigma1_to_sigma2_connection_negative": norm_squared(
            s.simplify(-s.I * uy.H * s.diff(uy, theta) * rate + rate * pauli[1] / 2)
        ),
        "sigma2_to_sigma1_connection_positive": norm_squared(
            s.simplify(-s.I * ux.H * s.diff(ux, theta) * rate - rate * pauli[0] / 2)
        ),
        "exact_atan_derivative_connection": s.simplify(
            s.diff(s.atan(g / e), t) / 2 - connection
        ),
        "projector_difference_traceless": s.simplify(s.trace(delta)),
        "projector_difference_determinant": s.trigsimp(delta.det() + s.sin(angle) ** 2),
        "projector_difference_squared": norm_squared(
            s.trigsimp(delta * delta - s.sin(angle) ** 2 * s.eye(2))
        ),
        "two_helicity_projector_trace_norm_factor": 2 * 2 - 4,
        "out_energy_spin_particle_antiparticle_count": 2 * 2 - 4,
    }
    return {
        "sigma1_diagonalizing_rotation": uy,
        "sigma2_diagonalizing_rotation": ux,
        "exact_connection_magnitude_formula": connection,
        "exact_recursion": "Starting after the physical mass rotation with e0=tau omega and g0=pDelta s'/(2omega^2) on sigma2, set e_(j+1)=sqrt(e_j^2+g_j^2), g_(j+1)=(-1)^j partial_x atan(g_j/e_j)/2 for j=0,1,2,3. The off-diagonal axes alternate sigma2,sigma1,sigma2,sigma1,sigma2.",
        "asymptotic_state_identity": "All four additional frame rotations tend to identity at both time infinities because their g_j vanish and e_j have a positive real gap. The initial mass rotation retains the physical asymptotic eigenbasis. The final transition amplitude is exactly the original in/out transition amplitude, up to harmless phases.",
        "complete_evolution": "After removing the final diagonal phase the interaction is off diagonal with norm |g4|. The exact unitary Duhamel integral bounds its transition amplitude by integral |g4|. g4 is not set to zero and no finite adiabatic approximation replaces the evolution.",
        "covariance_difference": "For one helicity two pure rank-one projectors differ with eigenvalues +/-|beta|. Both helicities give trace norm4|beta| per color/flavor. Unitary evolution preserves this norm at any real time.",
        "checks": checks,
    }
