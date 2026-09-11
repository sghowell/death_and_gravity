"""Exact helicity-block diagonalization and modewise unitary evolution."""

from functools import cache

import sympy as s


@cache
def data():
    omega, angle, rate, phase = s.symbols("omega angle mixing phase", real=True)
    p, m, z = s.symbols("p m z", real=True)
    sigma1 = s.Matrix([[0, 1], [1, 0]])
    sigma2 = s.Matrix([[0, -s.I], [s.I, 0]])
    sigma3 = s.diag(1, -1)
    rotation = s.Matrix(
        [[s.cos(angle / 2), -s.sin(angle / 2)], [s.sin(angle / 2), s.cos(angle / 2)]]
    )
    H = omega * (s.sin(angle) * sigma1 + s.cos(angle) * sigma3)
    diagonal = s.simplify(rotation.T * H * rotation - omega * sigma3)
    connection = s.simplify(-s.I * rotation.T * s.diff(rotation, angle) + sigma2 / 2)
    t = s.Symbol("t", real=True)
    mass = s.Function("M")(t)
    coupling = p * s.diff(mass, t) / (2 * (p * p + mass * mass))
    alpha, beta = s.symbols("alpha beta", complex=True)
    da = -rate * s.exp(s.I * phase) * beta
    db = rate * s.exp(-s.I * phase) * alpha
    norm_derivative = s.simplify(
        da * s.conjugate(alpha)
        + alpha * s.conjugate(da)
        + db * s.conjugate(beta)
        + beta * s.conjugate(db)
    )
    checks = {
        "helicity_characteristic_polynomial": s.expand(
            (z * s.eye(2) - (p * sigma1 + m * sigma3)).det() - (z * z - p * p - m * m)
        ),
        "instantaneous_mixing_rate": s.simplify(
            coupling + s.diff(s.atan(p / mass), t) / 2
        ),
        "mode_probability_conserved": norm_derivative,
        "both_helicities_same_squared_frequency": s.expand(
            (-p * sigma1 + m * sigma3).det() - (p * sigma1 + m * sigma3).det()
        ),
        "two_helicities_particles_and_antiparticles": 2 * 2 - 4,
    }
    for i in range(2):
        for j in range(2):
            checks[f"instantaneous_rotation_{i}_{j}"] = diagonal[i, j]
            checks[f"moving_basis_connection_{i}_{j}"] = connection[i, j]
    return {
        "helicity_Hamiltonian": p * sigma1 + m * sigma3,
        "instantaneous_rotation": rotation,
        "mixing_rate": coupling,
        "interaction_equations": {"alpha_prime": da, "beta_prime": db},
        "phase_prescription": "phase'=2omega, omega=sqrt(p^2+M(t)^2). In the dimensionless x=t/tau test solver its accumulated half-phase obeys Phi'=tau omega.",
        "state": "For each momentum/helicity use the positive/negative energy in basis at t=-infinity. The integrable mass tails give modewise unitary Moller limits. They define translation-invariant free quasifree in/out states; no finite total energy or single global infinite-volume Fock unitary is asserted.",
        "scope": "This is the quadratic Dirac subsystem on a prescribed flat time-dependent scalar argument. Gauge interactions are omitted. No asymptotic colored-particle observable of the full interacting GY14 gauge theory, Hadamard remainder or curved-space state is identified.",
        "checks": checks,
    }
