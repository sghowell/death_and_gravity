"""Literal nonzero-transfer ADM scalar action and regular coefficient-sector Hamiltonian."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_response import classical as tree
from p8_vacuum_affine_quantum_retuning import profile

t = tree.u
a = (1 + t * t) ** 2
H, delta, ell = tree.H, tree.delta, tree.ell
n, v, vd, sigma, sd, b, q = s.symbols("n v vdot sigma sigmadot b q", real=True)
pv, ps = s.symbols("p_v p_sigma", real=True)
J0, Theta, w, E, dJ, Tc, A = s.symbols("J Theta w E DeltaJ Tcorr A", real=True)
Z = s.Matrix([v, sigma, pv, ps])
JC = s.zeros(4)
JC[:2, 2:] = s.eye(2)
JC[2:, :2] = -s.eye(2)


def lagrangian(J=J0, theta=Theta, matter=w, charge=ell, gradient=E):
    return (
        -3 * vd**2
        + (J + matter**2 / 2 - 3 * theta**2) * n**2
        + 6 * theta * n * vd
        + sd**2 / 2
        + matter * n * sd
        - 3 * charge * vd * sigma
        + 2 * b * (vd - theta * n)
        + charge * b * sigma
        + q * v**2
        + 2 * gradient * q * n * v
        - q * sigma**2 / 2
    )


def lapse_numerator():
    return Theta * pv - w * ps + 3 * Theta * ell * sigma - (2 * E * q + 3 * Tc) * v


def reduced_hamiltonian():
    return (
        ps**2 / 2
        - ell * pv * sigma / 2
        + (q / 2 - 3 * ell**2 / 4) * sigma**2
        - q * v**2
        - s.Rational(9, 2) * A * v**2
        + lapse_numerator() ** 2 / (4 * (J0 + dJ))
    )


def actual_coefficients():
    base = tree.quadratic()
    aa, bb = -profile.P, -(profile.rho + profile.P) / 2
    return {
        J0: base["J"],
        Theta: base["theta"],
        w: base["w"],
        E: 1 - 3 * delta,
        A: aa,
        Tc: (1 + 3 * delta) * aa - 2 * bb,
        dJ: (21 * delta**2 - 3 * delta) * aa / 2 + (1 - 6 * delta) * bb,
    }


@cache
def spatial_data():
    rows = tree.clock_coefficients()["coefficient_clock_N_jets"]
    base = tree.quadratic()
    eps, vx, vxx, beta, sx = s.symbols("epsilon vx vxx beta sigmax", real=True)
    N = 1 + eps * n
    coef = {
        key: sum(vals[j] * (eps * n) ** j / s.factorial(j) for j in range(3))
        for key, vals in rows.items()
    }
    vol = 1 + 3 * eps * v + s.Rational(9, 2) * eps**2 * v * v
    M2 = -3 * coef["a"]
    C3 = s.Rational(1, 2) - s.Rational(3, 2) * delta * eps * n
    rate = H + eps * vd - eps**2 * beta * vx
    literal = vol * (
        -3 * M2 * rate**2 / N
        + 2 * M2 * rate * eps * b / N
        + coef["b"] * (3 * rate - eps * b)
        + N * coef["f"]
        + N
        * C3
        * (1 - 2 * eps * v + 2 * eps**2 * v * v)
        / a**2
        * (-4 * eps * vxx - 2 * eps**2 * vx**2)
        + coef["U"] * (ell + eps * sd - eps**2 * beta * sx) ** 2 / (2 * N)
        - eps**2 * sx * sx / (2 * a * a)
    )
    raw = s.expand(s.diff(literal, eps, 2).subs(eps, 0) / 2)
    vv = s.factor(raw.coeff(v, 1).coeff(vd, 1))
    normalized = s.expand(raw - vv * v * vd - tree.weighted(vv) * v * v / 2)
    shiftboundary = s.factor(normalized.coeff(beta, 1).coeff(vx, 1))
    checks = {
        "shift_metric_spatial_divergence": s.factor(
            normalized.coeff(v, 1).coeff(b, 1) - shiftboundary
        ),
        "metric_curvature_spatial_boundary": s.factor(
            normalized.coeff(v, 1).coeff(vxx, 1) + 2 / a**2
        ),
        "full_Theta_without_crossing_division": s.factor(
            base["theta"] - H * (1 + delta) - s.diff(delta, t)
        ),
        "Theta_explicit_regular_clock": s.factor(
            base["theta"] - H + t / (1 + t * t) ** 4
        ),
        "matter_weighted_charge_identity": s.factor(s.diff(ell, t) + 3 * H * ell),
    }
    normalized = s.expand(normalized - shiftboundary * (beta * vx + v * b))
    sb = s.factor(normalized.coeff(beta, 1).coeff(sx, 1))
    normalized += -sb * beta * sx - sb * b * sigma
    curvature = s.factor(s.expand(normalized).coeff(v, 1).coeff(vxx, 1))
    normalized += -curvature * v * vxx - curvature * vx * vx
    normalized += -3 * ell * v * sd - 3 * ell * vd * sigma
    normalized = s.expand(
        normalized.subs(
            {
                vx**2: a * a * q * v * v,
                sx**2: a * a * q * sigma * sigma,
                vxx: -a * a * q * v,
            },
            simultaneous=True,
        )
    )
    expected = lagrangian(base["J"], base["theta"], base["w"], ell, 1 - 3 * delta)
    checks["complete_literal_inhomogeneous_scalar_ADM_action"] = s.factor(
        normalized - expected
    )
    checks["nonzero_transfer_shift_constraint"] = s.diff(expected, b) - 2 * (
        vd - base["theta"] * n + ell * sigma / 2
    )
    return {
        "full_quadratic_per_kappa_a_cubed": expected,
        "theta": base["theta"],
        "J": base["J"],
        "w": base["w"],
        "boundary_convention": "Full S174 primitive, including I_NN, is already in S180 coefficient jets. Integrate spatial divergences, the weighted v vdot term, and 3ell v sigmadot to -3ell vdot sigma. These boundaries retain the original matter wave.",
        "nonzero_transfer": "hhat=a^2 exp(2v)I, beta=grad B/a^2, b=Delta B/a^2, q=|P|^2/a^2. Scalar isotropy reduces the literal quadratic spatial calculation to a Fourier axis; no homogeneous fixed-charge Routh constraint is imposed.",
        "checks": checks,
        "gates": {
            "complete_spatial_curvature_lapse_and_shift_terms": True,
            "original_matter_gradient_and_velocity_retained": True,
            "no_Theta_or_H_division": True,
            "literal_global_zero_mode_is_separate": True,
            "bare_action_is_only_comparison_input_until_retuning_added": True,
        },
    }


@cache
def retuning_data():
    D, Bc, N, eps = s.symbols("delta B N epsilon", real=True)
    R = 1 + 2 * D * (N**-2 - 1)
    F = A + Bc * (N**-2 - 1)
    density = N * R ** -s.Rational(3, 4) * F
    shift = (21 * D**2 - 3 * D) * A / 2 + (1 - 6 * D) * Bc
    tadpole = (1 + 3 * D) * A - 2 * Bc
    literal = s.exp(3 * eps * v) * density.subs(N, 1 + eps * n)
    quadratic = shift * n * n + 3 * tadpole * n * v + s.Rational(9, 2) * A * v * v
    checks = {
        "literal_fixed_profile_lapse_tadpole": s.factor(
            s.diff(density, N).subs(N, 1) - tadpole
        ),
        "literal_fixed_profile_lapse_pivot": s.factor(
            s.diff(density, N, 2).subs(N, 1) / 2 - shift
        ),
        "full_volume_and_lapse_retuning_quadratic": s.factor(
            s.diff(literal, eps, 2).subs(eps, 0) / 2 - quadratic
        ),
        "full_retuning_linear_volume_and_lapse": s.factor(
            s.diff(literal, eps).subs(eps, 0) - (3 * A * v + tadpole * n)
        ),
        "actual_tadpole_fixed_stress": s.factor(
            tadpole.subs({A: -profile.P, Bc: -(profile.rho + profile.P) / 2, D: delta})
            - profile.rho
            + 3 * delta * profile.P
        ),
        "actual_bounce_pivot_shift": s.factor(
            shift.subs(
                {A: -profile.P, Bc: -(profile.rho + profile.P) / 2, D: s.Rational(1, 2)}
            )
            - profile.rho
            + s.Rational(7, 8) * profile.P
        ),
    }
    return {
        "fixed_QG1_clock_profile": "The original full S182 switch profile is fixed once from rho_ref/kappa and P_ref/kappa. Clock X jets are A=-p, B=-(r+p)/2 and zero at orders2..1023. This local second variation does not replace that global profile or recompute stress on a varied history.",
        "quadratic_addition": quadratic,
        "DeltaJ": shift,
        "Tcorr": tadpole,
        "candidate_boundary": "This is the current QG1 classical coefficient sector. Its one-point terms alone do not vanish: they cancel the actual reference Gaussian one-point function in the retained combined equations. The Gaussian metric response must still be added, including chart contacts.",
        "checks": checks,
        "gates": {
            "full_fixed_retuning_not_only_lapse_pivot": True,
            "mixed_n_v_and_v_squared_terms_retained": True,
            "fixed_coefficient_not_live_stress_feedback": True,
            "no_old_isolated_homogeneous_quantum_inverse_used": True,
        },
    }


@cache
def hamiltonian_data():
    L = lagrangian() + dJ * n * n + 3 * Tc * n * v + s.Rational(9, 2) * A * v * v
    velocities = {
        vd: (-pv + 6 * Theta * n - 3 * ell * sigma + 2 * b) / 6,
        sd: ps - w * n,
    }
    Hfull = s.factor((pv * vd + ps * sd - L).subs(velocities, simultaneous=True))
    Hshift = s.factor(Hfull.subs(b, pv / 2))
    lc = lapse_numerator()
    reduced = reduced_hamiltonian()
    nsol = lc / (2 * (J0 + dJ))
    checks = {
        "first_canonical_momentum": s.factor(
            s.diff(L, vd) + 6 * vd - 6 * Theta * n + 3 * ell * sigma - 2 * b
        ),
        "second_canonical_momentum": s.factor(s.diff(L, sd) - sd - w * n),
        "regular_shift_pivot": s.factor(s.diff(Hfull, b) - (pv - 2 * b) / 3),
        "full_retuned_lapse_pivot_and_source": s.factor(
            s.diff(Hshift, n) + 2 * (J0 + dJ) * n - lc
        ),
        "full_shift_then_lapse_Legendre_reduction": s.factor(
            Hshift.subs(n, nsol) - reduced
        ),
        "scalar_velocity": s.factor(
            s.diff(reduced, pv) - Theta * nsol + ell * sigma / 2
        ),
        "matter_velocity": s.factor(s.diff(reduced, ps) - ps + w * nsol),
        "shift_constraint_recovered": s.factor(
            s.diff(reduced, pv) - Theta * nsol + ell * sigma / 2
        ),
        "bare_limit": s.factor(
            reduced.subs({dJ: 0, Tc: 0, A: 0})
            - (
                ps**2 / 2
                - ell * pv * sigma / 2
                + (q / 2 - 3 * ell**2 / 4) * sigma**2
                - q * v * v
                + (Theta * pv - w * ps + 3 * Theta * ell * sigma - 2 * E * q * v) ** 2
                / (4 * J0)
            )
        ),
        "canonical_symplectic_square": JC * JC + s.eye(4),
    }
    Hessian = s.hessian(reduced, Z)
    checks["quadratic_Hamiltonian_homogeneity"] = s.factor(
        (Z.T * Hessian * Z)[0] / 2 - reduced
    )
    return {
        "phase_coordinates": Z,
        "normalized_momenta": s.Matrix(
            [-6 * vd + 6 * Theta * n - 3 * ell * sigma + 2 * b, sd + w * n]
        ),
        "shift_solution": pv / 2,
        "lapse_solution": nsol,
        "reduced_coefficient_sector_Hamiltonian": reduced,
        "weighted_phase_equations": JC * Hessian * Z - 3 * H * s.diag(0, 0, 1, 1) * Z,
        "crossing": "Theta=t[4(1+t^2)^3-1]/(1+t^2)^4 vanishes at the bounce. The only lapse denominator is2(J+DeltaJ), whose original unit-slab lower bound is1/50. No crossing denominator is used.",
        "inverse_boundary": "This reduces the classical coefficient-sector quadratic comparison, not the full quantum constraints. The Gaussian current has nonlocal lapse response, and eliminating the full coupled constraint requires another inverse/Schur argument. A regular phase Hamiltonian is not a coercive energy or uniform continuum propagator.",
        "checks": checks,
        "gates": {
            "both_scalar_phase_pairs_retained": True,
            "nonzero_momentum_shift_eliminated_before_lapse": True,
            "current_QG1_fixed_retuning_fully_retained": True,
            "weighted_momentum_damping_not_omitted": True,
            "no_frozen_frequency_stability_inference": True,
            "no_full_quantum_constraint_elimination_claim": True,
        },
    }
