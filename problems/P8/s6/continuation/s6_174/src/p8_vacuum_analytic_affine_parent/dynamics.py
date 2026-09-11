"""Physical clock/vacuum quadratic blocks and regular nonlinear auxiliary pivots."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family


@cache
def clock_lapse():
    u, X = family.u, family.X
    v, N = s.symbols("s N", positive=True)
    h = (1 + u * u) ** 3
    H = 4 * u / (1 + u * u)
    # These are source-pinned local clock jets, not a replacement of the global action.
    R = 1 + (X - 1) / h
    B = -R / 2
    fphi = -s.diff(R, u) / 2
    ox = -s.diff(R, X) / (4 * R)
    op = -s.diff(R, u) / (4 * R)
    F = family.original.data()["original_retuned_tree_scalar"]
    volume = R ** -s.Rational(3, 4)
    linear = 4 * B * v * op + 2 * v * fphi
    potential = F + 6 * B * v * v * op * op + 6 * v * v * fphi * op
    primitive = volume * 12 * v * v * fphi * ox
    primitive = primitive.subs(X, v * v)
    IN = -primitive.subs(v, 1)
    INN = (s.diff(primitive, v) + 2 * primitive).subs(v, 1)
    I = IN * (N - 1) + INN * (N - 1) ** 2 / 2
    jet = lambda e: s.expand(
        sum(
            s.diff(e, N, j).subs(N, 1) * (N - 1) ** j / s.factorial(j) for j in range(3)
        )
    )
    vol = jet(volume.subs(X, N**-2))
    a = jet((s.Rational(2, 3) * volume * B).subs(X, N**-2))
    b = jet((volume * linear).subs({X: N**-2, v: 1 / N}, simultaneous=True) - I)
    f = jet(
        (volume * potential).subs({X: N**-2, v: 1 / N}, simultaneous=True)
        - s.diff(I, u) / N
    )
    w = 1 / (10 * (1 + u * u) ** 6)
    p = -2 * H
    Ham = N * ((p - b) ** 2 / (4 * a) - f + w * w / (2 * vol))
    first = s.factor(s.diff(Ham, N).subs(N, 1))
    second = s.factor(s.diff(Ham, N, 2).subs(N, 1))
    J = s.factor(-second / 2)
    numerator, denominator = s.fraction(J)
    z = s.symbols("nonnegative_u_squared", nonnegative=True)
    P = s.Poly(numerator, u)
    if any(power[0] % 2 for power in P.monoms()):
        raise ValueError("Expected an even positive-coefficient clock lapse numerator")
    Pz = s.Poly(sum(coeff * z ** (power[0] // 2) for power, coeff in P.terms()), z)
    # Joint temporal source vanishes with its first variation; W0=0 and j=0.
    return {
        "u": u,
        "J": J,
        "joint_lapse_temporal_Jacobian": s.diag(-2 * J, -1),
        "positive_lapse_numerator_in_u_squared": Pz.as_expr(),
        "lapse_numerator_coefficients": tuple(Pz.all_coeffs()),
        "lapse_denominator": denominator,
        "nonlinear_count": "At each finite clock point, continuity, the nonzero ten-velocity block and the nonsingular joint lapse/normal-vector Jacobian give a nonempty local canonical neighborhood with seven physical modes (2 tensor,2 light scalar,3 vector). No uniform neighborhood size or stability throughout it is supplied.",
        "gates": {
            "clock_lapse_polynomial_coefficients_positive": all(
                c > 0 for c in Pz.all_coeffs()
            ),
            "clock_lapse_denominator_positive": denominator.is_positive is True,
        },
        "checks": {
            "actual_target_clock_lapse_equation": first,
            "boundary_I_N_zero": s.simplify(IN),
            "boundary_I_NN_retained": s.factor(INN + 3 * s.diff(h, u) / h**3),
        },
    }


@cache
def data():
    r, H, Hd, q, m2 = s.symbols("r H Hdot q m_squared", real=True)
    A, Ad, A0, k, a, mu = s.symbols("A Ad A0 k a positive_mu", real=True)
    # Real Fourier longitudinal mode, with a harmless phase convention for A0.
    L = a * (Ad - k * A0) ** 2 / 2 + a**3 * mu * A0 * A0 / 2 - a * mu * A * A / 2
    A0sol = k * Ad / (k * k + a * a * mu)
    reduced = s.factor(L.subs(A0, A0sol))
    expected = a * mu * Ad * Ad / (2 * (k * k / a**2 + mu)) - a * mu * A * A / 2
    l = H * (1 + 2 * r)
    rdot = -2 * H * r * (1 - r)
    ld = Hd * (1 + 2 * r) + 2 * H * rdot
    omegaL = q + m2 - ld / 2 - l * l / 4
    omegaT = q + m2 - Hd / 2 - H * H / 4
    expectedL = (
        q
        + m2
        - (s.Rational(1, 2) + r) * Hd
        + H * H * (-s.Rational(1, 4) + r - 3 * r * r)
    )
    R, vol, trace, T, b, f, p, j = s.symbols("R volume K T b f p j", real=True)
    delta = R - 1
    aa = -vol * R / 3
    c = s.symbols("source_lower", real=True)
    full = aa * trace * trace + b * trace + f + vol * (T - delta * trace - c) ** 2 / 2
    ksol = (p - b - delta * j) / (2 * aa)
    tsol = delta * ksol + c - j / vol
    joint = s.factor(
        (p * trace - full - T * j).subs({trace: ksol, T: tsol}, simultaneous=True)
    )
    target = (p - b - delta * j) ** 2 / (4 * aa) - f + j * j / (2 * vol) - c * j
    effective = 1 + vol * delta * delta / (2 * aa)
    return {
        "canonical_vacuum_blocks": "At u=0,du=0,W=0,chi=constant: Einstein gravity plus mass-one canonical Phi=sqrt(kappa)u, a canonical massless chi, and canonical A=sqrt(kappa*zeta)W with mass_squared=1/zeta. There is no quadratic light-vector mixing.",
        "vacuum_source_first_degree": 4,
        "vacuum_leading_vector_light_vertex_degree": 5,
        "clock_quadratic_action": "Exactly the unchanged retuned CD/M1 light quadratic action plus a positive, minimally physical-metric-coupled Proca block. The additional source is quadratic in deviations; background and first variation vanish.",
        "vector_mass_squared": "1/zeta>=2000",
        "longitudinal_eliminated_positive_kinetic": a * mu / (k * k / a**2 + mu),
        "canonical_longitudinal_frequency": s.factor(omegaL),
        "canonical_transverse_frequency": omegaT,
        "uniform_frequency_floor": "For all real CD times and q=k_com^2/a^2>=0, Omega_L^2>=q+1/zeta-15>=q+1985; Omega_T^2>=q+1/zeta-3. This is a time-dependent oscillator floor, not a physical EFT cutoff or a stationary scattering gap.",
        "full_trace_temporal_Hamiltonian_after_elimination": target,
        "normalized_temporal_schur": effective,
        "checks": {
            "longitudinal_temporal_Euler": s.factor(s.diff(L, A0).subs(A0, A0sol)),
            "full_longitudinal_kinetic_after_constraint": s.factor(reduced - expected),
            "canonical_frequency_all_normalization_derivatives": s.factor(
                omegaL - expectedL
            ),
            "longitudinal_H2_floor_identity": s.factor(
                (-s.Rational(1, 4) + r - 3 * r * r)
                + s.Rational(9, 4)
                - (1 - r) * (3 * r + 2)
            ),
            "homogeneous_longitudinal_matches_transverse": s.factor(
                omegaL.subs(r, 0) - omegaT
            ),
            "joint_trace_Euler": s.factor(
                s.diff(full, trace).subs({trace: ksol, T: tsol}, simultaneous=True) - p
            ),
            "joint_temporal_Euler": s.factor(
                (-s.diff(full, T) - j).subs({trace: ksol, T: tsol}, simultaneous=True)
            ),
            "full_joint_trace_temporal_Legendre": s.factor(joint - target),
            "same_uniform_temporal_schur": s.factor(
                effective - 1 + 3 * (R - 1) ** 2 / (2 * R)
            ),
        },
    }
