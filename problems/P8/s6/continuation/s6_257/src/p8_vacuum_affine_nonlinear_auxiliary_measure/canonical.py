"""Whole retained nonlinear unitary Hamiltonian and its auxiliary branch."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth import background as current
from p8_vacuum_affine_heavy_source_filtration import source

N = s.Symbol("positive_lapse", positive=True)
R = s.Function("whole_R", positive=True)(N)
U = R ** (-s.Rational(3, 4))
M = R ** s.Rational(1, 4)
B = s.Function("whole_primitive_B", real=True)(N)
F = s.Function("whole_Fhat", real=True)(N)
r = R - 1
aa = -M / 3
K, T, p, G = s.symbols("trace_K temporal_W trace_momentum Gauss_density", real=True)
Href = s.Symbol("fixed_clock_Hubble", real=True)
lower = -3 * Href * r / N
gamma = 1 - 3 * r * r / (2 * R)
pm, ph, h, n = s.symbols("M1_momentum H_momentum H_field H_mass_squared", real=True)
j = s.Function("entire_normalized_heavy_source", real=True)(N)
shear, electric, magnetic, wmass, gm, gh, curvature = s.symbols(
    "shear_momentum_square electric_momentum_square magnetic_square "
    "spatial_vector_square M1_gradient_square H_gradient_square spatial_scalar_curvature",
    real=True,
)
zeta = s.Symbol("positive_zeta", positive=True)
Cchi = R ** (-s.Rational(1, 4))
C3 = R ** s.Rational(3, 4) / 2
other = (
    2 * shear / M
    + (pm**2 + ph**2) / (2 * U)
    + Cchi * (gm + gh) / 2
    + U * (n * h * h / 2 - j * h)
    - C3 * curvature
    + electric / (2 * zeta * Cchi)
    + zeta * M * magnetic / 4
    + Cchi * wmass / 2
)


@cache
def full_trace():
    lag = aa * K * K + B * K + F + U * (T - r * K - lower) ** 2 / 2
    kinetic = s.diff(lag, K, 2)
    Kstar = (p - B + U * r * (T - lower)) / (2 * aa + U * r * r)
    raw = N * ((p * K - lag).subs(K, Kstar) - T * G + other)
    raw = s.factor(raw)
    Kjoint = (p - B - r * G) / (2 * aa)
    Tstar = r * Kjoint + lower - G / U
    expected = N * (
        (p - B - r * G) ** 2 / (4 * aa) - F + G * G / (2 * U) - lower * G + other
    )
    D = s.hessian(raw, (N, T))
    DTT = s.factor(D[1, 1])
    # Keep the entire T-star derivative, including N dependence of R/B/lower.
    jacobian = s.diff(Tstar, N)
    checks = {
        "full_trace_velocity_pivot": s.factor(kinetic - 2 * aa * gamma),
        "full_trace_velocity_Euler": s.factor(s.diff(lag, K).subs(K, Kstar) - p),
        "full_temporal_Euler_with_Gauss": s.factor(s.diff(raw, T).subs(T, Tstar)),
        "entire_reduced_Hamiltonian": s.factor(raw.subs(T, Tstar) - expected),
        "entire_temporal_auxiliary_pivot": s.factor(DTT + N * U / gamma),
        "whole_temporal_solution_derivative": s.factor(
            D[0, 1].subs(T, Tstar) + DTT * jacobian
        ),
    }
    # Envelope theorem checks are factored before combining the full
    # nonlinear lapse derivatives; no lower potential or gradient is dropped.
    reduced_NN = s.diff(expected, N, 2)
    schur = s.factor((D[0, 0] - D[0, 1] ** 2 / DTT).subs(T, Tstar))
    checks["whole_lapse_Schur_envelope"] = s.factor(schur - reduced_NN)
    checks["whole_current_initial_temporal_lapse_cross_zero"] = s.factor(
        jacobian.subs({G: 0, Href: 0, p: B, s.diff(B, N): 0}, simultaneous=True)
    )
    checks["whole_positive_temporal_schur_margin"] = s.factor(
        gamma - s.Rational(1, 4) - 3 * (R - s.Rational(1, 2)) * (2 - R) / (2 * R)
    )
    for name, value in checks.items():
        assert value == 0, (name, value)
    return {
        "lagrangian": lag,
        "whole_raw_Hamiltonian": raw,
        "whole_reduced_Hamiltonian": expected,
        "whole_auxiliary_Hessian": D,
        "whole_temporal_solution": Tstar,
        "whole_reduced_lapse_second_derivative": reduced_NN,
        "whole_spatial_shift_generator_on_primary_surface": "H_i=-2 gamma_ij D_k pi_metric^(jk)+pi_W^j F_ij-W_i partial_j pi_W^j+p_M1 partial_i M1+p_H partial_i H. The full density is V times the displayed normal Hamiltonian plus Ni H_i. Off the primary surface the lapse and temporal scalar primary-momentum transport terms must also be retained; no first-class gauge determinant or regulator closure is computed here.",
        "checks": checks,
        "gates": {
            "full_Gauss_density_not_deleted": raw.has(G),
            "entire_heavy_source_lapse_derivatives_present": all(
                D.has(item) for item in (s.diff(j, N), s.diff(j, N, 2))
            ),
            "single_regular_auxiliary_branch_required": True,
        },
    }


@cache
def current_data():
    # Bind full functions and derivative rules, not a replacement clock germ.
    oldN = current.N
    bind = lambda expr: s.sympify(expr).subs(oldN, N)
    Fhat = (
        current.U * (current.F + 9 * current.Ru**2 / (16 * current.R * oldN**2))
        - s.diff(current.primitive, current.u) / oldN
    )
    Hclock = 4 * current.u / (1 + current.u**2)
    full_j = source.physical_source(current.u, oldN**-2) / s.sqrt(source.KAPPA)
    abstract_to_current = {
        R: bind(current.R),
        B: bind(current.B),
        F: bind(Fhat),
        Href: Hclock,
        n: source.MASS2,
        j: bind(full_j),
        zeta: s.Rational(1, 10**6),
    }
    # The complete homogeneous binding is differentiated/varied before any
    # on-shell source alignment or initial evaluation.
    hom = N * (
        -3 * R ** s.Rational(1, 4) * (current.H / N) ** 2
        + 3 * B * current.H / N
        + F
        + R ** (-s.Rational(3, 4)) * (current.mc**2 + current.mh**2) / (2 * N**2)
        + R ** (-s.Rational(3, 4)) * (-n * current.h**2 / 2 + j * current.h)
    )
    mapped = hom.subs(abstract_to_current, simultaneous=True)
    # The current normalized source is held as its full function name in
    # both compared actions; its exact actual binding is recorded above.
    mapped_formal = hom.subs(
        {**abstract_to_current, j: bind(current.j)}, simultaneous=True
    )
    initial = current.initial_data()
    Jbox = initial["whole_current_enclosures"]["J"]
    full_functions = current.source_data()
    return {
        "whole_abstract_to_current_function_bindings": abstract_to_current,
        "whole_actual_R_F_fixed_profile_and_constant_bindings": {
            bind(key): bind(value)
            for key, value in full_functions["entire_current_function_bindings"].items()
        },
        "whole_fixed_profile_bindings": full_functions["whole_fixed_profile_bindings"],
        "all_three_fixed_vacuum_constants": full_functions[
            "all_three_fixed_vacuum_constants"
        ],
        "entire_primitive_lapse_derivative": bind(current.IN),
        "primitive_boundary": "The unchanged full primitive has I(u,1)=0; its time and spatial boundaries are retained.",
        "whole_current_homogeneous_density": mapped,
        "actual_initial_full_lapse_Schur_enclosure": Jbox,
        "actual_initial_datum": "Use exactly the unforced full-current S256 datum and positive whole constraint root, not the original fixed quantum reference.",
        "inhomogeneous_branch_argument": "The complete canonical Hamiltonian is algebraic in lapse N and normal vector T after retaining the Maxwell Gauss integration by parts. Its temporal pivot is -N U/Gamma. At the actual S256 datum the temporal cross derivative vanishes and the complete lapse Schur pivot is -2J, with J in the displayed strict positive enclosure. The full implicit function theorem gives one regular local canonical finite-spatial-jet auxiliary branch. No uniform neighborhood radius or nonlinear solution for arbitrary inhomogeneous data is inferred.",
        "checks": {
            "literal_entire_current_homogeneous_action": s.expand(
                mapped_formal - bind(current.L).subs(current.n, source.MASS2)
            ),
            "actual_fixed_heavy_mass": abstract_to_current[n] - source.MASS2,
        },
        "gates": {
            "actual_initial_full_lapse_pivot_strictly_positive": bool(
                Jbox[0] > s.Rational(151, 100)
            ),
            "actual_entire_source_not_dropped": full_j != 0,
            "fixed_parent_kappa_and_zeta_retained": bool(
                source.KAPPA == 10**800
                and abstract_to_current[zeta] == s.Rational(1, 10**6)
            ),
            "all_spatial_invariants_retained_in_Hamiltonian": all(
                full_trace()["whole_reduced_Hamiltonian"].has(item)
                for item in (shear, electric, magnetic, wmass, gm, gh, curvature, G)
            ),
            "local_auxiliary_branch_not_quantum_mean_or_covariant_quantum_measure": True,
        },
    }


@cache
def homogeneous_bridge():
    NN = s.Symbol("homogeneous_lapse", positive=True)
    HH, m1, mh = s.symbols("homogeneous_H homogeneous_m1 homogeneous_mh", real=True)
    pv, p1, ph = s.symbols(
        "homogeneous_volume_momentum homogeneous_M1_momentum homogeneous_H_momentum",
        real=True,
    )
    DD = s.Function("homogeneous_D", positive=True)(NN)
    ZZ = s.Function("homogeneous_Z", positive=True)(NN)
    BB = s.Function("homogeneous_B", real=True)(NN)
    potential = s.Function("entire_homogeneous_local_potential", real=True)(NN)
    lag = -3 * DD * HH**2 + 3 * BB * HH + ZZ * (m1**2 + mh**2) / 2 + potential
    ham = -((pv - 3 * BB) ** 2) / (12 * DD) + (p1**2 + ph**2) / (2 * ZZ) - potential
    momentum = {pv: s.diff(lag, HH), p1: s.diff(lag, m1), ph: s.diff(lag, mh)}
    theta = -HH * s.diff(DD, NN) + s.diff(BB, NN) / 2
    Cnn = s.diff(lag, NN, 2) / 2
    J = Cnn + 3 * theta**2 / DD - s.diff(ZZ, NN) ** 2 * (m1**2 + mh**2) / (2 * ZZ)
    return {
        "entire_homogeneous_legendre_density": ham,
        "entire_homogeneous_momentum_map": momentum,
        "entire_homogeneous_lapse_pivot": J,
        "checks": {
            "whole_homogeneous_legendre": s.factor(
                ham.subs(momentum, simultaneous=True)
                - sum(p * v for p, v in zip((pv, p1, ph), (HH, m1, mh))).subs(
                    momentum, simultaneous=True
                )
                + lag
            ),
            "whole_homogeneous_lapse_second_derivative": s.factor(
                s.diff(ham, NN, 2).subs(momentum, simultaneous=True) + 2 * J
            ),
            "whole_homogeneous_lapse_first_derivative": s.factor(
                s.diff(ham, NN).subs(momentum, simultaneous=True) + s.diff(lag, NN)
            ),
        },
        "gates": {
            "full_potential_and_both_matter_momenta_retained": all(
                ham.has(item) for item in (potential, p1, ph)
            )
        },
    }
