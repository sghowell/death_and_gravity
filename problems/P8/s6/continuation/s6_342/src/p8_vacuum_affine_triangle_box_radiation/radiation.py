"""Literal mass-ordered TT loop densities; couplings and loop factor excluded."""

from functools import cache
from math import factorial

import sympy as s
from p8_vacuum_affine_physical_loop_contours import contour, parameters

ETA = s.diag(1, -1, -1, -1)


def dot(p, q):
    return (p.T * ETA * q)[0]


def setup(vertices, masses, k, split, gamma, heavy_mass):
    N = len(vertices)
    if N not in (3, 4) or len(masses) != N:
        raise ValueError("Require an original triangle or alternating box")
    masses = tuple(parameters.rational(value) for value in masses)
    if isinstance(split, bool) or not isinstance(split, int) or not 0 <= split < N:
        raise ValueError("Require a stated split line")
    gamma = parameters.rational(gamma)
    heavy_mass = parameters.heavy_mass(heavy_mass)
    if not 0 <= gamma <= 1:
        raise ValueError("Require the full physical split interval")
    if sum(vertices, s.zeros(4, 1)) + k != s.zeros(4, 1) or dot(k, k) != 0:
        raise ValueError(
            "Require exact radiative conservation and a null external graviton"
        )
    r = list(vertices[split:]) + list(vertices[:split])
    masses = list(masses[split:]) + list(masses[:split])
    if sorted(masses) != (
        [1, 1, heavy_mass] if N == 3 else [1, 1, heavy_mass, heavy_mass]
    ):
        raise ValueError("Require the stated original massive scalar assignment")
    qs = [s.zeros(4, 1)]
    for p in r[:-1]:
        qs.append(qs[-1] - p)
    eff = [q + gamma * k if j == 0 else q for j, q in enumerate(qs)]
    lights = [j for j, mass in enumerate(masses) if mass == 1]
    if N == 3:
        heavy = next(j for j, mass in enumerate(masses) if mass == heavy_mass)
        order = (lights[0], lights[1], heavy)
        v = dot(eff[order[0]] - eff[order[1]], eff[order[0]] - eff[order[1]])
        virtualities = tuple(
            dot(eff[j] - eff[heavy], eff[j] - eff[heavy]) for j in lights
        )
        w = s.S.Zero
        kind = "triangle"
    else:
        start = lights[0]
        order = tuple((start + j) % 4 for j in range(4))
        if any(
            masses[j] != (1 if pos % 2 == 0 else heavy_mass)
            for pos, j in enumerate(order)
        ):
            raise ValueError("Require the ordered alternating mass assignment")
        v = dot(eff[order[0]] - eff[order[2]], eff[order[0]] - eff[order[2]])
        w = dot(eff[order[1]] - eff[order[3]], eff[order[1]] - eff[order[3]])
        virtualities = tuple(
            dot(
                eff[order[j]] - eff[order[(j + 1) % 4]],
                eff[order[j]] - eff[order[(j + 1) % 4]],
            )
            for j in range(4)
        )
        kind = "ordered_box"
    parameters.domain(kind, v, virtualities, w, heavy_mass)
    return qs, eff, masses, order, v, w, virtualities, kind


def line_density(vertices, masses, k, eps, split, x, t, heavy_angle, gamma, heavy_mass):
    if (
        eps.shape != (4, 4)
        or eps != eps.T
        or s.trace(ETA * eps) != 0
        or eps * k != s.zeros(4, 1)
    ):
        raise ValueError("Require an exact physical TT polarization tensor")
    qs, eff, masses, order, v, w, aa, kind = setup(
        vertices, masses, k, split, gamma, heavy_mass
    )
    N = len(qs)
    for value in (x, t, heavy_angle):
        if not 0 <= parameters.rational(value) <= 1:
            raise ValueError("Require real unit-cube contour coordinates")
    xi = contour.light(x, timelike=bool(v > 4))
    z = contour.heavy(t)
    weights = [s.S.Zero] * N
    if N == 3:
        weights[order[0]] = (1 - z) * xi
        weights[order[1]] = (1 - z) * (1 - xi)
        weights[order[2]] = z
        e = xi * aa[0] + (1 - xi) * aa[1]
        b = s.S.Zero
        measure = 1 - z
    else:
        weights[order[0]] = (1 - z) * xi
        weights[order[2]] = (1 - z) * (1 - xi)
        weights[order[1]] = z * heavy_angle
        weights[order[3]] = z * (1 - heavy_angle)
        e = xi * (heavy_angle * aa[0] + (1 - heavy_angle) * aa[3]) + (1 - xi) * (
            heavy_angle * aa[1] + (1 - heavy_angle) * aa[2]
        )
        b = w * heavy_angle * (1 - heavy_angle)
        measure = z * (1 - z)
    bar = sum((a * q for a, q in zip(weights, qs)), s.zeros(4, 1))
    R = qs[0] - bar
    den = parameters.denominator(z, 1 - v * xi * (1 - xi), heavy_mass, 1 - e, b)
    numerator = -2 * factorial(N - 2) * weights[0] * (R.T * eps * R)[0]
    xjac = 1 + s.I * (1 - 6 * x * (1 - x)) if v > 4 else s.S.One
    zjac = 1 - s.I * (1 - 2 * t)
    density = measure * xjac * zjac * numerator / den ** (N - 1)
    return {
        "density": density,
        "denominator": den,
        "numerator": numerator,
        "weights": weights,
        "route": qs,
        "effective_route": eff,
        "masses": masses,
        "shift_vector": R,
        "order": order,
        "kind": kind,
        "light_invariant": v,
        "heavy_invariant": w,
        "virtualities": aa,
    }


@cache
def data():
    checks = {}
    gates = {}
    x, aa, v, u, Hpp = s.symbols("x a v u Hpp")
    M0 = 1 - v * x * (1 - x)
    M1 = 1 - u * x * (1 - x)
    y = s.Symbol("split_y")
    light0 = -2 * Hpp * x * x / (M0 - 2 * x * aa * y)
    anti0 = Hpp * x * s.log(M0 - 2 * x * aa * y) / aa
    checks["bubble_line0_literal_negative_kernel_antiderivative"] = s.factor(
        s.diff(anti0, y) - light0
    )
    endpoint0 = Hpp * x * (s.log(M1) - s.log(M0)) / aa
    endpoint1 = Hpp * (1 - x) * (s.log(M1) - s.log(M0)) / aa
    checks["both_bubble_lines_give_known_negative_DD"] = s.expand(
        endpoint0 + endpoint1 - Hpp * (s.log(M1) - s.log(M0)) / aa
    )
    checks["bubble_endpoint_denominator_routing"] = s.expand(
        (M0 - 2 * x * aa * (1 - x)).subs(aa, (u - v) / 2) - M1
    )
    gates["opposite_bubble_insertion_sign_fails"] = (
        s.expand(-(endpoint0 + endpoint1) - Hpp * (s.log(M1) - s.log(M0)) / aa) != 0
    )
    n = s.Symbol("n")
    checks["literal_heavy_resolvent_DD_sign"] = s.factor(
        -2 * Hpp * (1 / (n - u) - 1 / (n - v)) / (u - v) + 2 * Hpp / ((n - u) * (n - v))
    )
    J = s.symbols("J0:4")
    J0 = s.symbols("JB0:4")
    F = s.symbols("F0:4")
    Fb, Fsym, internal = s.symbols("Fborn Fsym internal")
    full = sum(j * (f - Fsym) for j, f in zip(J, F)) + internal
    leading = sum(J0) * (Fb - Fsym)
    remainder = (
        sum(j * (f - Fb) for j, f in zip(J, F))
        + (sum(J) - sum(J0)) * (Fb - Fsym)
        + internal
    )
    checks["whole_four_distinct_shifted_coefficients_soft_subtraction"] = s.expand(
        full - leading - remainder
    )

    return {
        "checks": {
            key: s.expand(value)
            if not value.is_rational_function()
            else s.cancel(value)
            for key, value in checks.items()
        },
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_two_separate_bubble_line_endpoints": (endpoint0, endpoint1),
        "whole_literal_insertion_sign": "For every internal scalar line the normalized D4 contour kernel is-2*Gamma(N-1)*xa*epsilon(Qa-barQ,Qa-barQ)/Delta_gamma^(N-1), with the complete original simplex and contour Jacobians. The sign is fixed by deltaG=-G deltaK G, the independent two-light bubble divided difference and literal S295 heavy exchange, not by a Ward-only transverse guess.",
        "whole_full_physical_TT_graph_sum": "M5=(sum_i J_i[F(p_i+k)-F_sym] +g^2/4 sum_perm[A(v_ab)*I_T(p_a+p_b,p_c,p_d) -2g^2 epsilon(P_ab,P_ab)*Cbar(u_ab;1,1)/((n-v_ab)(n-u_ab))] +g^4/4 sum_perm I_D(p_a,p_b,p_c,p_d))/(16pi^2 sqrt(kappa)). I_T and I_D integrate the sum of all literal line_density kernels over the complete unit cube, with no omitted scalar-line weights. u_ab=(P_ab+k)^2. The outer branch carries Cbar(u_ab), not Cbar(v_ab).",
        "whole_counterterm_and_UV_completion": "These triangle/box classes and their TT insertions are UV finite. Their fixed finite OS4 constant contributes-F_sym*sum J_i. The three old UV counterterms cancel the full factorized bubble class in S339 and are not double-counted here. Selected quadratic external attachments cancel as in S340 for arbitrary shifted hard coefficients. Pure-trace potential metric vertices vanish only after physical TT projection.",
        "whole_remainder_subtraction": remainder,
        "whole_soft_and_coincident_limits": "All internal contour denominators stay gapped, so the internal and outer-heavy kernels are finite as omega->0. The only leading1/omega term is the complete selected renormalized Born loop coefficient times the original soft current. Equal channel invariants use the ordinary gamma integral or continuous divided difference; no artificial pole is introduced.",
        "whole_boundary": "Complete selected minimal-matter triangle/box real-TT radiation and its nonleading remainder, not independent curved matching, internal-graviton loops, finite-gravity quantum decoupling or a full inclusive probability.",
    }
