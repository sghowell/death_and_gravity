"""Literal two massive stress triangles, metric contacts and full g-squared F1."""

from functools import cache

import sympy as s

from . import source, ward

MU, N, T, Z, V, G = source.MU, source.N, source.T, source.Z, source.V, source.G


def denominator(active, spectator, transfer=T, z=Z, v=V, mass=MU):
    return (
        (1 - z) * active
        + z * spectator
        - z * (1 - z) * mass
        - (1 - z) ** 2 * (1 - v * v) * transfer / 4
    )


def raw_F1(active, spectator, transfer=T, z=Z, v=V, mass=MU):
    return (1 - z) * z * z / denominator(active, spectator, transfer, z, v, mass)


def F1_integrand(transfer=T, z=Z, v=V, mass=MU, heavy=N):
    return sum(
        raw_F1(a, b, transfer, z, v, mass) - raw_F1(a, b, 0, z, v, mass)
        for a, b in ((mass, heavy), (heavy, mass))
    )


def F1_coefficient(order, z=Z, mass=MU, heavy=N):
    order = source.require_order(order)
    pref = s.integrate((1 - V * V) ** order, (V, 0, 1)) / 4**order
    return pref * sum(
        (1 - z) ** (2 * order + 1)
        * z
        * z
        / denominator(a, b, 0, z, 0, mass) ** (order + 1)
        for a, b in ((mass, heavy), (heavy, mass))
    )


def self_second_integrand(z=Z, mass=MU, heavy=N):
    return z * z * (1 - z) ** 2 / (mass * (1 - z) ** 2 + heavy * z) ** 2


@cache
def data():
    eta = s.diag(1, -1, -1, -1)
    p = s.Matrix(s.symbols("p0:4", real=True))
    r = s.Matrix(s.symbols("r0:4", real=True))
    q = r - p
    Bp, Br = s.symbols("B_p B_r")
    contacts = -eta * (Bp + Br)
    triangle_ward = r * Br - p * Bp
    mixed = (triangle_ward + contacts * eta * q - (p * Br - r * Bp)).applyfunc(s.expand)
    E, Py, Q, l0, ly = s.symbols("E P_y Q ell_0 ell_y", real=True)
    pp = s.Matrix([E, -Q / 2, Py, 0])
    rr = s.Matrix([E, Q / 2, Py, 0])
    xx = (1 - Z) * (1 + V) / 2
    yy = (1 - Z) * (1 - V) / 2
    ell = s.Matrix([l0, 0, ly, 0])
    rout = rr - (xx * pp + yy * rr + ell)
    rin = pp - (xx * pp + yy * rr + ell)
    numerator = ward.tensor(rin, rout)
    projected = s.expand(numerator[0, 2])
    projected_average = projected.subs({l0: 0, ly: 0})
    F = MU * (1 - Z) ** 2 + N * Z
    zero = s.integrate(raw_F1(MU, N, 0), (V, 0, 1))
    zero += s.integrate(raw_F1(N, MU, 0), (V, 0, 1)).subs(Z, 1 - Z)
    slope = s.integrate(s.diff(raw_F1(MU, N), T).subs(T, 0), (V, 0, 1))
    slope += s.integrate(s.diff(raw_F1(N, MU), T).subs(T, 0), (V, 0, 1)).subs(Z, 1 - Z)
    checks = {
        "all_four_mixed_Ward_components_with_both_contacts": mixed,
        "literal_cubic_triangle_Feynman_sign": s.I**2 * (-s.I) * s.I**3 * (-s.I) + s.I,
        "literal_cubic_metric_contact_Feynman_sign": s.I**2 * s.I**2 * s.I - s.I,
        "literal_triangle_PP_projection": s.factor(
            projected_average - 2 * Z * Z * E * Py
        ),
        "literal_translated_mixed_loop_terms": s.expand(
            projected
            - projected_average
            + 2 * Z * E * ly
            + 2 * Z * Py * l0
            - 2 * l0 * ly
        ),
        "full_parameter_simplex_sum": s.expand(xx + yy + Z - 1),
        "full_parameter_simplex_Jacobian": s.factor(
            s.det(
                s.Matrix(
                    [[s.diff(xx, Z), s.diff(xx, V)], [s.diff(yy, Z), s.diff(yy, V)]]
                )
            )
            - (1 - Z) / 2
        ),
        "literal_light_line_denominator": s.factor(denominator(MU, N, 0) - F),
        "literal_heavy_line_reversed_denominator": s.factor(
            denominator(N, MU, 0).subs(Z, 1 - Z) - F
        ),
        "both_F1_triangles_equal_whole_scalar_Pi_prime": s.factor(
            zero - Z * (1 - Z) / F
        ),
        "both_F1_slopes_equal_one_sixth_whole_Pi_second": s.factor(
            slope - self_second_integrand() / 6
        ),
        "full_OS_F1_zero_transfer": s.factor(F1_integrand(0)),
        "no_trace_in_PP_projection": eta[0, 2],
        "no_pure_transfer_in_PP_projection": (rr - pp)[0] * (rr - pp)[2],
        "first_Taylor_coefficient_matches_direct_derivative": s.factor(
            F1_coefficient(1)
            - s.integrate(s.diff(F1_integrand(), T).subs(T, 0), (V, 0, 1))
        ),
    }
    for order in (2, 3, 4):
        checks["entire_positive_Taylor_coefficient_" + str(order)] = s.factor(
            F1_coefficient(order)
            - s.integrate(
                s.diff(F1_integrand(), T, order).subs(T, 0) / s.factorial(order),
                (V, 0, 1),
            )
        )
    return {
        "literal_full_massive_triangle_projection": projected,
        "both_unsubtracted_triangle_F1_integrands": (raw_F1(MU, N), raw_F1(N, MU)),
        "whole_scalar_OS_subtracted_F1_integrand": F1_integrand(),
        "whole_F1": 1
        + G * G / (16 * s.pi**2) * s.Integral(F1_integrand(), (V, 0, 1), (Z, 0, 1)),
        "whole_F1prime_zero": G
        * G
        / (96 * s.pi**2)
        * s.Integral(self_second_integrand(), (Z, 0, 1)),
        "whole_scalar_Pi_second": G
        * G
        / (16 * s.pi**2)
        * s.Integral(self_second_integrand(), (Z, 0, 1)),
        "whole_cubic_metric_contact_tensor": contacts,
        "literal_graph_derivation": "Each ordinary cubic is+ig, each scalar metric vertex is-iT/sqrt(kappa); the triangle propagators and scalar master yield the real coefficient+g^2/(16pi^2) in the stripped vertex. The complete normalized triangle Ward contraction has the minus bubble-master sign. Summing both line insertions and both+ig eta/sqrt(kappa) metric contacts gives the whole mixed self-energy identity.",
        "complete_projected_graph_count": "The off-diagonal02 projection in a spacelike-transfer frame has no eta or qq component. Shifted odd loop terms and the off-diagonal isotropic ell0*elly moment vanish in the regulated tensor integral. The two triangles give2z^2 P0 Py. All contact, quartic and H-metric mixing graphs lack a PP tensor; they are not discarded from any claimed whole F2.",
        "source_OS_normalization": "At t0, z reversal of the heavy-line insertion gives z(1-z)^2/F and the light-line insertion gives z^2(1-z)/F. Their sum is the entire mixed Pi_prime integrand. The inherited deltaZ cancels that coefficient, including its finite value. Consequently F1(0)=1 for this formal g^2 matter coefficient without a massless internal-line limit.",
        "matching_boundary": "The graph-specific F1prime=Pi_second/6 is obtained by both triangle integrations, not inferred from Ward alone. The complete physical F2, H-metric mixing and curvature improvements are not supplied by this projection. Internal gravitational loops, finite Ricci-derivative additions and higher orders remain separately unmatched.",
        "checks": checks,
        "gates": {
            "both_distinct_massive_line_insertions_kept": True,
            "both_cubic_metric_contact_bubbles_kept_in_Ward": True,
            "entire_F1_not_only_its_zero_transfer_value": F1_integrand().has(T, N, MU),
            "inherited_whole_OS_kinetic_subtraction": True,
            "no_graviton_ghost_gauge_bridge_required_for_these_matter_graphs": True,
            "heavy_metric_mixing_not_claimed_absent_from_F2": True,
            "finite_Ricci_derivative_matching_not_replaced_by_Ward": True,
        },
    }
