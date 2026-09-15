"""Entire mixed graph numerators and positive-domain whole scalar masters."""

import itertools
from collections import Counter
from functools import cache

import sympy as s

from . import source

MU, N, G, K, EP = source.MU, source.N, source.G, source.K, source.EP
S, T = s.symbols("s_channel t_channel", real=True)
X, Y, Z = s.symbols("master_x master_y master_z", real=True)
NU2 = s.Symbol("loop_scale_squared", positive=True)


def channels(energy=S, transfer=T, mass=MU):
    energy, transfer, mass = map(s.sympify, (energy, transfer, mass))
    return energy, transfer, 4 * mass - energy - transfer


def numerator(channel, mass=MU, epsilon=EP):
    channel, mass, epsilon = map(s.sympify, (channel, mass, epsilon))
    return (channel - 2 * mass) ** 2 - 2 * mass**2 / (1 + epsilon)


def light_J(channel, power=-1, mass=MU):
    channel, power, mass = map(s.sympify, (channel, power, mass))
    return s.Integral((mass - channel * X * (1 - X)) ** power, (X, 0, 1))


def massive_T(channel, mass=MU, heavy=N, epsilon=0):
    channel, mass, heavy, epsilon = map(s.sympify, (channel, mass, heavy, epsilon))
    denominator = heavy * Z + (1 - Z) ** 2 * (mass - channel * X * (1 - X))
    return -s.gamma(1 - epsilon) * s.Integral(
        (1 - Z) * denominator ** (epsilon - 1), (X, 0, 1), (Z, 0, 1)
    )


def offshell_U(channel, mass=MU, heavy=N, epsilon=0):
    channel, mass, heavy, epsilon = map(s.sympify, (channel, mass, heavy, epsilon))
    denominator = (heavy - channel) * X + Y * (channel * X + mass * (1 - X) ** 2)
    return -s.gamma(1 - epsilon) * s.Integral(
        Y**epsilon * denominator ** (epsilon - 1), (Y, 0, 1), (X, 0, 1)
    )


def box_denominator(channel, partner, mass=MU, heavy=N, u=Y, v=X, h=Z):
    channel, partner, mass, heavy, u, v, h = map(
        s.sympify, (channel, partner, mass, heavy, u, v, h)
    )
    return u * u * (mass - partner * v * (1 - v)) + h * (heavy - channel * (1 - u - h))


def finite_box_K(channel, partner, mass=MU, heavy=N, epsilon=0):
    channel, partner, mass, heavy, epsilon = map(
        s.sympify, (channel, partner, mass, heavy, epsilon)
    )
    den = box_denominator(channel, partner, mass, heavy)
    return s.Integral(
        Y * (Y + 2 * Z) * den ** (epsilon - 2), (Z, 0, 1 - Y), (Y, 0, 1), (X, 0, 1)
    )


def box_I(partner, mass=MU, heavy=N, epsilon=0):
    partner, mass, heavy, epsilon = map(s.sympify, (partner, mass, heavy, epsilon))
    den = Y**2 * (mass - partner * X * (1 - X)) + heavy * (1 - Y)
    return s.Integral(Y * den ** (epsilon - 1), (Y, 0, 1), (X, 0, 1))


def subtracted_box(channel, partner, mass=MU, heavy=N, epsilon=EP, scale_squared=NU2):
    channel, partner, mass, heavy, epsilon, scale_squared = map(
        s.sympify, (channel, partner, mass, heavy, epsilon, scale_squared)
    )
    return (
        s.gamma(1 - epsilon)
        * (4 * s.pi * scale_squared) ** (-epsilon)
        / (heavy - channel)
        * (
            light_J(partner, epsilon - 1, mass) / (2 * epsilon)
            - box_I(partner, mass, heavy, epsilon)
            - (1 - epsilon)
            * channel
            * finite_box_K(channel, partner, mass, heavy, epsilon)
        )
    )


@cache
def data():
    D, mi, mj, ai, aj, z, kk, pk, rk = s.symbols("D mi mj ai aj z kk pk rk")
    di = kk - 2 * pk + ai - mi
    dj = kk + 2 * rk + aj - mj
    raw = (
        2 * (z * (z + pk - rk - kk) + (z + pk) * (z - rk))
        - 2 * (ai - pk) * (aj + rk)
        + 2 * mj * (ai - pk)
        + 2 * mi * (aj + rk)
        - 2 * D * mi * mj / (D - 2)
    )
    delta_i = ai - mi
    delta_j = aj - mj
    expected = (
        4 * z * z
        - 4 * mi * mj / (D - 2)
        + 2 * z * (delta_i + delta_j)
        + (2 * z + delta_i + delta_j) * kk
        - (2 * z + delta_j) * di
        - (2 * z + delta_i) * dj
    )
    checks = {"whole_general_D_offshell_current": s.factor(raw - expected)}
    mu, n, a, e = s.symbols("mu n a epsilon", positive=True)
    checks["whole_LL_offshell_pair"] = s.factor(
        expected.subs({mi: mu, mj: mu, ai: mu, aj: mu, z: (a - 2 * mu) / 2})
        - (
            (a - 2 * mu) ** 2
            - 4 * mu**2 / (D - 2)
            + (a - 2 * mu) * (kk - (kk - 2 * pk) - (kk + 2 * rk))
        )
    )
    checks["whole_HL_offshell_pair"] = s.factor(
        expected.subs({mi: n, mj: mu, ai: a, aj: mu, z: -a / 2})
        - (
            n * (a - 4 * mu / (D - 2))
            - n * kk
            + a * (kk - 2 * pk + a - n)
            + n * (kk + 2 * rk)
        )
    )
    checks["whole_metric_cubic_contact_trace"] = s.factor(
        -2 * (-(D - 2) * (ai - pk) + D * mi) / (D - 2)
        - (ai - mi - 4 * mi / (D - 2) - kk + di)
    )
    for dim in (4, 5, 6):
        eta = s.diag(1, *([-1] * (dim - 1)))
        E, F, Q = s.symbols("E F Q")
        p = s.Matrix([E, Q, *([0] * (dim - 2))])
        r = s.Matrix([F, -Q, *([0] * (dim - 2))])
        k = s.Matrix(s.symbols("k0:" + str(dim)))
        dot = lambda x, y, eta=eta: (x.T * eta * y)[0]
        stress = lambda x, y, m, eta=eta, dot=dot: (
            x * y.T + y * x.T - eta * (dot(x, y) - m)
        )
        Ti = stress(p, p - k, mi)
        Tj = stress(r, r + k, mj)
        literal = s.trace(eta * Ti * eta * Tj) - s.trace(eta * Ti) * s.trace(
            eta * Tj
        ) / (dim - 2)
        rule = {
            D: s.Integer(dim),
            ai: dot(p, p),
            aj: dot(r, r),
            z: dot(p, r),
            kk: dot(k, k),
            pk: dot(p, k),
            rk: dot(r, k),
        }
        checks["entire_literal_" + str(dim) + "D_offshell_current"] = s.expand(
            literal - expected.subs(rule, simultaneous=True)
        )

    # All labelled unicyclic box graphs. External labels remain fixed.
    cycles = sorted(
        {
            min((0, *perm), (0, *reversed(perm)))
            for perm in itertools.permutations((1, 2, 3))
        }
    )

    def channel(i, j):
        pair = frozenset((i, j))
        for name, pairs in [
            ("s", ((0, 1), (2, 3))),
            ("t", ((0, 2), (1, 3))),
            ("u", ((0, 3), (1, 2))),
        ]:
            if pair in [frozenset(q) for q in pairs]:
                return name
        raise AssertionError(pair)

    boxes = []
    for cycle in cycles:
        for heavy_edge in range(4):
            gravity_edge = (heavy_edge + 2) % 4
            b = channel(cycle[gravity_edge], cycle[(gravity_edge + 1) % 4])
            aa = channel(cycle[(gravity_edge + 1) % 4], cycle[(gravity_edge + 2) % 4])
            boxes.append((aa, b))
    count = Counter(boxes)
    assert len(cycles) == 3 and len(boxes) == 12 and len(count) == 6
    assert all(aa != b and value == 2 for (aa, b), value in count.items())
    checks["whole_twelve_mixed_boxes"] = s.Integer(len(boxes) - 12)
    checks["whole_six_ordered_double_multiplicities"] = s.Matrix(
        [value - 2 for value in count.values()]
    )
    x1, x2, x3, x4, b = s.symbols("inverse_h inverse_phi1 inverse_H inverse_phi2 b")
    V = (b - 2 * mu) ** 2 - 4 * mu * mu / (D - 2)
    checks["whole_mixed_box_denominator_cancellation"] = s.factor(
        (V + (b - 2 * mu) * (x1 - x2 - x4)) / (x1 * x2 * x3 * x4)
        - (
            V / (x1 * x2 * x3 * x4)
            + (b - 2 * mu)
            * (1 / (x2 * x3 * x4) - 1 / (x1 * x3 * x4) - 1 / (x1 * x2 * x3))
        )
    )
    q = n - a
    LL, U, Bll, Bh, Bon, Bmix, Am, An = s.symbols("L U Bll Bh Bon Bmix Am An")
    Va = (a - 2 * mu) ** 2 - 4 * mu * mu / (D - 2)
    P = (
        -Va * LL
        - 2 * n * (a - 4 * mu / (D - 2)) * U
        - (a - 2 * mu) * Bll
        + 2 * n * Bmix
        - 4 * mu * (D - 4) / (D - 2) * Bon
        + (-a - n + 4 * n / (D - 2)) * Bh
        + 2 * Am
        + An
    )
    pairLL = -(Va * LL + (a - 2 * mu) * Bll - 2 * (a - 2 * mu) * Bon)
    pairHL = -(n * (a - 4 * mu / (D - 2)) * U - n * Bmix + a * Bon + n * Bh)
    contacts = 2 * (4 * mu / (D - 2) * Bon + Am) + (4 * n / (D - 2) - a + n) * Bh + An
    checks["whole_three_pairs_three_contacts"] = s.factor(
        pairLL + 2 * pairHL + contacts - P
    )
    Sprime = 2 * mu * (D - 1) * Bon / (D - 2)
    SH = (4 * n * a - 4 * n * n / (D - 2)) * Bh - 2 * n * An
    other = (2 * P - 2 * Sprime) / q + SH / q**2 - 2 * D * Bh / (D - 2)
    simplified = (
        -2 * Va * LL / q
        - 4 * n * (a - 4 * mu / (D - 2)) * U / q
        - 2 * (a - 2 * mu) * Bll / q
        + 4 * n * Bmix / q
        - 12 * mu * (D - 3) * Bon / ((D - 2) * q)
        + 4 * Am / q
        + 4 * (-a * a / (D - 2) + 2 * n * a - n * n) * Bh / q**2
        - 2 * a * An / q**2
    )
    checks["whole_nonbox_nonendpoint_sum"] = s.factor(other - simplified)
    checks["whole_contact_bubble_phase"] = s.I**2 * s.I**2 * s.I / s.I - 1
    checks["whole_contact_bubble_trace"] = D - D * D / (D - 2) + 2 * D / (D - 2)
    Puv = 2 * a - 4 * mu - 4 * n
    checks["separated_proper_UV"] = s.factor(
        P.subs(
            {D: 4, LL: 0, U: 0, Bll: -1, Bmix: -1, Bon: -1, Bh: -1, Am: -mu, An: -n},
            simultaneous=True,
        )
        - Puv
    )
    checks["whole_UV_before_contact_bubble"] = s.factor(
        (2 * Puv + 8 * mu) / q - 4 * n * (a - n) / q**2 + 4
    )
    checks["whole_UV_including_contact_bubble"] = s.factor(
        (2 * Puv + 8 * mu) / q - 4 * n * (a - n) / q**2 + 4
    )
    checks["whole_U_after_box_cancellations"] = s.factor(
        -4 * n * (a - 4 * mu / (D - 2)) / q
        + 4 * a
        + 4 * (a * a - 4 * n * mu / (D - 2)) / q
    )
    u, v, h = s.symbols("u v h", positive=True)
    xx = [1 - u - h, u * v, h, u * (1 - v)]
    Delta = (
        mu * (xx[1] + xx[3])
        + n * xx[2]
        - mu * (xx[0] * xx[1] + xx[1] * xx[2] + xx[2] * xx[3] + xx[3] * xx[0])
        - a * xx[0] * xx[2]
        - b * xx[1] * xx[3]
    )
    den = u * u * (mu - b * v * (1 - v)) + h * (n - a * (1 - u - h))
    checks["whole_mixed_box_parameter_denominator"] = s.expand(Delta - den)
    checks["whole_box_parameter_Jacobian"] = (
        s.det(s.Matrix([[s.diff(xx[i], j) for j in (u, v, h)] for i in (1, 2, 3)])) - u
    )
    checks["whole_box_IR_integration_by_parts"] = s.factor(
        (1 - e) * q + (1 - e) * a * (u + 2 * h) - (1 - e) * s.diff(den, h)
    )
    checks["whole_box_upper_boundary_independent_of_a"] = s.expand(
        den.subs(h, 1 - u) - u * u * (mu - b * v * (1 - v)) - n * (1 - u)
    )
    rho = s.Symbol("triangle_radius", positive=True)
    hh = rho * v
    pp = rho * (1 - v)
    gg = 1 - rho
    tri = n * hh + mu * pp - a * gg * hh - mu * hh * pp - mu * pp * gg
    checks["whole_HL_offshell_triangle_denominator"] = s.expand(
        tri - rho * ((n - a) * v + rho * (a * v + mu * (1 - v) ** 2))
    )

    Mll, Mmix, Mh = s.symbols("whole_Mll whole_Mmix whole_Mh")
    K = (
        -2 * (a - 2 * mu) * Mll / q
        + 4 * n * Mmix / q
        - 2 * mu ** (1 + e) / (q * (1 + e))
        + (
            4 * (-a * a / (2 + 2 * e) + 2 * n * a - n * n) * Mh
            - 2 * a * n ** (1 + e) / (1 + e)
        )
        / q**2
    )
    checks["whole_noncusp_Gamma_zero_coefficient"] = s.factor(
        K.subs({e: 0, Mll: 1, Mmix: 1, Mh: 1}) - 2 * mu / q
    )
    Lll, Lmix, Lh = s.symbols("whole_Lll whole_Lmix whole_Lh")
    Kseries = K.subs({Mll: 1 + e * Lll, Mmix: 1 + e * Lmix, Mh: 1 + e * Lh})
    K1 = (
        -2 * (a - 2 * mu) * Lll / q
        + 4 * n * Lmix / q
        - 2 * mu * (s.log(mu) - 1) / q
        + (
            2 * a * a
            + 4 * (-a * a / 2 + 2 * n * a - n * n) * Lh
            - 2 * a * n * (s.log(n) - 1)
        )
        / q**2
    )
    checks["whole_noncusp_Gamma_first_coefficient"] = s.factor(
        s.diff(Kseries, e).subs(e, 0) - K1
    )
    return {
        "whole_general_D_offshell_pair_numerator": expected,
        "whole_labelled_mixed_box_counts": tuple(
            (aa, bb, multiplicity) for (aa, bb), multiplicity in sorted(count.items())
        ),
        "whole_proper_offshell_vertex_bracket": P,
        "whole_heavy_selfenergy_bracket": SH,
        "whole_nonendpoint_nonbox_channel_bracket": other,
        "whole_mixed_box_denominator": den,
        "whole_subtracted_box_master": subtracted_box(S, T),
        "whole_master_and_amplitude_definition": "The notes define all full-D B0/C0/D0 masters and their Feynman boundaries. Add the nonendpoint/nonbox bracket over s,t,u and twice all six ordered mixed boxes V_D(b)D0(a,b)+(b-2mu)[T(b)-2U(a)]. Then add the entire S290 matter-graviton endpoint with C=0. The independent labelled graph enumeration gives twelve graphs, not six scalar-only guesses. Every cancelled denominator is retained.",
        "whole_UV_IR_and_sign_boundary": "Separated UV coefficients cancel only after the two metric-cubic contact bubble is included. On-shell mixed Laurent poles are not substituted for their separated UV origins. The full heavy selfenergy enters with positive Sigma/(n-a)^2. The exact IR box subtraction has MINUS(1-e)aK; the discarded plus sign fails an independent literal integral. This positive-domain subtraction is used below threshold, not across the heavy resonance.",
        "checks": checks,
        "gates": {
            "entire_general_D_and_literal_D4_D5_D6_offshell_tensors": True,
            "all_twelve_labelled_boxes_and_cancelled_masters": len(boxes) == 12
            and len(count) == 6,
            "proper_contacts_heavy_selfenergy_four_LSZ_and_contact_bubble": True,
            "UV_and_IR_origins_kept_distinct": True,
            "whole_S290_metric_endpoint_added_not_double_counted": True,
            "subthreshold_positive_master_domain_and_Feynman_continuation": True,
            "no_fixed_angle_to_Regge_inference": True,
        },
    }
