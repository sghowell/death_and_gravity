"""Independent four-vector triangle and frozen full-kernel calibrations."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_minimal_gravity_radiation import vertices as original
from p8_vacuum_affine_radiative_curvature_matching import contact
from p8_vacuum_affine_triangle_box_radiation import radiation as frozen

from . import source

ETA = s.diag(1, -1, -1, -1)


def dot(p, q):
    return (p.T * ETA * q)[0]


def point_difference(vertices, k, eps, weights, order):
    gamma = s.Symbol("calibration_gamma")
    squares = tuple(dot(r, r) for r in vertices)
    soft = tuple(dot(k, r) for r in vertices)
    hs = tuple((r.T * eps * r)[0] for r in vertices)
    cs = (weights[0] * weights[1], weights[1] * weights[2], weights[0] * weights[2])
    U = sum(c * v for c, v in zip(cs, squares))
    lift = -2 * sum(
        h
        * sum(
            s.binomial(order, q) * U ** (order - q) * c**q * (2 * a) ** (q - 1)
            for q in range(1, order + 1)
        )
        for h, c, a in zip(hs, cs, soft)
    )
    total = 0
    for split in range(3):
        rr = vertices[split:] + vertices[:split]
        ws = weights[split:] + weights[:split]
        qs = (s.zeros(4, 1), -rr[0], -rr[0] - rr[1])
        effective = (gamma * k, qs[1], qs[2])
        bar = sum((w * q for w, q in zip(ws, qs)), s.zeros(4, 1))
        num = (bar.T * eps * bar)[0]
        denominator_jet = sum(
            ws[i]
            * ws[j]
            * dot(effective[i] - effective[j], effective[i] - effective[j])
            for i in range(3)
            for j in range(i + 1, 3)
        )
        terms = s.Poly(s.expand(denominator_jet ** (order - 1)), gamma)
        total -= (
            2
            * order
            * ws[0]
            * num
            * sum(co / (power[0] + 1) for power, co in terms.terms())
        )
    return s.factor(total - lift)


@cache
def data():
    checks, gates = {}, {}
    R = s.Rational
    weights = (R(1, 5), R(2, 7), R(18, 35))
    literal_count, kernel_count, nonzero = 0, 0, 0
    for state in range(3):
        ps, k, _ = original.sample(state)
        frame = tensor.frame(k)
        checks[f"original{state}_mass_shells"] = s.Matrix(
            [dot(p, p) - 1 for p in ps]
        ).applyfunc(s.factor)
        checks[f"original{state}_conservation"] = sum(ps, k.copy()).applyfunc(s.factor)
        for pol, (eps, norm) in enumerate(tensor.transverse_polarizations(frame, ETA)):
            curvature = contact.linear_curvature(k, eps)
            for pair, labels in enumerate(((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))):
                a, b, c, d = labels
                vertices = (ps[a] + ps[b], ps[c], ps[d])
                r1, r2 = vertices[:2]
                curvature_word = sum(
                    curvature[i, j, l, m] * r1[i] * r2[j] * r1[l] * r2[m]
                    for i, j, l, m in product(range(4), repeat=4)
                )
                edge = dot(k, r1) * r2 - dot(k, r2) * r1
                direct = (edge.T * eps * edge)[0]
                tag = f"state{state}_pol{pol}_pair{pair}"
                checks[tag + "_literal_curvature_tensor"] = s.factor(
                    curvature_word - direct
                )
                for order in (1, 2, 3):
                    difference = point_difference(vertices, k, eps, weights, order)
                    expected = (
                        0
                        if order < 3
                        else -8 * s.prod(w * w for w in weights) * curvature_word
                    )
                    checks[tag + f"_entire_vector_degree{2 * order}"] = s.factor(
                        difference - expected
                    )
                    literal_count += 1
                nonzero += curvature_word != 0
                # A common rescaling is only an OFF-SHELL calibration near the
                # analytic origin; the physical light mass in the source stays1.
                scale = R(1, 16)
                vv = tuple(scale * r for r in vertices)
                kk = scale * k
                z = R(1, 5) - s.I * R(4, 25)
                basew = ((1 - z) / 3, 2 * (1 - z) / 3, z)
                gamma = R(2, 5)
                for split in range(3):
                    rr = vv[split:] + vv[:split]
                    ws = basew[split:] + basew[:split]
                    masses = (1, 1, source.HEAVY_MASS2)
                    mm = masses[split:] + masses[:split]
                    first_light = next(i for i, m in enumerate(mm) if m == 1)
                    original_light = (split + first_light) % 3
                    angle = R(1, 3) if original_light == 0 else R(2, 3)
                    actual = frozen.line_density(
                        vv,
                        masses,
                        kk,
                        eps,
                        split,
                        angle,
                        R(1, 5),
                        R(1, 2),
                        gamma,
                        source.HEAVY_MASS2,
                    )
                    qs = (s.zeros(4, 1), -rr[0], -rr[0] - rr[1])
                    effective = (gamma * kk, qs[1], qs[2])
                    bar = sum((w * q for w, q in zip(ws, qs)), s.zeros(4, 1))
                    M = sum(w * m for w, m in zip(ws, mm))
                    U = sum(
                        ws[i]
                        * ws[j]
                        * dot(effective[i] - effective[j], effective[i] - effective[j])
                        for i in range(3)
                        for j in range(i + 1, 3)
                    )
                    den = M - U
                    num = -2 * ws[0] * (bar.T * eps * bar)[0]
                    jac = (1 - z) * (1 - s.I * R(3, 5))
                    checks[tag + f"_frozen_line{split}_denominator"] = s.factor(
                        actual["denominator"] - den
                    )
                    checks[tag + f"_frozen_line{split}_weighted_TT_density"] = s.factor(
                        actual["density"] * den**2 - jac * num
                    )
                    kernel_count += 1
    gates.update(
        {
            "three_original_massive_states_retained": True,
            "full_vector_orders_two_four_six_checked": literal_count == 54,
            "all54_original_mass_ordered_frozen_line_kernels_checked": kernel_count
            == 54,
            "known_curvature_response_nonzero_on_original_states": nonzero > 0,
            "off_shell_scaling_not_a_change_in_physical_light_mass": True,
            "uniform_and_generic_proof_not_inferred_from_calibrations": True,
        }
    )
    return {
        "checks": checks,
        "gates": gates,
        "whole_independent_counts": {
            "original_states": 3,
            "polarizations": 6,
            "vector_order_checks": literal_count,
            "frozen_kernel_checks": kernel_count,
            "nonzero_curvature_words": nonzero,
        },
        "whole_calibration_boundary": "Independent4-vector cyclic routing, gamma-polynomial integration and literal256-component curvature contractions on three original massive recoil states and both polarizations, with all three pairings. Separately all three frozen S342 triangle line kernels at the original heavy mass are checked after a stated common off-shell1/16 momentum rescaling into the analytic region. That rescaling tests local jets and does not redefine the source mass or approximate the above-threshold physical amplitude.",
    }
