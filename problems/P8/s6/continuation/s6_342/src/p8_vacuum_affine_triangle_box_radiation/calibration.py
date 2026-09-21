"""Original labelled-source, coincidence and every routed contour-density control."""

import itertools
from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree
from p8_vacuum_affine_minimal_gravity_radiation import vertices as actual

from . import radiation


@cache
def data():
    checks = {}
    gates = {}
    equal = wrong = nonzero = 0
    eta = s.diag(1, -1, -1, -1)
    n, g, C = germs.MASS2, germs.G, germs.CONTACT

    def dot(p, q):
        return (p.T * eta * q)[0]

    def ep(p, q, E):
        return (p.T * E * q)[0]

    def A(v):
        return C + g * g / (n - v)

    fixtures = [actual.sample(i) for i in range(3)]
    E, Ep = s.Rational(5, 4), s.Rational(461, 380)
    fixtures.append(
        actual.recoil(
            E, E - Ep * Ep / E, (1, 0, 0), (0, s.Rational(3, 5), s.Rational(4, 5))
        )
    )
    for index, (ps, k, born) in enumerate(fixtures):
        # Toy loop kernels Cbar=1,Dbar=1 test physical labelled normalization
        # and outer-branch placement. They are NOT evaluated actual loop values.
        for polidx, (eps, norm) in enumerate(
            tensor.transverse_polarizations(tensor.frame(k), eta)
        ):
            Js = [ep(p, p, eps) / dot(p, k) for p in ps]
            tri = 0
            box = 0
            outer = 0
            for a, b, c, d in itertools.permutations(range(4)):
                P = ps[a] + ps[b]
                v = dot(P, P)
                u = dot(P + k, P + k)
                for leg in range(4):
                    q = [p + k if j == leg else p for j, p in enumerate(ps)]
                    PQ = q[a] + q[b]
                    tri += g * g * Js[leg] * A(dot(PQ, PQ)) / 4
                    box += g**4 * Js[leg] / 4
                outer += -g * g * 2 * g * g * ep(P, P, eps) / ((n - v) * (n - u)) / 4
            tri += outer
            reftri = tree.whole_tensor(
                ps, k, mass=1, heavy=n, cubic=g, contact=3 * C, kappa=1
            )
            reft = (
                s.Rational(2)
                * g
                * g
                * sum(eps[i, j] * reftri[i, j] for i in range(4) for j in range(4))
            )
            refbox = tree.whole_tensor(
                ps, k, mass=1, heavy=n, cubic=0, contact=6 * g**4, kappa=1
            )
            refb = sum(eps[i, j] * refbox[i, j] for i in range(4) for j in range(4))
            checks[
                f"original{index}_pol{polidx}_complete_triangle_constant_kernel_vs26graph"
            ] = s.factor(tri - reft)
            checks[
                f"original{index}_pol{polidx}_complete_box_constant_kernel_vs_contact"
            ] = s.factor(box - refb)
            wrong += int(s.factor(tri - 2 * outer - reft) != 0)
            nonzero += int(tri != 0 and box != 0)
            # The actual triangle outer branch takes C(u), not C(v).
            # C_test(v)=v is a covariant two-endpoint kernel control.
            for left, right in actual.PARTS:
                P = ps[left[0]] + ps[left[1]]
                v, u = dot(P, P), dot(P + k, P + k)
                da = g * g / ((n - v) * (n - u))
                # A(v)*Cprime + Aprime*C(u) is DD[A(v)*v].
                dd = (
                    (A(u) * u - A(v) * v) / (u - v)
                    if u != v
                    else A(v) + v * g * g / (n - v) ** 2
                )
                checks[
                    f"original{index}_pol{polidx}_part{left}_ordered_outer_branch"
                ] = s.factor(-2 * ep(P, P, eps) * (A(v) + da * u - dd))
                if u == v:
                    equal += 1
                if u != v and ep(P, P, eps) != 0:
                    assert s.factor(A(v) + da * v - dd) != 0

    gates["physical_coincident_channels_present"] = equal == 4
    gates["wrong_outer_sign_fails_original_physical_cases"] = wrong == 8
    gates["nonzero_original_complete_TT_controls"] = nonzero == 8
    counts = {
        "routes": 0,
        "physical_TT": 0,
        "nonzero_TT": 0,
        "timelike": 0,
        "crossed": 0,
    }
    den_residual = s.S.Zero
    TT_residual = s.S.Zero
    eta = radiation.ETA
    n = germs.MASS2
    gam = s.Rational(2, 3)
    for index, (ps, k, born) in enumerate(fixtures):
        pols = tensor.transverse_polarizations(tensor.frame(k), eta)
        for a, b, c, d in itertools.permutations(range(4)):
            for vertices, masses in (
                ([ps[a] + ps[b], ps[c], ps[d]], [1, 1, n]),
                ([ps[a], ps[b], ps[c], ps[d]], [1, n, 1, n]),
            ):
                N = len(vertices)
                for split in range(N):
                    for polidx, (eps, norm) in enumerate(pols):
                        data = radiation.line_density(
                            vertices,
                            masses,
                            k,
                            eps,
                            split,
                            s.Rational(1, 5),
                            s.Rational(1, 7),
                            s.Rational(2, 5),
                            gam,
                            n,
                        )
                        weights = data["weights"]
                        eff = data["effective_route"]
                        if polidx == 0:
                            bar_eff = sum(
                                (u * q for u, q in zip(weights, eff)), s.zeros(4, 1)
                            )
                            direct = sum(
                                u * (mass - radiation.dot(q, q))
                                for u, mass, q in zip(weights, data["masses"], eff)
                            ) + radiation.dot(bar_eff, bar_eff)
                            residual = s.cancel(s.expand(direct - data["denominator"]))
                            assert residual == 0, (index, (a, b, c, d), N, split)
                            den_residual += residual * s.conjugate(residual)
                            assert s.im(s.expand(data["denominator"])) < 0
                            counts["routes"] += 1
                            counts[
                                "timelike" if data["light_invariant"] > 4 else "crossed"
                            ] += 1
                        R = data["shift_vector"]
                        yy = weights[0] * gam
                        literal = tree.tensor(
                            R - yy * k, R + (1 - yy) * k, data["masses"][0], eta
                        )
                        contracted = sum(
                            eps[i, j] * literal[i, j]
                            for i in range(4)
                            for j in range(4)
                        )
                        residual = s.cancel(
                            s.expand(contracted - 2 * (R.T * eps * R)[0])
                        )
                        assert residual == 0, (index, (a, b, c, d), N, split, polidx)
                        TT_residual += residual * s.conjugate(residual)
                        counts["physical_TT"] += 1
                        counts["nonzero_TT"] += int(
                            s.cancel(s.expand(data["numerator"])) != 0
                        )
    assert den_residual == TT_residual == 0
    assert counts["routes"] == 672 and counts["physical_TT"] == 1344
    assert counts["nonzero_TT"] > 0 and counts["timelike"] > 0 and counts["crossed"] > 0

    checks["all672_full_routed_denominators_sum_absolute_squares"] = den_residual
    checks["all1344_literal_TT_vertices_sum_absolute_squares"] = TT_residual
    gates["all_original_line_density_counts"] = counts == {
        "routes": 672,
        "physical_TT": 1344,
        "nonzero_TT": 1284,
        "timelike": 224,
        "crossed": 448,
    }
    return {
        "checks": {key: s.cancel(s.expand(value)) for key, value in checks.items()},
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_coincident_polarized_channels": equal,
        "whole_wrong_sign_rejections": wrong,
        "whole_nonzero_original_TT_controls": nonzero,
        "whole_literal_density_counts": counts,
        "whole_control_boundary": "Four original recoil states, original n,g,C and both real physical polarizations. Constant diagnostic loop kernels compare the complete labelled triangle/outer-heavy and box amplitudes with the independent original26graph tree engine; they are not claimed to be actual loop values. The v test kernel distinguishes the ordered outer branch, including four coincident polarized channels. Every actual triangle/box label and split line is then checked against the full original-mass Symanzik denominator and literal scalar stress at a complex contour point:672 routes and1344 TT contractions,1284 nonzero. Written operator identities and continuous contour bounds, not these finite samples, establish the full result.",
    }
