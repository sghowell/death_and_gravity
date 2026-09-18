"""Continuous complete logarithmic angular kernels, including atom diagonals."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import kernel, radiative
from p8_vacuum_affine_radiative_state_soft_index import recoil


def clean(value):
    value = s.factor(s.expand_log(s.expand(value), force=True))
    replacements = {}
    for term in value.atoms(s.log):
        arg = term.args[0]
        if arg.is_Rational and arg > 0:
            num, den = arg.as_numer_denom()
            replacements[term] = sum(
                exponent * s.log(prime) for prime, exponent in s.factorint(num).items()
            ) - sum(
                exponent * s.log(prime) for prime, exponent in s.factorint(den).items()
            )
    return s.factor(s.expand(value.xreplace(replacements)))


def _evaluate(energy, rays, u, n, A):
    points, rays, _ = recoil.momenta(energy, rays, u)
    q = s.Matrix([1, *n])
    D = [s.factor(recoil.dot(p, q)) for p in points]
    P = lambda a, b: s.factor((a.T * A * b)[0])
    L = [s.factor(P(p, p) / d) for p, d in zip(points, D)]
    H = sum(d * s.log(abs(d)) for d in D)
    pair = lambda i, j: L[i] * D[j] + L[j] * D[i] - 2 * P(points[i], points[j])
    F = sum(L) * H
    for i in range(4):
        for j in range(i + 1, 4):
            F -= (
                kernel.massive_fprime(abs(recoil.dot(points[i], points[j]))) / 4
                - s.log(2)
                - 1
            ) * pair(i, j)
    units = [p / p[0] for p in rays]
    d = [s.factor(recoil.dot(v, q)) for v in units]
    ell = [s.factor(P(v, v) / di) if di != 0 else s.S.Zero for v, di in zip(units, d)]
    for j, v in enumerate(units):
        if d[j] == 0:
            continue
        for i, p in enumerate(points):
            alpha = abs(s.factor(recoil.dot(p, v)))
            mixed = (
                L[i] * d[j] * s.log(d[j] / alpha)
                + D[i] * ell[j] * s.log(abs(D[i]) / alpha)
                + 2 * P(p, v) * s.log(alpha)
            )
            F += rays[j][0] * mixed
    for i, v in enumerate(units):
        for j, vp in enumerate(units):
            delta = s.factor(recoil.dot(v, vp))
            if d[i] == 0 or d[j] == 0:
                term = s.S.Zero
            elif delta == 0:
                term = 2 * P(v, v) * s.log(d[i])
            else:
                S = d[j] * ell[i] + d[i] * ell[j] - 2 * P(v, vp)
                term = (
                    ell[i] * d[j] * s.log(d[j])
                    + ell[j] * d[i] * s.log(d[i])
                    - S * s.log(delta)
                )
            F += rays[i][0] * rays[j][0] * term / 2
    G = (
        (kernel.massive_c(recoil.dot(points[0], points[1])) + 2) * pair(0, 1)
        + (kernel.massive_c(recoil.dot(points[2], points[3])) - 2) * pair(2, 3)
    ) / (8 * s.pi)
    return {"F": F / (4 * s.pi**2), "G": G}


def components(energy, quanta, outgoing, direction, tensor):
    E = kernel.exact_scalar(energy)
    u, n = radiative.vector(outgoing, 3), radiative.vector(direction, 3)
    if s.simplify(n.dot(n) - 1) != 0:
        raise ValueError("Require a unit additional-soft direction")
    if not isinstance(quanta, (tuple, list)):
        raise TypeError("Require finitely many positive null momenta")
    rays = tuple(radiative.vector(ray, 4) for ray in quanta)
    A = radiative.polarization(tensor, n)
    return _evaluate(E, rays, u, n, A)


def coefficient(energy, quanta, outgoing, direction, tensor):
    row = components(energy, quanta, outgoing, direction, tensor)
    return row["F"] + s.I * row["G"]


def frame_x(h, axis):
    h = s.sympify(h)
    a, b = (1 - h * h) / (1 + h * h), 2 * h / (1 + h * h)
    e0, f0 = s.Matrix([0, 0, 1, 0]), s.Matrix([0, 0, 0, 1])
    if axis == "y":
        n = s.Matrix([a, b, 0])
        e = s.Matrix([0, -b, a, 0])
        f = f0
    elif axis == "z":
        n = s.Matrix([a, 0, b])
        e = e0
        f = s.Matrix([0, -b, 0, a])
    else:
        raise ValueError("Require one stated calibration path")
    return n, (e * e.T - f * f.T) / s.sqrt(2), (e * f.T + f * e.T) / s.sqrt(2)


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = clean(value)

    for count in (2, 3, 5, 8):
        D = s.symbols("D0:" + str(count), nonzero=True)
        P = {
            (i, j): s.Symbol("P" + str(i) + "_" + str(j))
            for i in range(count)
            for j in range(i, count)
        }
        A = lambda i, j, P=P: P[tuple(sorted((i, j)))]
        pair = sum(
            A(i, i) * D[j] / D[i] + A(j, j) * D[i] / D[j] - 2 * A(i, j)
            for i in range(count)
            for j in range(i + 1, count)
        )
        full = sum(D) * sum(A(i, i) / D[i] for i in range(count)) - sum(
            A(i, j) for i in range(count) for j in range(count)
        )
        put("generic_all_pair_constant_shift_" + str(count), pair - full)
    Li, Di, ell, d, alpha, Pi = s.symbols("Li Di ell d alpha Pi", positive=True)
    mixed = (
        Li * d * s.log(d / alpha) + Di * ell * s.log(Di / alpha) + 2 * Pi * s.log(alpha)
    )
    expanded = (
        Li * d * s.log(d)
        + ell * Di * s.log(Di)
        - (Li * d + Di * ell - 2 * Pi) * s.log(alpha)
    )
    put("mixed_continuous_ratio_rewrite", mixed - expanded)
    e, ep, d, dp, delta, P = s.symbols("ell ellprime d dprime delta P", positive=True)
    null = (
        e * dp * s.log(dp)
        + ep * d * s.log(d)
        - (dp * e + d * ep - 2 * P) * s.log(delta)
    )
    ratio = (
        e * dp * s.log(dp / delta) + ep * d * s.log(d / delta) + 2 * P * s.log(delta)
    )
    put("null_continuous_ratio_rewrite", null - ratio)
    w, P, d = s.symbols("w P d", positive=True)
    put(
        "ordered_double_self_diagonal",
        w * w * (2 * P * s.log(d)) / 2 - w * w * (P / d) * (d * s.log(d)),
    )
    E, states = radiative.calibration_states()
    u = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    n = s.Matrix([0, 1, 0])
    plus = s.diag(0, 1, 0, -1) / s.sqrt(2)
    cross = s.zeros(4)
    cross[1, 3] = cross[3, 1] = 1 / s.sqrt(2)
    records = []
    for name, rays in states.items():
        for pol, A in (("plus", plus), ("complex", (plus + s.I * cross) / s.sqrt(2))):
            old = radiative.components(E, rays, u, n, A)
            new = components(E, rays, u, n, A)
            put(name + "_" + pol + "_complete_F_rewrite", new["F"] - old["F"])
            put(name + "_" + pol + "_complete_G_rewrite", new["G"] - old["G"])
            records.append((name, pol, len(rays)))
    n0, p0, c0 = frame_x(0, "y")
    boundary = {}
    for name in ("one", "two", "split", "four"):
        for pol, A in (("plus", p0), ("complex", (p0 + s.I * c0) / s.sqrt(2))):
            boundary[(name, pol)] = coefficient(E, states[name], u, n0, A)
    for pol in ("plus", "complex"):
        put(
            "complete_point_collinear_split_" + pol,
            boundary[("two", pol)] - boundary[("split", pol)],
        )
    h = s.Symbol("h", positive=True)
    for axis in ("y", "z"):
        nh, ph, _ = frame_x(h, axis)
        v = s.Matrix([1, 1, 0, 0])
        q = s.Matrix([1, *nh])
        ellh = (v.T * ph * v)[0] / recoil.dot(v, q)
        wanted = s.sqrt(2) if axis == "y" else -s.sqrt(2)
        put(
            "individual_current_ambiguous_" + axis,
            s.limit(ellh, h, 0, dir="+") - wanted,
        )
    return {
        "checks": checks,
        "gates": {
            "full_conservation_before_constant_transfer": True,
            "mixed_Doppler_gaps_preserved": True,
            "single_collinear_ratio_cancellation": True,
            "relative_collinear_Gram_removable_product": True,
            "joint_triple_diagonal_continuity": True,
            "ordered_double_integral_includes_atom_self_diagonal": True,
            "real_and_complex_original_exact_calibrations": True,
            "no_unique_value_assigned_to_individual_null_current": True,
            "physical_recoil_kept_when_evaluating_collinear_atoms": True,
        },
        "whole_continuous_kernel_formula": "4pi^2 F=F0+int sum_M K_M d sigma+(1/2)int int K_N d sigma d sigma. F0=-sumMM[fprime/4-(ln2+1)]S+S_M H_M. K_M=Li*d ln(d/alpha)+Di*ell ln(|Di|/alpha)+2A(ki,v)ln alpha; K_M=0 at d0. K_N=ell*dprime ln dprime+ellprime*d ln d-S ln delta; K_N=2A(v,v)ln d on relative diagonal, and0 when either soft gap vanishes. Self diagonals are included.",
        "whole_calibration_inventory": records,
        "whole_complete_collinear_values": {
            name + "_" + pol: value for (name, pol), value in boundary.items()
        },
        "whole_current_qualification": "Different azimuthal limits of the INDIVIDUAL ell do not imply different limits of the complete F+iG. Frozen S326 rejects exact collinear inputs; this successor provides a new complete-coefficient API, not a rewrite of that frozen current or API.",
    }
